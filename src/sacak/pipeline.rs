use std::marker::PhantomData;

use crossbeam::channel::{self, Receiver, Sender};
use scoped_threadpool::{Pool, Scope};

/// Pipeline infrastructure for induce sorting in parallel, which helps to improve data throughput.
pub struct Pipeline {
    pool: Option<Pool>,
}

impl Pipeline {
    pub fn new() -> Self {
        Pipeline { pool: None }
    }

    /// Start the induce pipeline for outer level SACA-K.
    pub fn outer_induce<'scope, S, T, FETCH, FLUSH, INDUCE>(&mut self, fetch: FETCH, flush: FLUSH, induce: INDUCE)
    where
        S: Send + 'scope,
        T: Send + 'scope,
        FETCH: Fn(S) -> S + 'scope + Send,
        FLUSH: Fn(T) -> T + 'scope + Send,
        INDUCE: FnOnce(Worker<'scope, S>, Worker<'scope, T>) + 'scope + Send,
    {
        if self.pool.is_none() {
            self.pool = Some(Pool::new(2));
        }
        if let Some(ref mut pool) = self.pool {
            pool.scoped(|scope| {
                let fetch_worker = Self::new_free_worker(scope, fetch);
                let flush_worker = Self::new_free_worker(scope, flush);
                induce(fetch_worker, flush_worker);
            });
        } else {
            unreachable!();
        }
    }

    /// Start the induce pipeline for inner level SACA-K.
    pub fn inner_induce<'scope, S, T, FETCH, FLUSH, INDUCE>(&mut self, fetch: FETCH, flush: FLUSH, induce: INDUCE)
    where
        S: Send + 'scope,
        T: Send + 'scope,
        FETCH: Fn(S) -> S + 'scope + Send,
        FLUSH: Fn(T) -> T + 'scope + Send,
        INDUCE: FnOnce(Worker<'scope, S>, Worker<'scope, T>) + 'scope + Send,
    {
        if self.pool.is_none() {
            self.pool = Some(Pool::new(2));
        }
        if let Some(ref mut pool) = self.pool {
            pool.scoped(|scope| {
                let (fetch_worker, flush_worker) = Self::new_worker_pair(scope, fetch, flush);
                induce(fetch_worker, flush_worker);
            });
        } else {
            unreachable!();
        }
    }

    /// Create new worker that iterates and alters the state using given action.
    fn new_free_worker<'pool, 'scope, S, F>(scope: &Scope<'pool, 'scope>, action: F) -> Worker<'scope, S>
    where
        S: Send + 'scope,
        F: Fn(S) -> S + 'scope + Send,
    {
        let (insend, inrecv) = channel::bounded(0);
        let (outsend, outrecv) = channel::bounded(0);

        scope.execute(move || {
            while let Ok(mut state) = inrecv.recv() {
                state = action(state);
                outsend.send(state).unwrap();
            }
        });

        Worker {
            input: insend,
            output: outrecv,
            _marker: PhantomData,
        }
    }

    /// Create new pair of workers that one happens before another.
    fn new_worker_pair<'pool, 'scope, S, T, FETCH, FLUSH>(
        scope: &Scope<'pool, 'scope>,
        fetch: FETCH,
        flush: FLUSH,
    ) -> (Worker<'scope, S>, Worker<'scope, T>)
    where
        S: Send + 'scope,
        T: Send + 'scope,
        FETCH: Fn(S) -> S + 'scope + Send,
        FLUSH: Fn(T) -> T + 'scope + Send,
    {
        let (flush_insend, flush_inrecv) = channel::bounded(0);
        let (flush_outsend, flush_outrecv) = channel::bounded(1);

        let (fetch_insend, fetch_inrecv) = channel::bounded(1);
        let (fetch_outsend, fetch_outrecv) = channel::bounded(0);

        let (blocksend, blockrecv) = channel::bounded(0);

        scope.execute(move || {
            while let Ok(mut state) = flush_inrecv.recv() {
                state = flush(state);
                flush_outsend.send(state).unwrap();
                blocksend.send(()).unwrap();
            }
        });

        scope.execute(move || {
            while let Ok(mut state) = fetch_inrecv.recv() {
                blockrecv.recv().unwrap();
                state = fetch(state);
                fetch_outsend.send(state).unwrap();
            }
        });

        let fetch_worker = Worker {
            input: fetch_insend,
            output: fetch_outrecv,
            _marker: PhantomData,
        };
        let flush_worker = Worker {
            input: flush_insend,
            output: flush_outrecv,
            _marker: PhantomData,
        };
        (fetch_worker, flush_worker)
    }
}

/// Worker thread controller.
pub struct Worker<'scope, S: Send + 'scope> {
    input: Sender<S>,
    output: Receiver<S>,
    _marker: PhantomData<&'scope S>,
}

impl<'scope, S: Send + 'scope> Worker<'scope, S> {
    /// Make a new task ready.
    pub fn start(&self, state: S) {
        self.input.send(state).unwrap();
    }

    /// Wait for the previous task.
    pub fn wait(&self) -> S {
        self.output.recv().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use std::marker::PhantomData;
    use std::mem::swap;

    use super::*;

    struct Buf(pub usize, PhantomData<Vec<u32>>);

    #[quickcheck]
    fn quickcheck_pipeline_outer(data: Vec<u8>) -> bool {
        let mut success = true;
        let mut pipeline = Pipeline::new();
        pipeline.outer_induce(
            |x| x ^ 0b11110000,
            |x| x ^ 0b00001111,
            |worker_a, worker_b| {
                for x in data.into_iter() {
                    worker_a.start(x);
                    worker_b.start(x);
                    if worker_a.wait() ^ worker_b.wait() != 0b11111111 {
                        success = false;
                        break;
                    }
                }
            },
        );
        success
    }

    #[quickcheck]
    fn quickcheck_pipeline_inner(n: usize) -> bool {
        let mut success = true;
        let mut pipeline = Pipeline::new();
        pipeline.inner_induce(
            |buf| buf,
            |buf| buf,
            |fetch, flush| {
                let mut rbuf0 = Buf(std::usize::MAX, PhantomData);
                let mut rbuf1 = Buf(std::usize::MAX, PhantomData);
                let mut wbuf0 = Buf(std::usize::MAX, PhantomData);
                let mut wbuf1 = Buf(std::usize::MAX, PhantomData);

                flush.start(wbuf1);

                rbuf1.0 = 0;
                fetch.start(rbuf1);

                wbuf1 = flush.wait();
                wbuf1.0 = 1;
                flush.start(wbuf1);

                for i in 0..n {
                    rbuf1 = fetch.wait();
                    swap(&mut rbuf0, &mut rbuf1);
                    rbuf1.0 = i + 1;
                    fetch.start(rbuf1);

                    if rbuf0.0 != i {
                        success = false;
                    }

                    wbuf1 = flush.wait();
                    swap(&mut wbuf0, &mut wbuf1);
                    wbuf1.0 = i + 2;
                    flush.start(wbuf1);

                    if wbuf0.0 != i + 1 {
                        success = false;
                    }

                    if !success {
                        break;
                    }
                }

                rbuf1 = fetch.wait();
                fetch.start(rbuf1);

                rbuf1 = fetch.wait();
                wbuf1 = flush.wait();
            },
        );
        success
    }
}
