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
    pub fn ready(&self, state: S) {
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

    impl Buf {
        pub fn new() -> Self {
            Buf(std::usize::MAX, PhantomData)
        }
    }
    #[quickcheck]
    fn quickcheck_pipeline_outer(n: usize) {
        if n == 0 {
            return;
        }
        let mut pipeline = Pipeline::new();

        // simulate the time sequence of outer level induce pipeline.
        pipeline.outer_induce(
            |buf| buf,
            |buf| buf,
            |fetch, flush| {
                let mut rbuf = Buf::new();
                let mut rbuf_fetch = Buf::new();
                let mut wbuf = Buf::new();
                let mut wbuf_flush;

                rbuf_fetch.0 = 0;
                fetch.ready(rbuf_fetch);

                for i in 0..n {
                    rbuf_fetch = fetch.wait();
                    swap(&mut rbuf, &mut rbuf_fetch);
                    rbuf_fetch.0 = i + 1;
                    fetch.ready(rbuf_fetch);

                    assert_eq!(rbuf.0, i);

                    if i == 0 {
                        wbuf_flush = Buf::new();
                    } else {
                        wbuf_flush = flush.wait();
                        assert_eq!(wbuf_flush.0, i + 1);
                    }
                    swap(&mut wbuf, &mut wbuf_flush);
                    wbuf_flush.0 = i + 2;
                    flush.ready(wbuf_flush);
                }

                fetch.wait();
                flush.wait();
            },
        );
    }

    #[quickcheck]
    fn quickcheck_pipeline_inner(n: usize) {
        if n == 0 {
            return;
        }
        let mut pipeline = Pipeline::new();

        // simulate the time sequence of inner level induce pipeline.
        let mut rbuf = Buf::new();
        let mut rbuf_fetch = Buf::new();
        let mut wbuf = Buf::new();
        let mut wbuf_flush = Buf::new();
        pipeline.inner_induce(
            |buf| buf,
            |buf| buf,
            |fetch, flush| {
                flush.ready(wbuf_flush);
                rbuf_fetch.0 = 0;
                fetch.ready(rbuf_fetch);
                wbuf_flush = flush.wait();
                wbuf_flush.0 = 1;
                flush.ready(wbuf_flush);

                for i in 0..n {
                    rbuf_fetch = fetch.wait();
                    swap(&mut rbuf, &mut rbuf_fetch);
                    rbuf_fetch.0 = i + 1;
                    fetch.ready(rbuf_fetch);

                    assert_eq!(rbuf.0, i);

                    wbuf_flush = flush.wait();
                    swap(&mut wbuf, &mut wbuf_flush);
                    wbuf_flush.0 = i + 2;
                    flush.ready(wbuf_flush);

                    assert_eq!(wbuf.0, i + 1);
                }

                fetch.ready(fetch.wait());
                rbuf_fetch = fetch.wait();
                wbuf_flush = flush.wait();
            },
        );
    }

    #[quickcheck]
    fn quickcheck_pipeline_reuse(n: usize) {
        let mut pipeline = Pipeline::new();
        pipeline.inner_induce(
            |x| x,
            |x| x,
            |worker_a, worker_b| {
                worker_a.ready(0);
                worker_b.ready(0);
                assert_eq!(worker_a.wait(), 0);
                assert_eq!(worker_b.wait(), 0);
            },
        );
        for i in 1..=n {
            pipeline.outer_induce(
                |x| x,
                |x| x,
                |worker_a, worker_b| {
                    worker_a.ready(i);
                    worker_b.ready(i);
                    assert_eq!(worker_a.wait(), i);
                    assert_eq!(worker_b.wait(), i);
                },
            );
        }
    }
}
