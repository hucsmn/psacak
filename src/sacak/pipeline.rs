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
        INDUCE: FnOnce(OuterWorker<'scope, S>, OuterWorker<'scope, T>) + 'scope + Send,
    {
        if self.pool.is_none() {
            self.pool = Some(Pool::new(2));
        }
        if let Some(ref mut pool) = self.pool {
            pool.scoped(|scope| {
                let fetch_worker = OuterWorker::new(scope, fetch);
                let flush_worker = OuterWorker::new(scope, flush);
                induce(fetch_worker, flush_worker);
            });
        } else {
            unreachable!();
        }
    }
}

/// Worker thread controller for outer level SACA-K.
pub struct OuterWorker<'scope, S: Send + 'scope> {
    input: Sender<S>,
    output: Receiver<S>,
    _marker: PhantomData<&'scope S>,
}

impl<'scope, S: Send + 'scope> OuterWorker<'scope, S> {
    fn new<'pool, F>(scope: &Scope<'pool, 'scope>, action: F) -> Self
    where
        F: Fn(S) -> S + 'scope + Send,
    {
        let (insend, inrecv) = channel::bounded(0);
        let (outsend, outrecv) = channel::bounded(0);
        scope.execute(move || {
            // channel send/recv and thread create/join would make synchronization.
            while let Ok(mut state) = inrecv.recv() {
                state = action(state);
                outsend.send(state).unwrap();
            }
        });
        OuterWorker {
            input: insend,
            output: outrecv,
            _marker: PhantomData,
        }
    }

    /// Start a new task.
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
    use super::*;

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
}
