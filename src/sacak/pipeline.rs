use std::marker::PhantomData;

use crossbeam::channel::{self, Receiver, Sender};
use scoped_threadpool::{Pool, Scope};

/// Pipeline infrastructure for induce sorting in parallel.
pub struct Pipeline {
    pool: Option<Pool>,
}

impl Pipeline {
    pub fn new() -> Self {
        Pipeline { pool: None }
    }

    /// Start the induce pipeline.
    pub fn begin<'scope, S, T, FETCH, FLUSH, INDUCE>(&mut self, fetch: FETCH, flush: FLUSH, induce: INDUCE)
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
                let fetch_worker = Worker::new(scope, fetch);
                let flush_worker = Worker::new(scope, flush);
                induce(fetch_worker, flush_worker);
            });
        } else {
            unreachable!();
        }
    }
}

/// Worker thread controller.
pub struct Worker<'scope, S: Send + 'scope> {
    input: Sender<S>,
    output: Receiver<S>,
    _marker: PhantomData<&'scope S>,
}

impl<'scope, S: Send + 'scope> Worker<'scope, S> {
    fn new<'pool, F>(scope: &Scope<'pool, 'scope>, action: F) -> Self
    where
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

    /// Start a new task.
    pub fn start(&self, state: S) {
        self.input.send(state).unwrap();
    }

    /// Wait for the previous task.
    pub fn wait(&self) -> S {
        self.output.recv().unwrap()
    }
}
