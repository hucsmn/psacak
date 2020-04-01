use std::marker::PhantomData;

use crossbeam::channel::{self, Receiver, Sender};
use scoped_threadpool::{Pool, Scope};

/// Pipeline infrastructure for paralleled induce sorting.
pub struct Pipeline {
    pool: Pool,
}

impl Pipeline {
    pub fn new() -> Self {
        let pool = Pool::new(2);
        Pipeline { pool }
    }

    /// Start the induce pipeline.
    pub fn begin<'scope, S, T, PF, FL, BODY>(&mut self, prefetch: PF, flush: FL, body: BODY)
    where
        S: Send + 'scope,
        T: Send + 'scope,
        PF: Fn(S) -> S + 'scope + Send,
        FL: Fn(T) -> T + 'scope + Send,
        BODY: FnOnce(Worker<'scope, S>, Worker<'scope, T>) + 'scope + Send,
    {
        self.pool.scoped(|scope| {
            let pf = Worker::new(scope, prefetch);
            let fl = Worker::new(scope, flush);
            body(pf, fl);
        });
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
