use std::mem::swap;
use std::ops::Range;

use crossbeam::channel;
use scoped_threadpool::{Pool, Scope};

use super::common::ceil_divide;

/// The read buffer from text.
pub trait ReadBuffer: Send {
    fn reset(&mut self, block_size: usize);
}

/// The write buffer into workspace.
pub trait WriteBuffer: Send {
    fn reset(&mut self);
}

/// The induce pipeline infrastructure.
pub struct Pipeline {
    pool: Option<Pool>,
}

impl Pipeline {
    /// Create an induce pipeline.
    #[inline(always)]
    pub fn new() -> Self {
        Pipeline { pool: None }
    }

    /// Start induce pipeline for the outer-level pSACAK.
    #[inline]
    pub fn induce_outer<'scope, RBUF, WBUF, FETCH, FLUSH, INDUCE>(
        &mut self,
        backward: bool,
        length: usize,
        block_size: usize,
        buffers: (RBUF, RBUF, WBUF, WBUF),
        fetch: FETCH,
        flush: FLUSH,
        mut induce: INDUCE,
    ) -> (RBUF, RBUF, WBUF, WBUF)
    where
        RBUF: ReadBuffer + 'scope,
        WBUF: WriteBuffer + 'scope,
        FETCH: Fn((Range<usize>, RBUF)) -> RBUF + Send + 'scope,
        FLUSH: Fn(WBUF) -> WBUF + Send + 'scope,
        INDUCE: FnMut(&mut InduceContext<RBUF, WBUF>) + 'scope,
    {
        let (rbuf0, rbuf1, wbuf0, wbuf1) = buffers;
        self.pool().scoped(|scope| {
            // setup the pipeiline.
            let fetch_worker = Self::spawn_independent(scope, fetch);
            let flush_worker = Self::spawn_independent(scope, flush);
            let mut ctx = InduceContext::new(backward, length, block_size, rbuf0, wbuf0, fetch_worker, flush_worker);

            // invoke the induce pipeline.
            ctx.pre_induce(rbuf1, wbuf1);
            for _ in 0..ctx.blocks() {
                ctx.sync_fetch();
                induce(&mut ctx);
                ctx.sync_flush();
            }

            // reuse buffers.
            let (rbuf1, wbuf1) = ctx.post_induce();
            let (rbuf0, wbuf0) = ctx.into_buffers();
            (rbuf0, rbuf1, wbuf0, wbuf1)
        })
    }

    /// Start induce pipeline for the inner-level pSACAK.
    #[inline]
    pub fn induce_inner<'scope, RBUF, WBUF, FETCH, FLUSH, INDUCE>(
        &mut self,
        backward: bool,
        length: usize,
        block_size: usize,
        buffers: (RBUF, RBUF, WBUF, WBUF),
        fetch: FETCH,
        flush: FLUSH,
        mut induce: INDUCE,
    ) -> (RBUF, RBUF, WBUF, WBUF)
    where
        RBUF: ReadBuffer + 'scope,
        WBUF: WriteBuffer + 'scope,
        FETCH: Fn((Range<usize>, RBUF)) -> RBUF + Send + 'scope,
        FLUSH: Fn(WBUF) -> WBUF + Send + 'scope,
        INDUCE: FnMut(&mut InduceContext<RBUF, WBUF>) + 'scope,
    {
        let (rbuf0, rbuf1, wbuf0, wbuf1) = buffers;
        self.pool().scoped(|scope| {
            // setup the pipeiline.
            // in the inner-level pSACAK, flush must happen before fetch to avoid data race.
            let (flush_worker, fetch_worker) = Self::spawn_sequential(scope, flush, fetch);
            let mut ctx = InduceContext::new(backward, length, block_size, rbuf0, wbuf0, fetch_worker, flush_worker);

            // invoke the induce pipeline.
            ctx.pre_induce(rbuf1, wbuf1);
            for _ in 0..ctx.blocks() {
                ctx.sync_fetch();
                induce(&mut ctx);
                ctx.sync_flush();
            }

            // reuse buffers.
            let (rbuf1, wbuf1) = ctx.post_induce();
            let (rbuf0, wbuf0) = ctx.into_buffers();
            (rbuf0, rbuf1, wbuf0, wbuf1)
        })
    }

    #[inline(always)]
    fn pool(&mut self) -> &mut Pool {
        self.pool.get_or_insert_with(|| Pool::new(2))
    }

    #[inline(always)]
    fn spawn_independent<'scope, S, T, ACTION>(scope: &Scope<'_, 'scope>, action: ACTION) -> Worker<S, T>
    where
        S: Send + 'scope,
        T: Send + 'scope,
        ACTION: Fn(S) -> T + Send + 'scope,
    {
        let (ins, inr) = channel::bounded(0);
        let (outs, outr) = channel::bounded(0);

        scope.execute(move || {
            while let Ok(input) = inr.recv() {
                outs.send(action(input)).unwrap();
            }
        });

        Worker::new(ins, outr)
    }

    #[inline(always)]
    fn spawn_sequential<'scope, S0, T0, S1, T1, BEFORE, AFTER>(
        scope: &Scope<'_, 'scope>,
        before: BEFORE,
        after: AFTER,
    ) -> (Worker<S0, T0>, Worker<S1, T1>)
    where
        S0: Send + 'scope,
        T0: Send + 'scope,
        S1: Send + 'scope,
        T1: Send + 'scope,
        BEFORE: Fn(S0) -> T0 + Send + 'scope,
        AFTER: Fn(S1) -> T1 + Send + 'scope,
    {
        let (ins0, inr0) = channel::bounded(1);
        let (outs0, outr0) = channel::bounded(1);
        let (ins1, inr1) = channel::bounded(1);
        let (outs1, outr1) = channel::bounded(0);

        scope.execute(move || {
            while let Ok(input0) = inr0.recv() {
                outs0.send(before(input0)).unwrap();
                let input1 = inr1.recv().unwrap();
                outs1.send(after(input1)).unwrap();
            }
        });

        (Worker::new(ins0, outr0), Worker::new(ins1, outr1))
    }
}

/// Handle for a worker routine.
struct Worker<S: Send, T: Send> {
    input: channel::Sender<S>,
    output: channel::Receiver<T>,
}

impl<S: Send, T: Send> Worker<S, T> {
    #[inline(always)]
    pub fn new(input: channel::Sender<S>, output: channel::Receiver<T>) -> Self {
        Worker { input, output }
    }

    /// Ready to start the next action.
    #[inline(always)]
    pub fn ready(&self, input: S) {
        self.input.send(input).unwrap();
    }

    /// Wait for the previous action.
    #[inline(always)]
    pub fn wait(&self) -> T {
        self.output.recv().unwrap()
    }
}

/// The induce context.
///
/// It assumes that fetch action should do nothing if the input range of text is empty,
/// and that flush action should do nothing if the input write buffer is reset.
pub struct InduceContext<RBUF: ReadBuffer, WBUF: WriteBuffer> {
    length: usize,
    block_size: usize,

    cur_start: usize,
    next_start: usize,
    cur_end: usize,
    next_end: usize,

    pub rbuf: RBUF,
    pub wbuf: WBUF,

    fetch: Worker<(Range<usize>, RBUF), RBUF>,
    flush: Worker<WBUF, WBUF>,

    backward: bool,
}

impl<RBUF: ReadBuffer, WBUF: WriteBuffer> InduceContext<RBUF, WBUF> {
    /// Create induce state.
    #[inline(always)]
    fn new(
        backward: bool,
        length: usize,
        block_size: usize,
        mut rbuf: RBUF,
        mut wbuf: WBUF,
        fetch: Worker<(Range<usize>, RBUF), RBUF>,
        flush: Worker<WBUF, WBUF>,
    ) -> Self {
        rbuf.reset(0);
        wbuf.reset();
        if !backward {
            InduceContext {
                length,
                block_size,
                cur_start: 0,
                next_start: 0,
                cur_end: 0,
                next_end: Ord::min(block_size, length),
                rbuf,
                wbuf,
                fetch,
                flush,
                backward: false,
            }
        } else {
            InduceContext {
                length,
                block_size,
                cur_start: length,
                next_start: length.saturating_sub(block_size),
                cur_end: length,
                next_end: length,
                rbuf,
                wbuf,
                fetch,
                flush,
                backward: true,
            }
        }
    }

    /// Take out buffers.
    #[inline(always)]
    fn into_buffers(self) -> (RBUF, WBUF) {
        (self.rbuf, self.wbuf)
    }

    /// Number of blocks.
    #[inline(always)]
    pub fn blocks(&self) -> usize {
        ceil_divide(self.length, self.block_size)
    }

    /// Test if position is in the current or the next block.
    #[inline(always)]
    pub fn contains(&self, i: usize) -> bool {
        if !self.backward {
            i >= self.cur_start && i < self.next_end
        } else {
            i >= self.next_start && i < self.cur_end
        }
    }

    /// Start of the current block.
    #[inline(always)]
    pub fn cur_start(&self) -> usize {
        self.cur_start
    }

    /// End of the current block.
    #[inline(always)]
    pub fn cur_end(&self) -> usize {
        self.cur_end
    }

    /// Start of the next block.
    #[inline(always)]
    pub fn next_start(&self) -> usize {
        self.next_start
    }

    /// End of the next block.
    #[inline(always)]
    pub fn next_end(&self) -> usize {
        self.next_end
    }

    /// Initial step that fetches the first block.
    #[inline]
    fn pre_induce(&mut self, mut rbuf: RBUF, mut wbuf: WBUF) {
        rbuf.reset(self.next_end - self.next_start);
        wbuf.reset();
        self.flush.ready(wbuf);
        self.fetch.ready((self.next_start..self.next_end, rbuf));
        wbuf = self.flush.wait();
        self.flush.ready(wbuf);
    }

    /// Final step that waits for the workers to be done.
    #[inline]
    fn post_induce(&mut self) -> (RBUF, WBUF) {
        let mut rbuf = self.fetch.wait();
        rbuf.reset(0);
        self.fetch.ready((self.next_start..self.next_end, rbuf));
        let wbuf = self.flush.wait();
        rbuf = self.fetch.wait();
        (rbuf, wbuf)
    }

    /// Get the read buffer then start the next fetch before inducing a block.
    #[inline]
    fn sync_fetch(&mut self) {
        // alter the block boundaries.
        self.cur_start = self.next_start;
        self.cur_end = self.next_end;
        if !self.backward {
            self.next_start = self.next_end;
            self.next_end = Ord::min(self.next_end.saturating_add(self.block_size), self.length);
        } else {
            self.next_end = self.next_start;
            self.next_start = self.next_start.saturating_sub(self.block_size);
        }

        // reload the fetch worker.
        let mut rbuf = self.fetch.wait();
        swap(&mut self.rbuf, &mut rbuf);
        rbuf.reset(self.next_end - self.next_start);
        self.fetch.ready((self.next_start..self.next_end, rbuf));
    }

    /// Start the next flush after inducing a block.
    #[inline]
    fn sync_flush(&mut self) {
        // reload the flush worker.
        let mut wbuf = self.flush.wait();
        swap(&mut self.wbuf, &mut wbuf);
        self.flush.ready(wbuf);
        self.wbuf.reset();
    }
}

#[cfg(test)]
mod tests {
    use std::marker::PhantomData;
    use std::ops::Range;

    use super::*;

    #[derive(Clone)]
    struct Buf(Range<usize>, PhantomData<Vec<u32>>);

    impl Buf {
        pub fn new(range: Range<usize>) -> Self {
            Buf(range, PhantomData)
        }

        pub fn get(&self) -> Range<usize> {
            self.0.clone()
        }

        pub fn set(&mut self, range: Range<usize>) {
            self.0 = range
        }
    }

    impl ReadBuffer for Buf {
        fn reset(&mut self, _: usize) {}
    }

    impl WriteBuffer for Buf {
        fn reset(&mut self) {}
    }

    #[quickcheck]
    fn quickcheck_pipeline_reuse(length: usize, block_size: usize, times: usize) {
        let times = Ord::min(times, 100);
        let (length, block_size, _) = get_checked_params(length, block_size).unwrap_or((1, 1, 1));

        let mut pipeline = Pipeline::new();
        let buffers = (Buf::new(0..0), Buf::new(0..0), Buf::new(0..0), Buf::new(0..0));

        pipeline.induce_outer(
            false,
            length,
            block_size,
            buffers.clone(),
            |(_, buf)| buf,
            |buf| buf,
            |_| {},
        );

        for _ in 0..times {
            pipeline.induce_inner(
                false,
                length,
                block_size,
                buffers.clone(),
                |(_, buf)| buf,
                |buf| buf,
                |_| {},
            );
        }
    }

    #[quickcheck]
    #[allow(unused_parens)]
    fn quickcheck_pipeline_outer(length: usize, block_size: usize) {
        let (length, block_size, last_block_size) = get_checked_params(length, block_size).unwrap_or((1, 1, 1));
        let mut pipeline = Pipeline::new();

        // simulate the first stage of outer-level induce.
        let mut i = 0;
        let mut p = 0..0;
        let mut pp = 0..0;
        let buffers = (Buf::new(0..0), Buf::new(0..0), Buf::new(0..0), Buf::new(0..0));
        pipeline.induce_outer(
            false,
            length,
            block_size,
            buffers,
            |(range, mut rbuf): (Range<usize>, Buf)| {
                rbuf.set(range);
                rbuf
            },
            |wbuf| wbuf,
            |ctx| {
                let start = Ord::min(i * block_size, length);
                let end = Ord::min((i + 1) * block_size, length);
                let cur_start = ctx.cur_start();
                let cur_end = ctx.cur_end();
                assert_eq!(cur_start..cur_end, start..end);
                assert_eq!(ctx.rbuf.get(), cur_start..cur_end);
                assert_eq!(ctx.wbuf.get(), pp);

                ctx.wbuf.set(cur_start..cur_end);
                pp = p.clone();
                p = cur_start..cur_end;
                i += 1;
            },
        );
        assert_ne!(p.end - p.start, 0);
        assert_eq!(p, length - last_block_size..length);

        // simulate the second stage of outer-level induce.
        let mut i = 0;
        let mut p = length..length;
        let mut pp = length..length;
        let buffers = (
            Buf::new(length..length),
            Buf::new(length..length),
            Buf::new(length..length),
            Buf::new(length..length),
        );
        pipeline.induce_outer(
            true,
            length,
            block_size,
            buffers,
            |(range, mut rbuf): (Range<usize>, Buf)| {
                rbuf.set(range);
                rbuf
            },
            |wbuf| wbuf,
            |ctx| {
                let start = length.saturating_sub((i + 1) * block_size);
                let end = length.saturating_sub(i * block_size);
                let cur_start = ctx.cur_start();
                let cur_end = ctx.cur_end();
                assert_eq!(cur_start..cur_end, start..end);
                assert_eq!(ctx.rbuf.get(), cur_start..cur_end);
                assert_eq!(ctx.wbuf.get(), pp);

                ctx.wbuf.set(cur_start..cur_end);
                pp = p.clone();
                p = cur_start..cur_end;
                i += 1;
            },
        );
        assert_ne!(p.end - p.start, 0);
        assert_eq!(p, 0..last_block_size);
    }

    #[quickcheck]
    #[allow(unused_parens)]
    fn quickcheck_pipeline_inner(length: usize, block_size: usize) {
        let (length, block_size, last_block_size) = get_checked_params(length, block_size).unwrap_or((1, 1, 1));
        let mut pipeline = Pipeline::new();

        // simulate the first stage of inner-level induce.
        let mut i = 0;
        let mut p = 0..0;
        let mut pp = 0..0;
        let buffers = (Buf::new(0..0), Buf::new(0..0), Buf::new(0..0), Buf::new(0..0));
        pipeline.induce_inner(
            false,
            length,
            block_size,
            buffers,
            |(range, mut rbuf): (Range<usize>, Buf)| {
                rbuf.set(range);
                rbuf
            },
            |wbuf| wbuf,
            |ctx| {
                let start = Ord::min(i * block_size, length);
                let end = Ord::min((i + 1) * block_size, length);
                let cur_start = ctx.cur_start();
                let cur_end = ctx.cur_end();
                assert_eq!(cur_start..cur_end, start..end);
                assert_eq!(ctx.rbuf.get(), cur_start..cur_end);
                assert_eq!(ctx.wbuf.get(), pp);

                ctx.wbuf.set(cur_start..cur_end);
                pp = p.clone();
                p = cur_start..cur_end;
                i += 1;
            },
        );
        assert_ne!(p.end - p.start, 0);
        assert_eq!(p, length - last_block_size..length);

        // simulate the second stage of inner-level induce.
        let mut i = 0;
        let mut p = length..length;
        let mut pp = length..length;
        let buffers = (
            Buf::new(length..length),
            Buf::new(length..length),
            Buf::new(length..length),
            Buf::new(length..length),
        );
        pipeline.induce_inner(
            true,
            length,
            block_size,
            buffers,
            |(range, mut rbuf): (Range<usize>, Buf)| {
                rbuf.set(range);
                rbuf
            },
            |wbuf| wbuf,
            |ctx| {
                let start = length.saturating_sub((i + 1) * block_size);
                let end = length.saturating_sub(i * block_size);
                let cur_start = ctx.cur_start();
                let cur_end = ctx.cur_end();
                assert_eq!(cur_start..cur_end, start..end);
                assert_eq!(ctx.rbuf.get(), cur_start..cur_end);
                assert_eq!(ctx.wbuf.get(), pp);

                ctx.wbuf.set(cur_start..cur_end);
                pp = p.clone();
                p = cur_start..cur_end;
                i += 1;
            },
        );
        assert_ne!(p.end - p.start, 0);
        assert_eq!(p, 0..last_block_size);
    }

    // helper functions.

    fn get_checked_params(length: usize, mut block_size: usize) -> Option<(usize, usize, usize)> {
        if length == 0 {
            return None;
        }
        block_size = if block_size % length != 0 {
            block_size % length
        } else {
            length
        };
        let last_block_size = if length % block_size != 0 {
            length % block_size
        } else {
            block_size
        };
        Some((length, block_size, last_block_size))
    }
}
