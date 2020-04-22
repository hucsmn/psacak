use std::cell::UnsafeCell;
use std::ops::Range;

use rayon::prelude::*;

use super::common::*;
use super::naming::*;
use super::pipeline::*;
use super::types::*;

/// Empty symbol in workspace.
const EMPTY: u32 = 1 << 31;

/// Block size for induce sorting in parallel.
const BLOCK_SIZE: usize = 64 * 1024;

/// Threshold to enable induce sorting in parallel.
const THRESHOLD_PARALLEL_INDUCE: usize = 2 * BLOCK_SIZE;

/// The inner level SACA-K algorithm for unsigned integer strings.
///
/// Assumes that characters are correctly translated to corresponding bucket pointers.
#[inline]
pub fn sacak32(text: &[u32], suf: &mut [u32], pipeline: &mut Pipeline) {
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // induce sort lms-substrings.
    put_lmscharacters(text, suf);
    induce_sort(text, suf, pipeline, BLOCK_SIZE, true);

    // construct the subproblem, compute its suffix array, and get sorted lms-suffixes.
    let n = compact_left_range(suf, 1..EMPTY);
    let k = name_lmssubstrings(text, suf, n);
    if k < n {
        // need to solve the subproblem recursively.
        let (suf1, text1) = suf.split_at_mut(suf.len() - n);
        sacak32(text1, suf1, pipeline);
    } else {
        // the subproblem itself is the inversed suffix array.
        let (suf1, text1) = suf.split_at_mut(suf.len() - n);
        for i in 0..n {
            suf1[text1[i].as_index()] = i as u32;
        }
    }
    permutate_lmssuffixes(text, suf, n);

    // induce sort the suffix array from sorted lms-suffixes.
    put_lmssuffixes(text, suf, n);
    induce_sort(text, suf, pipeline, BLOCK_SIZE, false);
}

/// Put lms-characters to their corresponding bucket tails, in arbitary order.
#[inline]
fn put_lmscharacters(text: &[u32], suf: &mut [u32]) {
    reset_slice(suf, EMPTY);

    foreach_lmschars(text, |i, c| {
        let p = c as usize;
        if suf[p] < EMPTY {
            // right shift the borrowed chunk.
            let n = suf[p..].iter().take_while(|&&x| x < EMPTY).count();
            suf.copy_within(p..p + n, p + 1);
            suf[p] = EMPTY;
        }

        if suf[p] == EMPTY {
            if p > 0 && suf[p - 1] == EMPTY {
                suf[p] = encode_counter(1);
                suf[p - 1] = i as u32;
            } else {
                suf[p] = i as u32;
            }
        } else {
            let q = p - decode_counter(suf[p]);
            if q > 0 && suf[q - 1] == EMPTY {
                suf[p] -= 1;
                suf[q - 1] = i as u32;
            } else {
                suf.copy_within(q..p, q + 1);
                suf[q] = i as u32;
            }
        }
    });

    // clean up counters.
    for i in (1..suf.len()).rev() {
        if suf[i] > EMPTY {
            let n = decode_counter(suf[i]);
            suf.copy_within(i - n..i, i - n + 1);
            suf[i - n] = EMPTY;
        }
    }
}

/// Put the sorted lms-suffixes, originally located in head of workspace, to their corresponding bucket tails.
#[inline]
fn put_lmssuffixes(text: &[u32], suf: &mut [u32], n: usize) {
    reset_slice(&mut suf[n..], EMPTY);

    // copy sorted lms-suffixes bucket by bucket.
    let mut prev = text.len();
    let mut m = 0;
    for i in (0..n).rev() {
        let x = suf[i] as usize;
        suf[i] = EMPTY;

        let p = text[x] as usize;
        if p != prev {
            m = 0;
        }
        suf[p - m] = x as u32;
        prev = p;
        m += 1;
    }
}

/// Induce sort all the suffixes (or lms-substrings) from the sorted lms-suffixes (or lms-characters).
#[inline]
fn induce_sort(text: &[u32], suf: &mut [u32], pipeline: &mut Pipeline, block_size: usize, left_most: bool) {
    if text.len() >= THRESHOLD_PARALLEL_INDUCE {
        // induce sort in parallel using pipeline.
        par_induce_sort(text, suf, pipeline, block_size, left_most);
    } else {
        // induce sort in serial.
        nonpar_induce_sort(text, suf, left_most);
    }
}

/// Induce sort in serial.
#[inline]
fn nonpar_induce_sort(text: &[u32], suf: &mut [u32], left_most: bool) {
    // stage 1. induce l (or lml) from lms.

    // induce from the sentinel.
    let p = text[text.len() - 1] as usize;
    if p + 1 < suf.len() && suf[p + 1] == EMPTY {
        suf[p] = encode_counter(1);
        suf[p + 1] = (text.len() - 1) as u32;
    } else {
        suf[p] = (text.len() - 1) as u32;
    }

    let mut i = 0;
    while i < suf.len() {
        if suf[i] > 0 && suf[i] < EMPTY {
            // c1 is non-empty, and has a preceding character c0.
            let j = (suf[i] - 1) as usize;
            let c0 = text[j];
            let c1 = text[j + 1];
            let p = c0 as usize;
            if c0 >= c1 {
                // c0 is l-type.
                if suf[p] < EMPTY {
                    // left shift the borrowed chunk.
                    let n = suf[..p + 1].iter().rev().take_while(|&&x| x <= EMPTY).count();
                    suf.copy_within(p + 1 - n..p + 1, p - n);
                    suf[p] = EMPTY;
                    if i > p - n {
                        // shift cursor.
                        i -= 1;
                    }
                }

                if suf[p] == EMPTY {
                    if p + 1 < suf.len() && suf[p + 1] == EMPTY {
                        suf[p] = encode_counter(1);
                        suf[p + 1] = j as u32;
                    } else {
                        suf[p] = j as u32;
                    }
                } else {
                    let q = p + 1 + decode_counter(suf[p]);
                    if q < suf.len() && suf[q] == EMPTY {
                        suf[p] -= 1;
                        suf[q] = j as u32;
                    } else {
                        suf.copy_within(p + 1..q, p);
                        suf[q - 1] = j as u32;
                        if i > p {
                            // shift cursor.
                            i -= 1;
                        }
                    }
                }

                if left_most {
                    // clear non-lml.
                    suf[i] = EMPTY;
                } else if j + 2 < text.len() {
                    // clear lms.
                    let c2 = text[j + 2];
                    if !(c1 > c2 || (c1 == c2 && i > c1 as usize)) {
                        suf[i] = EMPTY;
                    }
                }
            }
        }
        i += 1;
    }

    // clean up bucket counters.
    for i in 0..suf.len() {
        if suf[i] > EMPTY {
            let n = decode_counter(suf[i]);
            suf.copy_within(i + 1..i + 1 + n, i);
            suf[i + n] = EMPTY;
        }
    }

    // stage 2. induce s (or lms) from l (or lml).

    let mut i = text.len() - 1;
    loop {
        if suf[i] > 0 && suf[i] < EMPTY {
            // c1 is non-empty, and has a preceding character c0.
            let j = (suf[i] - 1) as usize;
            let c0 = text[j];
            let c1 = text[j + 1];
            let p = c0 as usize;
            if c0 < c1 || (c0 == c1 && i < p) {
                // c0 is s-type.
                if suf[p] < EMPTY {
                    // right shift the borrowed chunk.
                    let n = suf[p..].iter().take_while(|&&x| x < EMPTY).count();
                    suf.copy_within(p..p + n, p + 1);
                    suf[p] = EMPTY;
                    if i < p + n {
                        // shift cursor.
                        i += 1;
                    }
                }

                if suf[p] == EMPTY {
                    if p > 0 && suf[p - 1] == EMPTY {
                        suf[p] = encode_counter(1);
                        suf[p - 1] = j as u32;
                    } else {
                        suf[p] = j as u32;
                    }
                } else {
                    let q = p - decode_counter(suf[p]);
                    if q > 0 && suf[q - 1] == EMPTY {
                        suf[p] -= 1;
                        suf[q - 1] = j as u32;
                    } else {
                        suf.copy_within(q..p, q + 1);
                        suf[q] = j as u32;
                        if i < p {
                            // shift cursor.
                            i += 1;
                        }
                    }
                }

                if left_most {
                    // clear non-lms.
                    suf[i] = EMPTY;
                }
            }
        }

        // avoid underflow.
        if let Some(next_i) = i.checked_sub(1) {
            i = next_i;
        } else {
            break;
        }
    }
}

/// Induce sort in parallel using pipeline that avoids most of the data races, keeps correct
/// and improves the throughput.
///
/// Comparing to the outer-level induce pipeline, the inner-level induce pipeline is designed to
/// avoid the fetch worker and the flush worker running at the same time.
///
/// In addition, the induce routine only operates on buffers other than the actual workspace
/// which is mainly maintained by the flush worker.
///
/// To avoid data displacement on the shared read buffer, the induced data in read buffer is
/// placing as if all the bucket counters were finally eliminated.
#[inline]
fn par_induce_sort(text: &[u32], suf: &mut [u32], pipeline: &mut Pipeline, block_size: usize, left_most: bool) {
    let suf = AtomicSlice::new(suf);

    // stage 1. induce l (or lml) from lms.

    // set all the bucket counters of l-type to -1.
    foreach_typedchars(text, |_, t, c| {
        if !t.stype {
            unsafe {
                suf.set(c as usize, encode_counter(1));
            }
        }
    });

    // induce from the sentinel.
    unsafe {
        let pos = (text.len() - 1) as u32;
        let chr = text[text.len() - 1] as usize;
        if chr + 1 < suf.len() && suf.get(chr + 1) == EMPTY {
            suf.decrease(chr);
            suf.set(chr + 1, pos);
        } else {
            suf.set(chr, pos);
        }
    }

    // setup buffers.
    let owned_rbufs = OwnedRBufs::new(block_size);
    let (rbuf0, rbuf1) = owned_rbufs.as_rbufs();
    let mut buffers = (
        rbuf0,
        rbuf1,
        WBuf::with_capacity(block_size),
        WBuf::with_capacity(block_size),
    );

    // pipeline of the forward scanning.
    buffers = pipeline.induce_inner(
        false,
        suf.len(),
        block_size,
        buffers,
        |(range, rbuf): (Range<usize>, RBuf<'_>)| {
            // the fetch worker.
            unsafe {
                rbuf.fetch_ltype(text, &suf, range, left_most);
            }
            rbuf
        },
        |mut wbuf: WBuf| {
            // the flush worker.
            unsafe {
                wbuf.flush_ltype(&suf);
            }
            wbuf
        },
        |ctx| {
            // the induce routine.
            let cur_start = ctx.cur_start();
            let cur_end = ctx.cur_end();
            let next_start = ctx.next_start();
            let next_end = ctx.next_end();
            for i in cur_start..cur_end {
                if let Some((chr, pos)) = unsafe { ctx.rbuf.get(i - cur_start) } {
                    // the preceding character is already known as l-type. get its actual position,
                    // then defer a write operation.
                    // access to the atomic bucket counter here is free of conflict.
                    let offset = decode_counter(unsafe { suf.decrease(chr as usize) }) as u32;
                    let idx = chr + offset;
                    ctx.wbuf.defer_write(idx, pos);

                    // when the preceding character of the preceding character might finally located
                    // in the current or the next block, we need to test if it is l-type.
                    let final_idx = idx as usize - 1;
                    if ctx.contains(final_idx) && pos > 0 {
                        let j = (pos - 1) as usize;
                        let c0 = text[pos as usize - 1];
                        let c1 = chr;
                        if c0 >= c1 {
                            // c0 is l-type, and it finally belongs to the current or the next block.
                            unsafe {
                                // write to the shared read buffer here is free of data race, because
                                // there is at most once write for each element, and the induce
                                // routine itself is in serial.
                                if final_idx >= cur_start && final_idx < cur_end {
                                    ctx.rbuf.set_cur(final_idx - cur_start, c0, pos - 1);
                                } else if final_idx >= next_start && final_idx < next_end {
                                    ctx.rbuf.set_next(final_idx - next_start, c0, pos - 1);
                                }
                            }

                            // defer a clear operation for non-lml.
                            if left_most {
                                ctx.wbuf.defer_clear(idx - 1, pos);
                            }
                        }
                    }
                }
            }
        },
    );

    // clean up the bucket counters.
    for i in 0..suf.len() {
        let x = unsafe { suf.get(i) };
        if x > EMPTY {
            let mut j = i + decode_counter(x) as usize;
            let mut y = EMPTY;
            while j > i {
                j -= 1;
                let x = unsafe { suf.get(j) };
                unsafe {
                    suf.set(j, y);
                }
                y = x;
            }
        }
    }

    // stage 2. induce s (or lms) from l (or lml).

    // set all the bucket counters of s-type to -1.
    foreach_typedchars(text, |_, t, c| {
        if t.stype {
            unsafe {
                suf.set(c as usize, encode_counter(1));
            }
        }
    });

    // pipeline of the backward scanning.
    pipeline.induce_inner(
        true,
        suf.len(),
        block_size,
        buffers,
        |(range, rbuf): (Range<usize>, RBuf<'_>)| {
            // the fetch worker.
            unsafe {
                rbuf.fetch_stype(text, &suf, range, left_most);
            }
            rbuf
        },
        |mut wbuf: WBuf| {
            // the flush worker.
            unsafe {
                wbuf.flush_stype(&suf);
            }
            wbuf
        },
        |ctx| {
            // the induce routine.
            let cur_start = ctx.cur_start();
            let cur_end = ctx.cur_end();
            let next_start = ctx.next_start();
            let next_end = ctx.next_end();
            for i in (cur_start..cur_end).rev() {
                if let Some((chr, pos)) = unsafe { ctx.rbuf.get(i - cur_start) } {
                    // the preceding character is already known as s-type. get its actual position
                    // right shifted by one to avoid underflow, and defer a write operation.
                    // access to the atomic bucket counter here is free of conflict.
                    let offset = decode_counter(unsafe { suf.decrease(chr as usize) }) as u32;
                    let idx = chr + 1 - offset;
                    ctx.wbuf.defer_write(idx, pos);

                    // when the preceding character of the preceding character might finally located
                    // in the current or the next block, we need to test if it is s-type.
                    let final_idx = idx as usize;
                    if ctx.contains(final_idx) && pos > 0 {
                        let j = (pos - 1) as usize;
                        let c0 = text[pos as usize - 1];
                        let c1 = chr;
                        if c0 < c1 || (c0 == c1 && idx <= c0) {
                            // c0 is s-type, and it finally belongs to the current or the next block.
                            unsafe {
                                // write to the shared read buffer here is free of data race, because
                                // there is at most once write for each element, and the induce
                                // routine itself is in serial.
                                if final_idx >= cur_start && final_idx < cur_end {
                                    ctx.rbuf.set_cur(final_idx - cur_start, c0, pos - 1);
                                } else if final_idx >= next_start && final_idx < next_end {
                                    ctx.rbuf.set_next(final_idx - next_start, c0, pos - 1);
                                }
                            }

                            // defer a clear operation for non-lms.
                            if left_most {
                                ctx.wbuf.defer_clear(idx, pos);
                            }
                        }
                    }
                }
            }
        },
    );
}

/// Decode offset from bucket counter.
#[inline(always)]
fn decode_counter(i: u32) -> usize {
    (-(i as i32)) as u32 as usize
}

/// Encode offset to bucket counter.
#[inline(always)]
fn encode_counter(n: usize) -> u32 {
    (-(n as u32 as i32)) as u32
}

/// The shared read buffer that stores the l/s-typed preceding characters
/// induced by the fetch worker or the induce routine.
struct RBuf<'a> {
    cur: &'a UnsafeRBuf,
    next: &'a UnsafeRBuf,
}

impl<'a> ReadBuffer for RBuf<'a> {
    /// Clear read buffer of the current block.
    #[inline(always)]
    fn reset(&mut self, block_size: usize) {
        unsafe {
            self.cur.reset(block_size);
        }
    }
}

impl<'a> RBuf<'a> {
    /// Get element of the current block non-atomically.
    ///
    /// At most one thread could access this block before reset.
    #[inline(always)]
    pub unsafe fn get(&self, i: usize) -> Option<(u32, u32)> {
        let (chr, pos) = self.cur.get(i);
        if pos < EMPTY {
            Some((chr, pos))
        } else {
            None
        }
    }

    /// Set element of the current block non-atomically.
    ///
    /// There should be at most once write to every element in the buffer after reset.
    #[inline(always)]
    pub unsafe fn set_cur(&self, i: usize, chr: u32, pos: u32) {
        self.cur.set(i, chr, pos);
    }

    /// Set element of the next block non-atomically.
    ///
    /// There should be at most once write to every element in the buffer after reset.
    #[inline(always)]
    pub unsafe fn set_next(&self, i: usize, chr: u32, pos: u32) {
        self.next.set(i, chr, pos);
    }

    /// The fetch action that loads all the l-typed preceding charaters that finally located in the current block.
    ///
    /// It would always clear lms, and clear non-lml if needed.
    #[inline]
    pub unsafe fn fetch_ltype(&self, text: &[u32], suf: &AtomicSlice<'_, u32>, range: Range<usize>, left_most: bool) {
        let start = range.start;
        let end = range.end.saturating_add(1).min(suf.len());
        let block_size = range.end - range.start;

        suf.slice(start..end).par_iter().enumerate().for_each(|(i, x)| {
            if x == 0 || x >= EMPTY {
                return;
            }

            let j = (x - 1) as usize;
            let c0 = text[j];
            let c1 = text[j + 1];
            if c0 >= c1 {
                // c0 is l-type, test if it finally belongs to the current block.
                let left_shift = suf.get(c1 as usize) > EMPTY;
                if left_shift && i > 0 {
                    self.cur.set(i - 1, c0, x - 1);
                } else if !left_shift && i < block_size {
                    self.cur.set(i, c0, x - 1);
                } else {
                    return;
                }

                // clear non-lml and lms.
                // here is no data race, because the flush worker is not running.
                if left_most {
                    suf.set(i + start, EMPTY);
                } else if j + 2 < text.len() {
                    let c2 = text[j + 2];
                    if !(c1 > c2 || (c1 == c2 && i + start > c1 as usize)) {
                        suf.set(i + start, EMPTY);
                    }
                }
            }
        });
    }

    /// The fetch action that loads all the s-typed preceding charaters that finally located in the current block in parallel.
    ///
    /// It would clear non-lms if needed.
    #[inline]
    pub unsafe fn fetch_stype(&self, text: &[u32], suf: &AtomicSlice<'_, u32>, range: Range<usize>, left_most: bool) {
        let start = range.start.saturating_sub(1);
        let end = range.end;
        let block_size = range.end - range.start;
        let o_point = if range.start == 0 { 0 } else { 1 };

        suf.slice(start..end).par_iter().enumerate().for_each(|(i, x)| {
            if x == 0 || x >= EMPTY {
                return;
            }

            let j = (x - 1) as usize;
            let c0 = text[j];
            let c1 = text[j + 1];
            if c0 < c1 || (c0 == c1 && i + start < c0 as usize) {
                // c0 is s-type, test if it finally belongs to the current block.
                let right_shift = suf.get(c1 as usize) > EMPTY;
                if right_shift && i + 1 < block_size + o_point {
                    self.cur.set(i + 1 - o_point, c0, x - 1);
                } else if !right_shift && i >= o_point {
                    self.cur.set(i - o_point, c0, x - 1);
                } else {
                    return;
                }

                // clear non-lms directly.
                // here is no data race, because the flush worker is not running.
                if left_most {
                    suf.set(i + start, EMPTY);
                }
            }
        });
    }
}

/// Owned read buffers.
struct OwnedRBufs {
    rbuf0: UnsafeRBuf,
    rbuf1: UnsafeRBuf,
}

impl OwnedRBufs {
    /// Create owned read buffers as the backend for shared read buffers.
    #[inline]
    pub fn new(cap: usize) -> Self {
        OwnedRBufs {
            rbuf0: UnsafeRBuf::with_capacity(cap),
            rbuf1: UnsafeRBuf::with_capacity(cap),
        }
    }

    /// Get shared read buffers.
    #[inline]
    pub fn as_rbufs(&self) -> (RBuf<'_>, RBuf<'_>) {
        (
            RBuf {
                cur: &self.rbuf0,
                next: &self.rbuf1,
            },
            RBuf {
                cur: &self.rbuf1,
                next: &self.rbuf0,
            },
        )
    }
}

/// Interior mutable read buffer for a single block in the workspace.
struct UnsafeRBuf {
    buf: UnsafeCell<Vec<(u32, u32)>>,
}

unsafe impl Sync for UnsafeRBuf {}

impl UnsafeRBuf {
    pub fn with_capacity(cap: usize) -> Self {
        UnsafeRBuf {
            buf: UnsafeCell::new(Vec::with_capacity(cap)),
        }
    }

    #[inline(always)]
    pub unsafe fn reset(&self, block_size: usize) {
        let buf = self.buf.get().as_mut().unwrap();
        buf.truncate(0);
        buf.resize(block_size, (u32::MAX, EMPTY));
    }

    #[inline(always)]
    pub unsafe fn len(&self) -> usize {
        let buf = self.buf.get().as_mut().unwrap();
        buf.len()
    }

    #[inline(always)]
    pub unsafe fn get(&self, i: usize) -> (u32, u32) {
        let buf = self.buf.get().as_mut().unwrap();
        buf[i]
    }

    #[inline(always)]
    pub unsafe fn set(&self, i: usize, chr: u32, pos: u32) {
        let buf = self.buf.get().as_mut().unwrap();
        debug_assert_eq!(buf[i], (u32::MAX, EMPTY));
        buf[i] = (chr, pos);
    }
}

/// Write buffer and clear buffer for induce sorting in parallel.
struct WBuf {
    wbuf: Vec<(u32, u32)>,
    cbuf: Vec<(u32, u32)>,
}

impl WriteBuffer for WBuf {
    /// Clear the write buffer and the clear buffer.
    #[inline(always)]
    fn reset(&mut self) {
        self.wbuf.truncate(0);
        self.cbuf.truncate(0);
    }
}

impl WBuf {
    /// Create a new write buffer.
    #[inline]
    pub fn with_capacity(cap: usize) -> Self {
        WBuf {
            wbuf: Vec::with_capacity(cap),
            cbuf: Vec::with_capacity(cap),
        }
    }

    /// Defer a write operation at the actual position when flush l-type,
    /// or at the actual position right shifted by one when flush s-type (to avoid underflow).
    #[inline(always)]
    pub fn defer_write(&mut self, i: u32, x: u32) {
        self.wbuf.push((i, x));
    }

    /// Defer a clear operation at the final position after all counters eliminated.
    #[inline(always)]
    pub fn defer_clear(&mut self, i: u32, x: u32) {
        self.cbuf.push((i, x));
    }

    /// Flush write buffer of l-typed characters in parallel.
    #[inline]
    pub unsafe fn flush_ltype(&mut self, suf: &AtomicSlice<'_, u32>) {
        // flush except conflicts.
        self.wbuf.par_iter_mut().for_each(|(i, x)| {
            let idx = *i as usize;
            let pos = *x;
            if idx < suf.len() && suf.get(idx) == EMPTY {
                // write and mark done.
                suf.set(idx, pos);
                *i = u32::MAX;
            }
        });

        // shift buckets.
        self.wbuf.par_iter().for_each(|&(i, x)| {
            if i != u32::MAX {
                let mut i = i as usize;
                let mut y = x;
                while i > 0 {
                    i -= 1;
                    let x = suf.get(i);
                    suf.set(i, y);
                    y = x;
                    if x > EMPTY {
                        break;
                    }
                }
            }
        });

        // clear non-lml.
        self.cbuf.par_iter().for_each(|&(i, x)| {
            let i = i as usize;
            if suf.get(i) == x {
                suf.set(i, EMPTY);
            } else if i + 1 < suf.len() && suf.get(i + 1) == x {
                suf.set(i + 1, EMPTY);
            }
        });
    }

    /// Flush write buffer of s-typed characters in parallel.
    #[inline]
    pub unsafe fn flush_stype(&mut self, suf: &AtomicSlice<'_, u32>) {
        // flush except conflicts.
        self.wbuf.par_iter_mut().for_each(|(i, x)| {
            let idx = *i as usize;
            let pos = *x;
            if idx > 0 && suf.get(idx - 1) == EMPTY {
                // write and mark done.
                suf.set(idx - 1, pos);
                *i = u32::MAX;
            }
        });

        // shift buckets.
        self.wbuf.par_iter().for_each(|&(i, x)| {
            if i != u32::MAX {
                let mut i = i as usize;
                let mut y = x;
                while i < suf.len() {
                    let x = suf.get(i);
                    suf.set(i, y);
                    y = x;
                    if x > EMPTY {
                        break;
                    }
                    i += 1;
                }
            }
        });

        // clear non-lms.
        self.cbuf.par_iter().for_each(|&(i, x)| {
            let i = i as usize;
            if suf.get(i) == x {
                suf.set(i, EMPTY);
            } else if i > 0 && suf.get(i - 1) == x {
                suf.set(i - 1, EMPTY);
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::super::common::*;
    use super::super::pipeline::*;
    use super::super::types::*;
    use super::*;

    #[test]
    fn tablecheck_sacak32() {
        let texts: &[&[u32]] = &[
            &[0, 0, 0, 0, 0, 0],
            &[0, 0, 0, 0, 0, 1],
            &[5, 4, 3, 2, 1, 0],
            &[3, 4, 5, 2, 0, 1],
            &[2, 0, 2, 0, 2, 1, 4, 3],
            &[3, 2, 1, 3, 2, 3, 2, 1, 0, 1],
            &[2, 1, 4, 1, 1, 4, 1, 3, 1],
            &[2, 1, 1, 3, 3, 1, 1, 3, 3, 1, 2, 1],
            &[2, 2, 1, 4, 4, 1, 4, 4, 1, 3, 3, 1, 1],
            &[6, 8, 9, 5, 2, 4, 3, 0, 0, 7, 1, 2],
            &[
                1, 2, 2, 1, 1, 0, 0, 1, 1, 2, 2, 0, 0, 2, 2, 0, 1, 0, 2, 0, 1, 1, 1, 1, 2, 2, 0, 0, 2, 1, 2, 1, 1, 0,
                2, 1, 2, 2, 0, 2, 1, 1, 2, 2, 2, 1, 2, 0, 0, 1, 2, 0, 0, 0, 1, 2, 2, 2, 1, 1, 1, 1, 2, 0, 2, 1, 1, 1,
                2, 1, 0, 1,
            ],
        ];

        for &text in texts.iter() {
            assert_eq!(calc_sacak32(text), calc_naive(text));
        }
    }

    #[quickcheck]
    fn quickcheck_sacak32(text: Vec<u32>) -> bool {
        calc_sacak32(&text[..]) == calc_naive(&text[..])
    }

    #[quickcheck]
    fn quickcheck_par_induce32(text: Vec<u32>, block_size: usize) {
        if text.len() == 0 || block_size == 0 {
            return;
        }
        if let Some(text) = translate_sample(&text[..]) {
            assert_eq!(
                calc_nonpar_lms_induce(&text[..], false),
                calc_par_lms_induce(&text[..], block_size, false)
            );
            assert_eq!(
                calc_nonpar_lms_induce(&text[..], true),
                calc_par_lms_induce(&text[..], block_size, true)
            );
        }
    }

    // helper functions.

    fn calc_sacak32(text: &[u32]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        if let Some(text) = translate_sample(text) {
            sacak32(&text[..], &mut suf[..], &mut Pipeline::new());
        } else {
            saca_tiny(&text[..], &mut suf[..]);
        }
        suf
    }

    fn calc_naive(text: &[u32]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        saca_tiny(text, &mut suf[..]);
        suf
    }

    fn calc_nonpar_lms_induce(text: &[u32], left_most: bool) -> Vec<u32> {
        let mut suf = vec![0; text.len()];
        put_lmscharacters(text, &mut suf[..]);
        nonpar_induce_sort(text, &mut suf[..], left_most);
        if left_most {
            let n = compact_left_range(&mut suf[..], 0..EMPTY);
            reset_slice(&mut suf[n..], EMPTY);
        }
        suf
    }

    fn calc_par_lms_induce(text: &[u32], block_size: usize, left_most: bool) -> Vec<u32> {
        let mut suf = vec![0; text.len()];
        let mut pipeline = Pipeline::new();
        put_lmscharacters(text, &mut suf[..]);
        par_induce_sort(text, &mut suf[..], &mut pipeline, block_size, left_most);
        if left_most {
            let n = compact_left_range(&mut suf[..], 0..EMPTY);
            reset_slice(&mut suf[n..], EMPTY);
        }
        suf
    }

    fn translate_sample(text: &[u32]) -> Option<Vec<u32>> {
        let mut bkt: BTreeMap<u32, (u32, u32)> = BTreeMap::new();
        for c in text.iter().cloned() {
            if let Some((count, _)) = bkt.get_mut(&c) {
                *count += 1;
            } else {
                bkt.insert(c, (1, 0));
            }
        }
        if bkt.len() >= text.len() {
            return None;
        }

        let mut sum = 0;
        for (_, (left, right)) in &mut bkt {
            let count = *left;
            *left = sum;
            sum += count;
            *right = sum;
        }

        let mut text = Vec::from(text);
        foreach_typedchars_mut(&mut text[..], |i, t, p| {
            let c = *p;
            if !t.stype {
                *p = bkt.get(&c).unwrap().0;
            } else if t.stype {
                *p = bkt.get(&c).unwrap().1 - 1;
            }
        });
        Some(text)
    }
}
