use std::mem::swap;
use std::ops::{Index, IndexMut, Range};

use rayon::prelude::*;

use super::common::*;
use super::naming::*;
use super::pipeline::*;
use super::sacak32::sacak32;
use super::types::*;

/// Block size for induce sorting in parallel.
const BLOCK_SIZE: usize = 128 * 1024;

/// Threshold to enable induce sorting in parallel.
const THRESHOLD_PARALLEL_INDUCE: usize = 2 * BLOCK_SIZE;

/// The outer level SACA-K algorithm for byte strings.
#[inline]
pub fn sacak8(text: &[u8], suf: &mut [u32]) {
    assert!(text.len() <= u32::MAX as usize);
    assert!(text.len() <= suf.len());

    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // induce sort lms-substrings.
    let mut bkt = Buckets::new(text);
    let mut pipeline = Pipeline::new();
    put_lmscharacters(text, suf, &mut bkt);
    induce_sort(text, suf, &mut bkt, &mut pipeline, BLOCK_SIZE, true);

    // construct the subproblem, compute its suffix array, and get sorted lms-suffixes.
    let n = compact_left(suf, 0);
    let k = name_lmssubstrings(text, suf, n);
    if k < n {
        // need to solve the subproblem recursively.
        let (suf1, text1) = suf.split_at_mut(suf.len() - n);
        sacak32(text1, suf1, &mut pipeline);
    } else {
        // the subproblem itself is the inversed suffix array.
        let (suf1, text1) = suf.split_at_mut(suf.len() - n);
        for i in 0..n {
            suf1[text1[i].as_index()] = i as u32;
        }
    }
    permutate_lmssuffixes(text, suf, n);

    // induce sort the suffix array from sorted lms-suffixes.
    put_lmssuffixes(text, suf, &mut bkt, n);
    induce_sort(text, suf, &mut bkt, &mut pipeline, BLOCK_SIZE, false);
}

/// Put lms-characters to their corresponding bucket tails, in arbitary order.
#[inline]
fn put_lmscharacters(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) {
    bkt.set_tail();
    reset_slice(suf, 0);

    let mut c_prev = text[text.len() - 1];
    let mut p = bkt[c_prev] as usize;
    foreach_lmschars(text, |i, c| {
        if c != c_prev {
            bkt[c_prev] = p as u32;
            p = bkt[c] as usize;
            c_prev = c;
        }
        p -= 1;
        suf[p] = i as u32;
    });
}

/// Put the sorted lms-suffixes, originally located in head of workspace, to their corresponding bucket tails.
#[inline]
fn put_lmssuffixes(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, mut n: usize) {
    bkt.set_tail();
    reset_slice(&mut suf[n..], 0);

    for c in (0..=255).rev() {
        let m = bkt.get_lms_count(c);
        let src = n - m..n;
        let dest = bkt[c] as usize - m;
        move_within(suf, src, dest, 0);
        n -= m;
    }
}

/// Copy within slice, and reset source area to given value.
#[inline(always)]
fn move_within<T: Uint>(slice: &mut [T], src: Range<usize>, dest: usize, reset: T) {
    let (i, j, k) = (src.start, src.end, dest);
    slice.copy_within(src, dest);
    let leave = if dest > j {
        i..j
    } else if dest > i {
        i..k
    } else if k + (j - i) > i {
        k + (j - i)..j
    } else {
        i..j
    };
    reset_slice(&mut slice[leave], reset);
}

/// Induce sort all the suffixes (or lms-substrings) from the sorted lms-suffixes (or lms-characters).
#[inline]
fn induce_sort(
    text: &[u8],
    suf: &mut [u32],
    bkt: &mut Buckets,
    pipeline: &mut Pipeline,
    block_size: usize,
    left_most: bool,
) {
    if text.len() >= THRESHOLD_PARALLEL_INDUCE {
        // induce sort in parallel using pipeline.
        par_induce_sort(text, suf, bkt, pipeline, block_size, left_most);
    } else {
        // induce sort in serial.
        nonpar_induce_sort(text, suf, bkt, left_most);
    }
}

/// Induce sort in serial.
#[inline]
fn nonpar_induce_sort(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
    // stage 1. induce l (or lml) from lms.

    // induce from the sentinel.
    bkt.set_head();
    let mut prev_c0 = text[text.len() - 1];
    let mut p = bkt[prev_c0] as usize;
    suf[p] = (text.len() - 1) as u32;
    p += 1;

    for i in 0..suf.len() {
        if suf[i] > 0 {
            // c1 is non-empty, and has a preceding character.
            let j = (suf[i] - 1) as usize;
            let c0 = text[j];
            let c1 = text[j + 1];
            if c0 != prev_c0 {
                // reduce cache misses for repetitive data.
                bkt[prev_c0] = p as u32;
                p = bkt[c0] as usize;
                prev_c0 = c0;
            }
            if c0 >= c1 {
                // c0 is l-type.
                suf[p] = j as u32;
                p += 1;
                if left_most {
                    // only keep lml.
                    suf[i] = 0;
                }
            }
        }
    }

    // stage 2. induce s or (lms) from l (or lml).

    bkt.set_tail();
    let mut prev_c0 = 255;
    let mut p = bkt[prev_c0] as usize;

    for i in (0..suf.len()).rev() {
        if suf[i] > 0 {
            // c1 is non-empty, and has a preceding character c0.
            let j = (suf[i] - 1) as usize;
            let c0 = text[j];
            let c1 = text[j + 1];
            if c0 != prev_c0 {
                // reduce cache misses for repetitive data.
                bkt[prev_c0] = p as u32;
                p = bkt[c0] as usize;
                prev_c0 = c0;
            }
            if c0 <= c1 && p <= i {
                // c0 is s-type.
                p -= 1;
                suf[p] = j as u32;
                if left_most {
                    // only keep lms.
                    suf[i] = 0;
                }
            }
        }
    }
}

/// Induce sort in parallel using pipeline.
///
/// For each stage of induce, every element in the workspace can be proved:
///   1. read only once (sequentially, considering no prefetch)
///   2. write at most once (before readed)
///   3. clear at most once (immediately after readed)
///
/// Based on these, a pipeline is designed to improve throughput and avoid data races.
#[inline]
fn par_induce_sort(
    text: &[u8],
    suf: &mut [u32],
    bkt: &mut Buckets,
    pipeline: &mut Pipeline,
    block_size: usize,
    left_most: bool,
) {
    let suf = AtomicSlice::new(suf);

    // stage 1. induce l (or lml) from lms.

    // induce from the sentinel, and setup states.
    bkt.set_head();
    let mut prev_c0 = text[text.len() - 1];
    let mut p = bkt[prev_c0] as usize;
    unsafe {
        suf.set(p, (text.len() - 1) as u32);
    }
    p += 1;

    // setup buffers.
    let mut buffers = (
        RBuf::with_capacity(block_size),
        RBuf::with_capacity(block_size),
        WBuf::with_capacity(block_size),
        WBuf::with_capacity(block_size),
    );

    // pipeline of the forward scanning.
    buffers = pipeline.induce_outer(
        false,
        suf.len(),
        block_size,
        buffers,
        |(range, mut rbuf): (Range<usize>, RBuf)| {
            // the fetch worker.
            rbuf.fetch(text, &suf, range.start..range.end);
            rbuf
        },
        |mut wbuf: WBuf| {
            // the flush worker.
            wbuf.flush(&suf);
            wbuf
        },
        |ctx| {
            // the induce routine.
            let Range { start, end } = ctx.cur_block();
            for i in start..end {
                let x = unsafe { suf.get(i) };
                if x > 0 {
                    // query the read buffer, or fallback to manually fetch text.
                    let j = (x - 1) as usize;
                    let (c0, c1) = ctx
                        .rbuf
                        .get_matched(i - start, j)
                        .unwrap_or_else(|| (text[j], text[j + 1]));

                    if c0 != prev_c0 {
                        bkt[prev_c0] = p as u32;
                        p = bkt[c0] as usize;
                        prev_c0 = c0;
                    }

                    if c0 >= c1 {
                        if ctx.contains(p) {
                            // directly write to B[i] and B[i+1].
                            unsafe { suf.set(p, j as u32) }
                        } else {
                            // push to the write buffer.
                            ctx.wbuf.defer(p as u32, j as u32);
                        }
                        p += 1;
                        if left_most {
                            unsafe { suf.set(i, 0) }
                        }
                    }
                }
            }
        },
    );

    // stage 2. induce s or (lms) from l (or lml).

    // setup states.
    bkt.set_tail();
    let mut prev_c0 = 255;
    let mut p = bkt[prev_c0] as usize;

    // pipeline of the backward scanning.
    buffers = pipeline.induce_outer(
        true,
        suf.len(),
        block_size,
        buffers,
        |(range, mut rbuf): (Range<usize>, RBuf)| {
            // the fetch worker.
            rbuf.fetch(text, &suf, range.start..range.end);
            rbuf
        },
        |mut wbuf: WBuf| {
            // the flush worker.
            wbuf.flush(&suf);
            wbuf
        },
        |ctx| {
            // the induce routine.
            let Range { start, end } = ctx.cur_block();
            for i in (start..end).rev() {
                let x = unsafe { suf.get(i) };
                if x > 0 {
                    // query the read buffer, or fallback to manually fetch text.
                    let j = (x - 1) as usize;
                    let (c0, c1) = ctx
                        .rbuf
                        .get_matched(i - start, j)
                        .unwrap_or_else(|| (text[j], text[j + 1]));

                    if c0 != prev_c0 {
                        bkt[prev_c0] = p as u32;
                        p = bkt[c0] as usize;
                        prev_c0 = c0;
                    }

                    if c0 <= c1 && p <= i {
                        p -= 1;
                        if ctx.contains(p) {
                            // directly write revB[i] and revB[i+1].
                            unsafe { suf.set(p, j as u32) }
                        } else {
                            // push to the write buffer.
                            ctx.wbuf.defer(p as u32, j as u32);
                        }
                        if left_most {
                            unsafe { suf.set(i, 0) }
                        }
                    }
                }
            }
        },
    );
}

/// Bucket pointers and lms-character counters for byte string.
struct Buckets {
    ptrs: [u32; 256],
    bounds: [u32; 257],
    lmscnts: [u32; 256],
}

impl Buckets {
    #[inline(always)]
    pub fn new(text: &[u8]) -> Self {
        let mut bkt = Buckets {
            ptrs: [0; 256],
            bounds: [0; 257],
            lmscnts: [0; 256],
        };
        foreach_typedchars(text, |_, t, c| {
            bkt.bounds[c as usize + 1] += 1;
            if t.is_lms() {
                bkt.lmscnts[c as usize] += 1;
            }
        });
        let mut p = 0;
        for i in 1..257 {
            let cnt = bkt.bounds[i];
            bkt.bounds[i] += p;
            p += cnt;
        }
        bkt
    }

    #[inline(always)]
    pub fn set_head(&mut self) {
        self.ptrs.copy_from_slice(&self.bounds[..256]);
    }

    #[inline(always)]
    pub fn set_tail(&mut self) {
        self.ptrs.copy_from_slice(&self.bounds[1..257]);
    }

    #[inline(always)]
    pub fn get_lms_count(&self, c: u8) -> usize {
        self.lmscnts[c as usize] as usize
    }
}

impl Index<u8> for Buckets {
    type Output = u32;

    #[inline(always)]
    fn index(&self, c: u8) -> &Self::Output {
        &self.ptrs[c as usize]
    }
}

impl IndexMut<u8> for Buckets {
    #[inline(always)]
    fn index_mut(&mut self, c: u8) -> &mut Self::Output {
        &mut self.ptrs[c as usize]
    }
}

/// Read buffer for induce sorting in parallel.
#[derive(Debug)]
struct RBuf {
    buf: Vec<(u8, u8, u32)>,
}

impl RBuf {
    #[inline(always)]
    pub fn with_capacity(cap: usize) -> Self {
        RBuf {
            buf: Vec::with_capacity(cap),
        }
    }

    #[inline(always)]
    pub fn get_matched(&mut self, i: usize, j_to_match: usize) -> Option<(u8, u8)> {
        if i >= self.buf.len() {
            return None;
        }
        let (c0, c1, j) = self.buf[i];
        if j as usize == j_to_match {
            Some((c0, c1))
        } else {
            None
        }
    }

    /// Fetch text in parallel to the read buffer, using cached block of non-zero indices.
    #[inline]
    pub fn fetch<'a>(&mut self, text: &[u8], suf: &AtomicSlice<'a, u32>, range: Range<usize>) {
        unsafe {
            // here is a potential data race against the induce routine and flush worker.
            // however, it is correct because the atomic integers won't be teared into bad values.
            // the only bad effect is a relatively small amount of read buffer misses.
            self.buf.par_extend(suf.slice(range).par_iter().map(|x| {
                if x > 0 {
                    let j = (x - 1) as usize;
                    let c0 = text[j];
                    let c1 = text[j + 1];
                    (c0, c1, x - 1)
                } else {
                    // treat u32::MAX as empty mark, considering j is always less than u32::MAX.
                    (0, 0, u32::MAX)
                }
            }));
        }
    }
}

impl ReadBuffer for RBuf {
    #[inline(always)]
    fn reset(&mut self, _: usize) {
        self.buf.truncate(0);
    }
}

/// Write buffer for induce sorting in parallel.
#[derive(Debug)]
struct WBuf {
    buf: Vec<(u32, u32)>,
}

impl WBuf {
    pub fn with_capacity(cap: usize) -> Self {
        WBuf {
            buf: Vec::with_capacity(cap),
        }
    }

    #[inline(always)]
    pub fn defer(&mut self, i: u32, x: u32) {
        self.buf.push((i, x));
    }

    /// Flush the write buffer in parallel.
    #[inline]
    pub fn flush<'a>(&mut self, suf: &AtomicSlice<'a, u32>) {
        // the method that spilt and flush chunks of buffer in parallel is slower than
        // simply utilizing the work-stealing thread pool to randomly flush the buffer.
        self.buf.par_iter().for_each(|&(i, x)| unsafe {
            suf.set(i as usize, x);
        });
    }
}

impl WriteBuffer for WBuf {
    #[inline(always)]
    fn reset(&mut self) {
        self.buf.truncate(0);
    }
}

#[cfg(test)]
mod tests {
    use super::super::common::saca_tiny;
    use super::super::pipeline::Pipeline;
    use super::super::types::*;
    use super::*;

    #[test]
    fn tablecheck_sacak8() {
        let texts: &[&[u8]] = &[
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
            assert_eq!(calc_sacak8(text), calc_naive8(text));
        }
    }

    #[quickcheck]
    fn quickcheck_sacak8(text: Vec<u8>) -> bool {
        calc_sacak8(&text[..]) == calc_naive8(&text[..])
    }

    #[quickcheck]
    fn quickcheck_par_induce8(text: Vec<u8>, block_size: usize) {
        if text.len() == 0 || block_size == 0 {
            return;
        }
        assert_eq!(
            calc_nonpar_lms_induce(&text[..], false),
            calc_par_lms_induce(&text[..], block_size, false)
        );
        assert_eq!(
            calc_nonpar_lms_induce(&text[..], true),
            calc_par_lms_induce(&text[..], block_size, true)
        );
    }

    // helper functions.

    fn calc_sacak8(text: &[u8]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        sacak8(text, &mut suf[..]);
        suf
    }

    fn calc_naive8(text: &[u8]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        saca_tiny(text, &mut suf[..]);
        suf
    }

    fn calc_nonpar_lms_induce(text: &[u8], left_most: bool) -> Vec<u32> {
        let mut suf = vec![0; text.len()];
        let mut bkt = Buckets::new(text);
        put_lmscharacters(text, &mut suf[..], &mut bkt);
        nonpar_induce_sort(text, &mut suf[..], &mut bkt, left_most);
        if left_most {
            // only the order of lms should be correct.
            let n = compact_left(&mut suf[..], 0);
            reset_slice(&mut suf[n..], 0);
        }
        suf
    }

    fn calc_par_lms_induce(text: &[u8], block_size: usize, left_most: bool) -> Vec<u32> {
        let mut suf = vec![0; text.len()];
        let mut bkt = Buckets::new(text);
        let mut pipeline = Pipeline::new();
        put_lmscharacters(text, &mut suf[..], &mut bkt);
        par_induce_sort(text, &mut suf[..], &mut bkt, &mut pipeline, block_size, left_most);
        if left_most {
            // only the order of lms should be correct.
            let n = compact_left(&mut suf[..], 0);
            reset_slice(&mut suf[n..], 0);
        }
        suf
    }
}
