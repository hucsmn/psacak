use std::mem::swap;
use std::ops::{Index, IndexMut};
use std::sync::atomic::Ordering;

use rayon::prelude::*;

use super::common::*;
use super::naming::*;
use super::pipeline::*;
use super::sacak32::sacak32;
use super::types::*;

/// Block size for paralled induce sorting.
const BLOCK_SIZE: usize = 128 * 1024;

/// Threshold to enable paralled induce sorting.
const THRESHOLD_PARALLEL_INDUCE: usize = BLOCK_SIZE * 8;

/// Initial level of SACA-K for byte strings.
#[inline]
pub fn sacak8(text: &[u8], suf: &mut [u32]) {
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // induce sort lms-substrings.
    let mut bkt = Buckets::new(text);
    let mut lazy_pipeline = None;
    put_lmscharacters(text, suf, &mut bkt);
    induce_sort(text, suf, &mut bkt, &mut lazy_pipeline, BLOCK_SIZE, true);

    // construct subproblem, compute its suffix array, and get sorted lms-suffixes.
    let n = compact_exclude(suf, 0, false);
    let k = name_lmssubstrings(text, suf, n);
    if k < n {
        // need to solve the subproblem recursively.
        let (suf1, text1) = suf.split_at_mut(suf.len() - n);
        sacak32(text1, suf1);
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
    induce_sort(text, suf, &mut bkt, &mut lazy_pipeline, BLOCK_SIZE, false);
}

/// Put lms-characters to their corresponding bucket tails, in arbitary order.
#[inline]
fn put_lmscharacters(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) {
    bkt.set_tail();
    suf.iter_mut().for_each(|p| *p = 0);

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
    suf[n..].iter_mut().for_each(|p| *p = 0);

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
fn move_within<T: Copy>(slice: &mut [T], src: std::ops::Range<usize>, dest: usize, reset: T) {
    let (i, j, k) = (src.start, src.end, dest);
    slice.copy_within(src, dest);
    let blank = if dest > j {
        i..j
    } else if dest > i {
        i..k
    } else if k + (j - i) > i {
        k + (j - i)..j
    } else {
        i..j
    };
    slice[blank].iter_mut().for_each(|p| *p = reset);
}

/// Induce sort all the suffixes (or lms-substrings) from the sorted lms-suffixes (or lms-characters).
#[inline]
fn induce_sort(
    text: &[u8],
    suf: &mut [u32],
    bkt: &mut Buckets,
    lazy_pipeline: &mut Option<Pipeline>,
    block_size: usize,
    left_most: bool,
) {
    if text.len() > THRESHOLD_PARALLEL_INDUCE {
        if lazy_pipeline.is_none() {
            *lazy_pipeline = Some(Pipeline::new());
        }
        if let Some(pipeline) = lazy_pipeline {
            par_induce_sort(text, suf, bkt, pipeline, block_size, left_most);
        }
    } else {
        nonpar_induce_sort(text, suf, bkt, left_most);
    }
}

/// Induce sort in serial.
#[inline]
fn nonpar_induce_sort(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
    // stage 1. induce l (or lml) from lms.

    bkt.set_head();

    // the sentinel.
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
                    // only keep lml-suffixes.
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
                    // only keep lms-suffixes.
                    suf[i] = 0;
                }
            }
        }
    }
}

/// Induce sort in parallel.
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
    let prefetch_worker = |(block, mut rbuf)| {
        prefetch_procedure(text, &block, &mut rbuf);
        (block, rbuf)
    };
    let flush_worker = |mut wbuf| {
        flush_procedure(&suf, &mut wbuf);
        wbuf
    };
    pipeline.begin(prefetch_worker, flush_worker, |prefetch, flush| {
        // stage 1. induce l (or lml) from lms.

        // the sentinel.
        bkt.set_head();
        let mut prev_c0 = text[text.len() - 1];
        let mut p = bkt[prev_c0] as usize;
        unsafe {
            suf.set(p, (text.len() - 1) as u32);
            suf.fence(Ordering::SeqCst);
        }
        p += 1;

        // initialize worker states.
        let mut rbuf = ReadBuffer::with_capacity(block_size);
        let mut wbuf = WriteBuffer::with_capacity(block_size);
        let mut prefetch_state;
        let mut flush_state;

        // start prefetch B[0].
        prefetch_state = (
            unsafe { suf.slice(..block_size.min(suf.len())).into_vec() },
            ReadBuffer::with_capacity(block_size),
        );
        prefetch.start(prefetch_state);

        // start the pipeline.
        for i in 0..ceil_divide(suf.len(), block_size) {
            let start = i * block_size;
            let end = start.saturating_add(block_size).min(suf.len());
            let bound = end.saturating_add(block_size).min(suf.len());

            // wait prefetch B[i], start prefetch B[i+1].
            prefetch_state = prefetch.wait();
            swap(&mut rbuf, &mut prefetch_state.1);
            unsafe {
                suf.slice(end..bound).copy_to_vec(&mut prefetch_state.0);
            }
            prefetch.start(prefetch_state);

            // induce_l B[i].
            for i in start..end {
                let x = unsafe { suf.get(i) };
                if x > 0 {
                    // query the prefetched buffer.
                    let j = (x - 1) as usize;
                    let (c0, c1) = rbuf.get(j).unwrap_or_else(|| (text[j], text[j + 1]));

                    if c0 != prev_c0 {
                        bkt[prev_c0] = p as u32;
                        p = bkt[c0] as usize;
                        prev_c0 = c0;
                    }

                    if c0 >= c1 {
                        if p < bound {
                            // direct writes to B[i] and B[i+1].
                            unsafe { suf.set(p, j as u32) }
                        } else {
                            // deferred writes to B[i+2..].
                            wbuf.push(p as u32, j as u32);
                        }
                        p += 1;
                        if left_most {
                            unsafe { suf.set(i, 0) }
                        }
                    }
                }
            }

            // wait flush B[i+1..], start flush B[i+2..].
            if i == 0 {
                flush_state = WriteBuffer::with_capacity(block_size);
            } else {
                flush_state = flush.wait();
            }
            swap(&mut wbuf, &mut flush_state);
            flush.start(flush_state);
        }

        // wait prefetch B[$], wait flush B[$].
        prefetch_state = prefetch.wait();
        flush.wait();

        // stage 2. induce s or (lms) from l (or lml).

        // start prefetch B_rev[0].
        unsafe {
            suf.slice(suf.len().saturating_sub(block_size)..)
                .copy_to_vec(&mut prefetch_state.0);
        }
        prefetch.start(prefetch_state);

        // start the pipeline.
        bkt.set_tail();
        let mut prev_c0 = 255;
        let mut p = bkt[prev_c0] as usize;
        for i in 0..ceil_divide(suf.len(), block_size) {
            let start = suf.len() - (i * block_size);
            let end = start.saturating_sub(block_size);
            let bound = end.saturating_sub(block_size);

            // wait prefetch B_rev[i], start prefetch B_rev[i+1].
            prefetch_state = prefetch.wait();
            swap(&mut rbuf, &mut prefetch_state.1);
            unsafe {
                suf.slice(bound..end).copy_to_vec(&mut prefetch_state.0);
            }
            prefetch.start(prefetch_state);

            // induce_s B_rev[i].
            for i in (end..start).rev() {
                let x = unsafe { suf.get(i) };
                if x > 0 {
                    // query the prefetched buffer.
                    let j = (x - 1) as usize;
                    let (c0, c1) = rbuf.get(j).unwrap_or_else(|| (text[j], text[j + 1]));

                    if c0 != prev_c0 {
                        bkt[prev_c0] = p as u32;
                        p = bkt[c0] as usize;
                        prev_c0 = c0;
                    }

                    if c0 <= c1 && p <= i {
                        p -= 1;
                        if p >= bound {
                            // direct writes to B_rev[i] and B_rev[i+1].
                            unsafe { suf.set(p, j as u32) }
                        } else {
                            // deferred writes to B_rev[..i+2].
                            wbuf.push(p as u32, j as u32);
                        }
                        if left_most {
                            unsafe { suf.set(i, 0) }
                        }
                    }
                }
            }

            // wait flush B_rev[..i+1], start flush B_rev[..i+2].
            if i == 0 {
                flush_state = WriteBuffer::with_capacity(block_size);
            } else {
                flush_state = flush.wait();
            }
            swap(&mut wbuf, &mut flush_state);
            flush.start(flush_state);
        }

        // wait prefetch B[^], wait flush B[^].
        prefetch.wait();
        flush.wait();
    });
}

/// Prefetch procedure for paralleled induce sorting.
#[inline(always)]
fn prefetch_procedure(text: &[u8], block: &Vec<u32>, rbuf: &mut ReadBuffer) {
    rbuf.reset();
    rbuf.buf.par_extend(
        block
            .par_iter()
            .filter(|&&i| i > 0)
            .map(|&i| (i - 1, text[i as usize - 1], text[i as usize])),
    );
}

/// Flush procedure for paralleled induce sorting.
#[inline(always)]
fn flush_procedure<'a>(suf: &AtomicSlice<'a, u32>, wbuf: &mut WriteBuffer) {
    if wbuf.buf.len() < rayon::current_num_threads() {
        for (i, x) in wbuf.buf.iter().cloned() {
            unsafe {
                suf.set(i as usize, x);
            }
        }
        suf.fence(Ordering::SeqCst);
    } else {
        wbuf.buf
            .par_chunks(ceil_divide(wbuf.buf.len(), rayon::current_num_threads()))
            .for_each(|chunk| {
                for (i, x) in chunk.iter().cloned() {
                    unsafe {
                        suf.set(i as usize, x);
                    }
                }
            });
    }
    wbuf.reset();
}

/// Read buffer for paralleled induce sorting.
#[derive(Debug)]
struct ReadBuffer {
    buf: Vec<(u32, u8, u8)>,
    pos: usize,
}

impl ReadBuffer {
    #[inline(always)]
    pub fn with_capacity(cap: usize) -> Self {
        ReadBuffer {
            buf: Vec::with_capacity(cap),
            pos: 0,
        }
    }

    #[inline(always)]
    pub fn reset(&mut self) {
        self.buf.truncate(0);
        self.pos = 0;
    }

    #[inline(always)]
    pub fn push(&mut self, i: u32, c0: u8, c1: u8) {
        self.buf.push((i, c0, c1));
    }

    #[inline(always)]
    pub fn get(&mut self, i: usize) -> Option<(u8, u8)> {
        if self.pos >= self.buf.len() {
            return None;
        }
        let record = self.buf[self.pos];
        if i as u32 == record.0 {
            self.pos += 1;
            Some((record.1, record.2))
        } else {
            None
        }
    }
}

/// Random write buffer for parallel induce sorting.
#[derive(Debug)]
struct WriteBuffer {
    buf: Vec<(u32, u32)>,
}

impl WriteBuffer {
    pub fn with_capacity(cap: usize) -> Self {
        WriteBuffer {
            buf: Vec::with_capacity(cap),
        }
    }

    #[inline(always)]
    pub fn reset(&mut self) {
        self.buf.truncate(0);
    }

    #[inline(always)]
    pub fn push(&mut self, i: u32, x: u32) {
        self.buf.push((i, x));
    }
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

// Simple sacak8 tests.
#[cfg(test)]
mod tests {
    use super::super::common::saca_tiny;
    use super::super::pipeline::Pipeline;
    use super::super::types::*;
    use super::{nonpar_induce_sort, par_induce_sort, put_lmscharacters, sacak8, Buckets};

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

    #[test]
    fn tablecheck_induce8() {
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
            for block_size in 1..(2 * text.len()) {
                assert_eq!(calc_nonpar_lms_induce(text), calc_par_lms_induce(text, block_size));
            }
        }
    }

    #[quickcheck]
    fn quickcheck_induce8(text: Vec<u8>, block_size: usize) -> bool {
        if text.len() < 3 || block_size < 1 {
            return true;
        }
        calc_nonpar_lms_induce(&text[..]) == calc_par_lms_induce(&text[..], block_size)
    }

    // helper funtions.

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

    fn calc_nonpar_lms_induce(text: &[u8]) -> Vec<u32> {
        let mut suf = vec![0; text.len()];
        let mut bkt = Buckets::new(text);
        put_lmscharacters(text, &mut suf[..], &mut bkt);
        nonpar_induce_sort(text, &mut suf[..], &mut bkt, false);
        suf
    }

    fn calc_par_lms_induce(text: &[u8], block_size: usize) -> Vec<u32> {
        let mut suf = vec![0; text.len()];
        let mut bkt = Buckets::new(text);
        let mut pipeline = Pipeline::new();
        put_lmscharacters(text, &mut suf[..], &mut bkt);
        par_induce_sort(text, &mut suf[..], &mut bkt, &mut pipeline, block_size, false);
        suf
    }
}
