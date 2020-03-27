use std::ops::{Index, IndexMut};

use super::common::*;
use super::naming::*;
use super::sacak32::sacak32;
use super::types::*;

/// Initial level of SACA-K for byte strings.
#[inline]
pub fn sacak8(text: &[u8], suf: &mut [u32]) {
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // make buckets.
    let mut bkt = Buckets::new(text);

    // induce sort lms-substrings.
    put_lmscharacters(text, suf, &mut bkt);
    induce_sort(text, suf, &mut bkt, true);

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
    induce_sort(text, suf, &mut bkt, false);
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

/// Induce sort all the suffixes (or lms-substrings) from the sorted lms-suffixes (or lms-characters).
#[inline]
fn induce_sort(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
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
            let c1 = text[j];
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
    use super::super::types::*;
    use super::sacak8;

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
}
