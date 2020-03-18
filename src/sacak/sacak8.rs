use std::ops::{Index, IndexMut};

use super::common::*;
use super::ranking::rank_lmssubs;
use super::sacak32::sacak32;
use super::types::*;

const FAST_INDUCE_SIZE: usize = !(1u32 << 31) as usize;

/// Sort suffix array for byte string.
#[inline]
pub fn sacak8(text: &[u8], suf: &mut [u32]) {
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // induce sort lms-substrings.
    let mut bkt = Buckets::new(text);
    let n = sort_lmssubs(text, suf, &mut bkt);

    // get ranks of lms-substrings into the tail of workspace.
    let k = rank_lmssubs(text, suf, n);

    if k < n {
        // order of lms-suffixes != order of lms-substrings
        {
            let (subsuf, subtext) = suf.split_at_mut(suf.len() - n);
            sacak32(subtext, subsuf, k);
        }

        // get the original problem.
        let mut p = suf.len();
        foreach_lmschars(text, |i, _| {
            p -= 1;
            suf[p] = u32::from_index(i);
        });

        // permutate lms-substrings by the suffix array of subproblem.
        for i in 0..n {
            let j = suf[i].as_index();
            suf[i] = suf[suf.len() - n + j];
        }
    }

    // induce sort the suffix array from sorted lms-suffixes.
    put_lmssufs(text, suf, &mut bkt, n);
    induce_all(text, suf, &mut bkt);
}

/// Induce sort all the lms-substrings into the head of workspace.
fn sort_lmssubs(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) -> usize {
    // induce sort lms-substrings from lms-characters.
    put_lmschars(text, suf, bkt);
    if text.len() <= FAST_INDUCE_SIZE {
        fast_induce_lchars(text, suf, bkt);
        fast_induce_schars(text, suf, bkt);
    } else {
        induce_lchars(text, suf, bkt, true);
        induce_schars(text, suf, bkt, true);
    }

    // collect sorted lms-substrings into the head of workspace.
    let mut n = 0;
    if text.len() <= FAST_INDUCE_SIZE {
        for i in 0..suf.len() {
            if suf[i] > 0 && suf[i] < (1 << 31) {
                suf[n] = suf[i];
                n += 1;
            }
        }
    } else {
        for i in 0..suf.len() {
            if suf[i] > 0 {
                suf[n] = suf[i];
                n += 1;
            }
        }
    }
    n
}

/// Induce sort all the suffixes from sorted lms-suffixes.
fn induce_all(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) {
    // induce sort lms-substrings from lms-characters.
    if text.len() <= FAST_INDUCE_SIZE {
        fast_induce_lchars(text, suf, bkt);
        fast_induce_schars(text, suf, bkt);

        // recover non left most mark.
        for p in suf.iter_mut() {
            if *p >= (1 << 31) {
                *p = p.wrapping_neg();
            }
        }
    } else {
        induce_lchars(text, suf, bkt, false);
        induce_schars(text, suf, bkt, false);
    }
}

/// Put lms-characters to their corresponding bucket tails, in arbitary order.
#[inline]
fn put_lmschars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) {
    bkt.set_tail();
    suf.iter_mut().for_each(|p| *p = 0);

    foreach_lmschars(text, |i, c| {
        let p = &mut bkt[c];
        *p -= 1;
        suf[p.as_index()] = u32::from_index(i);
    });
}

/// Put the sorted lms-suffixes, originally located in the head of workspace,
/// to their corresponding bucket tails.
#[inline]
fn put_lmssufs(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, n: usize) {
    bkt.set_tail();
    suf[n..].iter_mut().for_each(|p| *p = 0);

    for i in (0..n).rev() {
        let x = suf[i];
        suf[i] = 0;

        let p = &mut bkt[text[x.as_index()]];
        *p -= 1;
        suf[p.as_index()] = x;
    }
}

/// Induce l-suffixes (or lml-suffixes) from sorted lms-suffixes.
///
/// Assumes that non lms-suffixes among the input `suf` have been marked as zero.
///
/// Outputs the induced l-suffixes together with the input lms-suffixes,
/// or the induced lml-suffixes with any other suffixes marked as zero.
#[inline]
fn induce_lchars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
    bkt.set_head();

    // the sentinel.
    let p = &mut bkt[text[text.len() - 1]];
    suf[p.as_index()] = u32::from_index(text.len() - 1);
    *p += 1;

    for i in 0..suf.len() {
        if suf[i] > 0 {
            // non-empty, and has a preceding character.
            let j = (suf[i] - 1).as_index();
            let p = &mut bkt[text[j]];
            if text[j] >= text[j + 1] {
                // the preceding character is l-type.
                suf[p.as_index()] = u32::from_index(j);
                *p += 1;
                if left_most {
                    // only keep lml-suffixes.
                    suf[i] = 0;
                }
            }
        }
    }
}

/// Induce s-suffixes (or lms-suffixes) from sorted l-suffixes (or lml-suffixes).
///
/// Assumes that non l-suffixes (or non lml-suffixes) among input `suf`
/// have been marked as zero.
///
/// Outputs the induced s-suffixes together with the input l-suffixes,
/// or the induced lms-suffixes with any other suffixes marked as zero.
#[inline]
fn induce_schars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
    bkt.set_tail();

    for i in (0..suf.len()).rev() {
        if suf[i] > 0 {
            // non-empty, and has a preceding character.
            let j = (suf[i] - 1).as_index();
            let p = &mut bkt[text[j]];
            if text[j] <= text[j + 1] && p.as_index() <= i {
                // the preceding character is s-type.
                *p -= 1;
                suf[p.as_index()] = u32::from_index(j);
                if left_most {
                    // only keep lms-suffixes.
                    suf[i] = 0;
                }
            }
        }
    }
}

/// Fast induce l-suffixes from s-suffixes.
///
/// Assumes that `text.len() < (1<<31)`,
/// and non lms-suffixes among input `suf` have been marked as zero or negative.
///
/// Outputs the induced lml/l-suffixes, non lml-suffixes are marked as zero or negative.
#[inline]
fn fast_induce_lchars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) {
    bkt.set_head();

    // the sentinel.
    let p = &mut bkt[text[text.len() - 1]];
    suf[p.as_index()] = u32::from_index(text.len() - 1);
    *p += 1;

    for i in 0..suf.len() {
        if suf[i] > 0 && suf[i] < (1 << 31) {
            // non-empty, left most, and has a preceding character.
            let j = (suf[i] - 1).as_index();
            let p = &mut bkt[text[j]];
            if p.as_index() > i {
                // the preceding character is l-type.
                suf[p.as_index()] = u32::from_index(j);
                *p += 1;
                // mark non left most.
                suf[i] = suf[i].wrapping_neg();
            }
        }
    }
}

/// Fast induce s-suffixes from l-suffixes.
///
/// Assumes that `text.len() < (1<<31)`,
/// and non lml-suffixes among input `suf` have been marked as zero or negative.
///
/// Outputs the induced lms/s-suffixes, non lms-suffixes are marked as zero or negative.
#[inline]
fn fast_induce_schars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) {
    bkt.set_tail();

    for i in (0..suf.len()).rev() {
        if suf[i] > 0 && suf[i] < (1 << 31) {
            // non-empty, left most, and has a preceding character.
            let j = (suf[i] - 1).as_index();
            let p = &mut bkt[text[j]];
            if p.as_index() <= i {
                // the preceding character is s-type.
                *p -= 1;
                suf[p.as_index()] = u32::from_index(j);
                // mark non left most.
                suf[i] = suf[i].wrapping_neg();
            }
        }
    }
}

/// Byte string bucket pointers.
struct Buckets {
    ptrs: [u32; 256],
    cache: [u32; 257],
}

impl Buckets {
    #[inline(always)]
    pub fn new(text: &[u8]) -> Self {
        let mut bkt = Buckets {
            ptrs: [0; 256],
            cache: [0; 257],
        };
        text.iter().for_each(|&c| bkt.cache[c.as_index() + 1] += 1);
        bkt.cache[1..].iter_mut().fold(0, |sum, p| {
            *p += sum;
            *p
        });
        bkt
    }

    #[inline(always)]
    pub fn set_head(&mut self) {
        self.ptrs.copy_from_slice(&self.cache[..256]);
    }

    #[inline(always)]
    pub fn set_tail(&mut self) {
        self.ptrs.copy_from_slice(&self.cache[1..257]);
    }
}

impl Index<u8> for Buckets {
    type Output = u32;

    #[inline(always)]
    fn index(&self, i: u8) -> &Self::Output {
        &self.ptrs[i as usize]
    }
}

impl IndexMut<u8> for Buckets {
    #[inline(always)]
    fn index_mut(&mut self, i: u8) -> &mut Self::Output {
        &mut self.ptrs[i as usize]
    }
}

// tests for sacak8.
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

    // helpers.

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
