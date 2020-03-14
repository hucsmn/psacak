use std::ops::{Index, IndexMut};

use super::common::*;
use super::sacak32::sacak32;
use super::types::*;

/// Sort lms-substrings by inducing.
const LMS_INDUCE: bool = false;

/// Sort suffix array for byte string.
#[inline(never)]
pub fn sacak8(text: &[u8], suf: &mut [u32]) {
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    let mut bkt = Buckets::new();

    // sort and rank lms-substrings.
    let (n, k) = if LMS_INDUCE {
        induce_sort_lmssubs(text, suf, &mut bkt)
    } else {
        direct_sort_lmssubs(text, suf, &mut bkt)
    };

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
    induce_lchars(text, suf, &mut bkt, false);
    induce_schars(text, suf, &mut bkt, false);
}

/// Induce sort and rank lms-substrings.
///
/// Returns count of lms-substrings (placed in the head of workspace),
/// and upper bound or alphabet size of ranks (placed in the tail of workspace).
fn induce_sort_lmssubs(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) -> (usize, usize) {
    // induce sort lms-substrings.
    bkt.init(text);
    put_lmschars(text, suf, bkt);
    induce_lchars(text, suf, bkt, true);
    induce_schars(text, suf, bkt, true);

    // collect sorted lms-substrings into the head of workspace.
    let mut n = 0;
    for i in 0..suf.len() {
        if suf[i] > 0 {
            suf[n] = suf[i];
            n += 1;
        }
    }

    // get ranks of lms-substrings into the tail of workspace.
    let k = rank_sorted_lmssubs(text, suf, n);
    (n, k)
}

/// Direct sort and rank lms-substrings.
///
/// Returns count of lms-substrings (placed in the head of workspace),
/// and upper bound or alphabet size of ranks (placed in the tail of workspace).
fn direct_sort_lmssubs(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) -> (usize, usize) {
    // make buckets for lms-substrings.
    let mut n = 0;
    for c in 1..=255u8 {
        bkt[c] = 0;
    }
    foreach_lmschars(text, |_, c| {
        bkt[c] += 1;
        n += 1;
    });
    bkt[0] *= 2;
    for c in 1..=255u8 {
        bkt[c] = bkt[c] * 2 + bkt[c - 1];
    }

    // place lms-substrings together with lengths in buckets.
    // the layout is like: [idx0, len0, idx1, len1, ...].
    let mut j = text.len();
    foreach_lmschars(text, |i, c| {
        let p = &mut bkt[c];
        *p -= 2;
        suf[p.as_index()] = i as u32; // lms-substring index.
        suf[p.as_index() + 1] = (j - i + 1) as u32; // lms-substring length.
        j = i;
    });
    bkt.init(text);

    // sort and rank lms-substrings with lengths in suf[..2*n].
    let k = sortnrank_lengthed_lmssubs(text, suf, n);
    (n, k)
}

/// Put lms-characters.
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

/// Put the sorted lms-suffixes in head of workspace to the right place.
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

/// Induce (left most) l-typed characters from left most s-type characters.
#[inline]
fn induce_lchars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
    bkt.set_head();

    // sentinel.
    let p = &mut bkt[text[text.len() - 1]];
    suf[p.as_index()] = u32::from_index(text.len() - 1);
    *p += 1;

    for i in 0..suf.len() {
        if suf[i] > 0 {
            // non-empty, and has preceding character.
            let j = (suf[i] - 1).as_index();
            let p = &mut bkt[text[j]];
            if text[j] >= text[j + 1] {
                // preceding character is l-type.
                suf[p.as_index()] = u32::from_index(j);
                *p += 1;
                if left_most || is_schar(text, suf[i].as_index()) {
                    // clean up lms-characters.
                    // when left_most is toggled, leave lml-type characters only.
                    suf[i] = 0;
                }
            }
        }
    }
}

/// Induce (left most) s-typed characters from (left most) l-type characters.
#[inline]
fn induce_schars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
    bkt.set_tail();

    for i in (0..suf.len()).rev() {
        if suf[i] > 0 {
            // non-empty, and has preceding character.
            let j = (suf[i] - 1).as_index();
            let p = &mut bkt[text[j]];
            if text[j] <= text[j + 1] && p.as_index() <= i {
                // preceding character is s-type.
                *p -= 1;
                suf[p.as_index()] = u32::from_index(j);
                if left_most {
                    // leave lms-characters only.
                    suf[i] = 0;
                }
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
    #[inline]
    pub fn new() -> Self {
        Buckets {
            ptrs: [0; 256],
            cache: [0; 257],
        }
    }

    #[inline]
    pub fn init(&mut self, text: &[u8]) {
        text.iter().for_each(|&c| self.cache[c.as_index() + 1] += 1);
        self.cache[1..].iter_mut().fold(0, |sum, p| {
            *p += sum;
            *p
        });
    }

    #[inline]
    pub fn set_head(&mut self) {
        for i in 0..256 {
            self.ptrs[i] = self.cache[i]
        }
    }

    #[inline]
    pub fn set_tail(&mut self) {
        for i in 0..256 {
            self.ptrs[i] = self.cache[i + 1]
        }
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
            &[],
            &[0],
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
                1, 2, 2, 1, 1, 0, 0, 1, 1, 2, 2, 0, 0, 2, 2, 0, 1, 0, 2, 0, 1, 1, 1, 1, 2, 2, 0, 0,
                2, 1, 2, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 1, 2, 2, 2, 1, 2, 0, 0, 1, 2, 0, 0, 0, 1, 2,
                2, 2, 1, 1, 1, 1, 2, 0, 2, 1, 1, 1, 2, 1, 0, 1,
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
