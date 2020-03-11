use std::ops::{Index, IndexMut};

use super::common::*;
use super::sacak32::sacak32;
use super::types::*;

/// Sort suffix array for byte string.
pub fn sacak8(text: &[u8], suf: &mut [u32]) {
    debug_assert!(text.len() <= suf.len());
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // induce sort lms-substrings.
    let mut bkt = Buckets::new(text);
    put_lmschars(text, suf, &mut bkt);
    induce_lchars(text, suf, &mut bkt, true);
    induce_schars(text, suf, &mut bkt, true);

    // collect sorted lms-substrings into the head of workspace.
    let mut n = 0;
    for i in 0..suf.len() {
        if suf[i] > 0 {
            suf[n] = suf[i];
            n += 1;
        }
    }

    // sort lms-substrings lexicographically in the head of workspace.
    let k = make_subproblem(text, suf, n);
    if k < n {
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

/// Put lms-characters.
fn put_lmschars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets) {
    bkt.reset_tail();
    suf.iter_mut().for_each(|p| *p = 0);

    foreach_lmschars(text, |i, c| {
        let p = &mut bkt[c.as_index()];
        *p -= 1;
        suf[*p] = u32::from_index(i);
    });
}

/// Put the sorted lms-suffixes in head of workspace to the right place.
fn put_lmssufs(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, n: usize) {
    bkt.reset_tail();
    suf[n..].iter_mut().for_each(|p| *p = 0);

    for i in (0..n).rev() {
        let x = suf[i];
        suf[i] = 0;

        let p = &mut bkt[text[x.as_index()].as_index()];
        *p -= 1;
        suf[*p] = x;
    }
}

/// Induce (left most) l-typed characters from left most s-type characters.
fn induce_lchars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
    bkt.reset_head();

    // sentinel.
    let p = &mut bkt[text[text.len() - 1].as_index()];
    suf[*p] = u32::from_index(text.len() - 1);
    *p += 1;

    for i in 0..suf.len() {
        if suf[i] > 0 {
            // non-empty, and has preceding character.
            let j = (suf[i] - 1).as_index();
            let p = &mut bkt[text[j].as_index()];
            if text[j] >= text[j + 1] {
                // preceding character is l-type.
                suf[*p] = u32::from_index(j);
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
fn induce_schars(text: &[u8], suf: &mut [u32], bkt: &mut Buckets, left_most: bool) {
    bkt.reset_tail();

    for i in (0..suf.len()).rev() {
        if suf[i] > 0 {
            // non-empty, and has preceding character.
            let j = (suf[i] - 1).as_index();
            let p = &mut bkt[text[j].as_index()];
            if text[j] <= text[j + 1] {
                //&& *p <= i
                // preceding character is s-type.
                *p -= 1;
                suf[*p] = u32::from_index(j);
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
    ptrs: [usize; 256],
    cache: [usize; 257],
}

impl Buckets {
    #[inline]
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

    #[inline]
    pub fn reset_head(&mut self) {
        for i in 0..256 {
            self.ptrs[i] = self.cache[i]
        }
    }

    #[inline]
    pub fn reset_tail(&mut self) {
        for i in 0..256 {
            self.ptrs[i] = self.cache[i + 1]
        }
    }
}

impl Index<usize> for Buckets {
    type Output = usize;

    #[inline(always)]
    fn index(&self, i: usize) -> &Self::Output {
        &self.ptrs[i]
    }
}

impl IndexMut<usize> for Buckets {
    #[inline(always)]
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.ptrs[i]
    }
}

#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::super::common::saca_tiny;
    use super::sacak8;

    #[test]
    fn tablecheck_sacak8() {
        let texts: &[&[u8]] = &[
            &[],
            &[0],
            &[0, 0, 0, 0, 0, 0],
            &[2, 0, 2, 0, 2, 1, 4, 3],
            &[3, 2, 1, 3, 2, 3, 2, 1, 0, 1],
            &[2, 1, 4, 1, 1, 4, 1, 3, 1],
            &[2, 1, 1, 3, 3, 1, 1, 3, 3, 1, 2, 1],
            &[2, 2, 1, 4, 4, 1, 4, 4, 1, 3, 3, 1, 1],
            &[
                1, 2, 2, 1, 1, 0, 0, 1, 1, 2, 2, 0, 0, 2, 2, 0, 1, 0, 2, 0, 1, 1, 1, 1, 2, 2, 0, 0, 2,
                1, 2, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 1, 2, 2, 2, 1, 2, 0, 0, 1, 2, 0, 0, 0, 1, 2, 2, 2,
                1, 1, 1, 1, 2, 0, 2, 1, 1, 1, 2, 1, 0, 1,
            ],
        ];

        for &text in texts.iter() {
            assert_eq!(sacak(text), naive(text));
        }
    }

    quickcheck! {
        fn quickcheck_sacak8(text: Vec<u8>) -> bool {
            naive(&text[..]) == sacak(&text[..])
        }
    }

    fn sacak(text: &[u8]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        sacak8(text, &mut suf[..]);
        suf
    }

    fn naive(text: &[u8]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        saca_tiny(text, &mut suf[..]);
        suf
    }
}
