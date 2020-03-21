use rayon::prelude::*;

use super::common::*;
use super::types::*;

/// Threshold to disable parallel ranking for little amount of lms-substrings.
const PARALLEL_RANK_THRESHOLD: usize = 32768;

/// Get ranks of the sorted lms-substrings originally located in the head.
///
/// Then place these ranks to the tail of workspace.
///
/// Assumes that `n <= suf.len() / 2 < I::MAX / 2`.
#[inline]
pub fn rank_lmssubs<C, I>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex,
{
    rank_lmssubs_using::<C, I, DefaultLmsFingerprint>(text, suf, n)
}

/// Fingerprints for comparison of lms-substrings.
pub trait Fingerprint<C: SacaChar>: Copy + Eq {
    fn get(text: &[C], i: usize) -> Self;
    fn equals(text: &[C], i: usize, fpi: Self, j: usize, fpj: Self) -> bool;
}

/// Defualt fingerprint for lms-substring comparisons.
pub type DefaultLmsFingerprint = UintFingerprint<u128>;

/// Get ranks using customed fingerprint originally located in the head.
///
/// Then place these ranks to the tail of workspace.
///
/// Assumes that `n <= suf.len() / 2 < I::MAX / 2`.
#[inline]
pub fn rank_lmssubs_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex,
    FP: Fingerprint<C>,
{
    par_rank_lmssubs_using::<C, I, FP>(text, suf, n)
}

/// Get ranks in parallel.
#[inline]
fn par_rank_lmssubs_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex,
    FP: Fingerprint<C>,
{
    if n < PARALLEL_RANK_THRESHOLD {
        return nonpar_rank_lmssubs_using::<C, I, FP>(text, suf, n);
    }

    // compare lms-substrings.
    let (lmssubs, work) = suf.split_at_mut(n);
    let jobs = rayon::current_num_threads();
    let chunk_size = if (n - 1) % jobs == 0 {
        (n - 1) / jobs
    } else {
        (n - 1) / jobs + 1
    };
    work.iter_mut().for_each(|p| *p = I::LOWER_BITS);
    lmssubs[1..]
        .par_chunks(chunk_size)
        .enumerate()
        .zip(work[1..n].par_chunks_mut(chunk_size))
        .for_each(|((i, lmssubs_chunk), work_chunk)| {
            let mut x0 = lmssubs[i * chunk_size].as_index();
            let mut fp0 = FP::get(text, x0);
            for j in 0..lmssubs_chunk.len() {
                let x1 = lmssubs_chunk[j].as_index();
                let fp1 = FP::get(text, x1);
                if !FP::equals(text, x0, fp0, x1, fp1) {
                    // set highest bit if two neighboring lms-substrings are different.
                    work_chunk[j] |= I::HIGHEST_BIT;
                }
                x0 = x1;
                fp0 = fp1;
            }
        });

    // insert ranks of lms-substrings to work.
    let mut k = 1;
    for i in 0..n {
        let x = lmssubs[i].as_index();
        let equals = work[i] & I::HIGHEST_BIT == I::ZERO;
        if !equals {
            k += 1;
        }
        if x / 2 < n {
            // keep comparison results;
            work[x / 2] = (work[x / 2] & I::HIGHEST_BIT) | (I::from_index(k - 1) & I::LOWER_BITS);
        } else {
            work[x / 2] = I::from_index(k - 1);
        }
    }

    // compact ranks to the tail if needed.
    if k < n {
        let mut p = work.len();
        for i in (0..work.len()).rev() {
            let x = work[i] & I::LOWER_BITS;
            if x != I::LOWER_BITS {
                p -= 1;
                work[p] = x;
            }
        }
    }
    k
}

/// Get ranks with parallel disabled.
#[inline]
fn nonpar_rank_lmssubs_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex,
    FP: Fingerprint<C>,
{
    if n == 0 {
        return 0;
    }

    // insert ranks of lms-substrings to work.
    let (lmssubs, work) = suf.split_at_mut(n);
    let mut k = 1;
    let mut x0 = lmssubs[0].as_index();
    let mut fp0 = FP::get(text, x0);
    work.iter_mut().for_each(|p| *p = I::MAX);
    work[lmssubs[0].as_index() / 2] = I::ZERO;
    for i in 1..n {
        let x1 = lmssubs[i].as_index();
        let fp1 = FP::get(text, x1);
        if !FP::equals(text, x0, fp0, x1, fp1) {
            k += 1;
        }
        work[x1 / 2] = I::from_index(k - 1);
        x0 = x1;
        fp0 = fp1;
    }

    // compact ranks to the tail if needed.
    if k < n {
        let mut p = work.len();
        for i in (0..work.len()).rev() {
            if work[i] != I::MAX {
                p -= 1;
                work[p] = work[i];
            }
        }
    }
    k
}

// Length of lms-substring as fingerprint.
impl<C: SacaChar> Fingerprint<C> for usize {
    /// Calculate the length of lms-substring (probably contains sentinel).
    #[inline(always)]
    fn get(text: &[C], i: usize) -> usize {
        if i == text.len() {
            return 1;
        }

        let mut n = 1;

        // upslope and plateau.
        while i + n < text.len() && text[i + n] >= text[i + n - 1] {
            n += 1;
        }

        // downslope and valley.
        while i + n < text.len() && text[i + n] <= text[i + n - 1] {
            n += 1;
        }

        // include sentinel.
        if i + n == text.len() {
            return n + 1;
        }

        // exclude the trailing part of valley.
        while n > 0 && text[i + n - 1] == text[i + n - 2] {
            n -= 1;
        }

        n
    }

    /// Test if two lms-substrings of known length are equal.
    #[inline(always)]
    fn equals(text: &[C], i: usize, m: usize, j: usize, n: usize) -> bool {
        let p = i + m;
        let q = j + n;
        if i == j {
            return true;
        }
        if m != n || p > text.len() || q > text.len() {
            return false;
        }

        &text[i..p] == &text[j..q]
    }
}

/// Fingerprint of a short prefix of the lms-substring storing in a big integer,
/// together with lms-substring length.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct UintFingerprint<X: Uint> {
    /// Length of lms-substring, including sentinel.
    len: usize,

    /// Fingerprint value of a short prefix of lms-substring.
    fp: X,
}

impl<C: SacaChar + As<X>, X: Uint> Fingerprint<C> for UintFingerprint<X> {
    /// Calculate the fingerprint of lms-substring.
    #[inline(always)]
    fn get(text: &[C], i: usize) -> Self {
        if i == text.len() {
            return UintFingerprint {
                len: 1,
                fp: C::MAX.r#as(),
            };
        }

        let mut n = 1;
        let mut fp = text[i].r#as();

        // upslope and plateau.
        while i + n < text.len() && text[i + n] >= text[i + n - 1] {
            if n < X::SIZE / C::SIZE {
                fp <<= C::BIT_WIDTH;
                fp |= text[i + n].r#as();
            }
            n += 1;
        }

        // downslope and valley.
        while i + n < text.len() && text[i + n] <= text[i + n - 1] {
            if n < X::SIZE / C::SIZE {
                fp <<= C::BIT_WIDTH;
                fp |= text[i + n].r#as();
            }
            n += 1;
        }

        // include sentinel.
        if i + n == text.len() {
            n += 1;
            if n < X::SIZE / C::SIZE {
                fp <<= C::BIT_WIDTH;
                fp |= C::MAX.r#as(); // C::MAX cannot be s-type.
            }
        } else {
            // exclude the trailing part of valley.
            while n > 0 && text[i + n - 1] == text[i + n - 2] {
                n -= 1;
                if n < X::SIZE / C::SIZE {
                    fp >>= C::BIT_WIDTH;
                }
            }
        }

        UintFingerprint { len: n, fp }
    }

    /// Test if two lms-substrings of known fingerprints are equal.
    #[inline(always)]
    fn equals(text: &[C], mut i: usize, fpi: Self, mut j: usize, fpj: Self) -> bool
    where
        C: SacaChar + As<X>,
        X: Uint,
    {
        if i == j {
            return true;
        }
        let p = i + fpi.len;
        let q = j + fpj.len;
        if fpi != fpj || p > text.len() || q > text.len() {
            return false;
        }

        if fpi.len <= X::SIZE / C::SIZE {
            return true;
        }
        i += X::SIZE / C::SIZE;
        j += X::SIZE / C::SIZE;
        &text[i..p] == &text[j..q]
    }
}

/// Unrank sorted lms-suffixes in place, from the suffix array of subproblem in the head of workspace.
#[inline]
pub fn unrank_lmssufs<C, I>(text: &[C], suf: &mut [I], n: usize)
where
    C: SacaChar,
    I: SacaIndex,
{
    // get the original problem in the tail of workspace.
    let mut p = suf.len();
    foreach_lmschars(text, |i, _| {
        p -= 1;
        suf[p] = I::from_index(i);
    });

    // permutate lms-substrings in place, using the suffix array of subproblem.
    for i in 0..n {
        let j = suf[i].as_index();
        suf[i] = suf[suf.len() - n + j];
    }
}
