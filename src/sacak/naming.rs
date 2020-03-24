use rayon::prelude::*;

use super::common::*;
use super::types::*;

/// Threshold for lms-substrings parallel ranking.
const PARALLEL_RANK_THRESHOLD: usize = 256;

/// Construct subproblem in tail of workspace from sorted lms-substrings in head.
#[inline]
pub fn name_lmssubs<C, I>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
{
    name_lmssubs_using::<C, I, DefaultLmsFingerprint>(text, suf, n)
}

/// Construct subproblem from sorted lms-substringsin tail of workspace from sorted lms-substrings in head,
/// using customed fingerprint comparator.
#[inline]
pub fn name_lmssubs_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
    FP: Fingerprint<C>,
{
    if n < PARALLEL_RANK_THRESHOLD {
        nonpar_make_subproblem_using::<C, I, FP>(text, suf, n)
    } else {
        par_make_subproblem_using::<C, I, FP>(text, suf, n)
    }
}

/// Construct subproblem from sorted lms-substrings in serial, using customed fingerprint comparator.
#[inline]
fn nonpar_make_subproblem_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
    FP: Fingerprint<C>,
{
    if n == 0 {
        return 0;
    }
    let (lmssubs, work) = suf.split_at_mut(n);
    work.iter_mut().for_each(|p| *p = I::MAX);

    // compare neighboring lms-substrings,
    // and compute the interim subproblem and its bucket boundaries in place.
    let mut k = 0;
    let mut h = 0;
    let mut m = 1;
    let mut i0 = lmssubs[0].as_index();
    let mut fp0 = FP::get(text, i0);
    work[i0 / 2] = I::ZERO;
    for i in 1..n {
        let i1 = lmssubs[i].as_index();
        let fp1 = FP::get(text, i1);
        if !FP::equals(text, i0, fp0, i1, fp1) {
            // reuse workspace to store bucket tails.
            lmssubs[k] = I::from_index(h + m);
            k += 1;
            h += m;
            m = 0;
        }
        m += 1;
        i0 = i1;
        fp0 = fp1;
        work[i1 / 2] = I::from_index(k);
    }
    lmssubs[k] = I::from_index(h + m);
    k += 1;

    // compact the interim subproblem into tail of workspace.
    let mut p = work.len();
    for i in (0..work.len()).rev() {
        if work[i] != I::MAX {
            p -= 1;
            work[p] = work[i];
        }
    }

    // translate characters to bucket pointers.
    foreach_typedchars_mut(&mut work[p..], |i, t, p| {
        let c = p.as_index();
        if !t.stype {
            // l-type => bucket head.
            *p = if c == 0 { I::ZERO } else { lmssubs[c - 1] };
        } else if t.stype {
            // s-type => bucket tail - 1.
            *p = lmssubs[c] - I::ONE;
        }
    });
    k
}

/// Construct subproblem from sorted lms-substrings in parallel, using customed fingerprint comparator.
#[inline]
fn par_make_subproblem_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
    FP: Fingerprint<C>,
{
    let jobs = rayon::current_num_threads();
    let chunk_size = ceil_divide(n - 1, jobs);

    // compare lms-substrings in parallel.
    let (lmssubs, work) = suf.split_at_mut(n);
    work.par_chunks_mut(ceil_divide(n, jobs))
        .for_each(|chunk| chunk.iter_mut().for_each(|p| *p = I::LOWER_BITS));
    lmssubs[1..]
        .par_chunks(chunk_size)
        .zip(work[1..n].par_chunks_mut(chunk_size))
        .enumerate()
        .for_each(|(i, (lmssubs_chunk, work_chunk))| {
            let mut i0 = lmssubs[i * chunk_size].as_index();
            let mut fp0 = FP::get(text, i0);
            for j in 0..lmssubs_chunk.len() {
                let i1 = lmssubs_chunk[j].as_index();
                let fp1 = FP::get(text, i1);
                if !FP::equals(text, i0, fp0, i1, fp1) {
                    // set highest bit if two neighboring lms-substrings are different.
                    work_chunk[j] |= I::HIGHEST_BIT;
                }
                i0 = i1;
                fp0 = fp1;
            }
        });

    // compute the interim subproblem and its bucket boundaries in place.
    let mut k = 0;
    let mut h = 0;
    let mut m = 0;
    for i in 0..n {
        let x = lmssubs[i].as_index();
        if work[i] & I::HIGHEST_BIT != I::ZERO {
            // reuse workspace to store bucket tails.
            lmssubs[k] = I::from_index(h + m);
            k += 1;
            h += m;
            m = 0;
        }
        m += 1;
        let j = x / 2;
        if j < n {
            // preserve equals bit.
            work[j] = (work[j] & I::HIGHEST_BIT) | (I::from_index(k) & I::LOWER_BITS);
        } else {
            work[j] = I::from_index(k);
        }
    }
    lmssubs[k] = I::from_index(h + m);
    k += 1;

    // compact ranks to tail.
    let mut p = work.len();
    for i in (0..work.len()).rev() {
        let x = work[i] & I::LOWER_BITS;
        // k cannot be I::LOWER_BITS, if the initial level of SACA-K use u32 or u64 index.
        if x != I::LOWER_BITS {
            p -= 1;
            work[p] = x;
        }
    }

    // translate characters to bucket pointers.
    foreach_typedchars_mut(&mut work[p..], |i, t, p| {
        let c = p.as_index();
        if !t.stype {
            // l-type => bucket head.
            *p = if c == 0 { I::ZERO } else { lmssubs[c - 1] };
        } else if t.stype {
            // s-type => bucket tail - 1.
            *p = lmssubs[c] - I::ONE;
        }
    });
    k
}

/// Calculate `⌈x/y⌉`.
#[inline(always)]
fn ceil_divide(x: usize, y: usize) -> usize {
    if x != 0 {
        1 + ((x - 1) / y)
    } else {
        0
    }
}

/// Fingerprints for comparison of lms-substrings.
pub trait Fingerprint<C: SacaChar>: Copy + Eq {
    fn get(text: &[C], i: usize) -> Self;
    fn equals(text: &[C], i: usize, fpi: Self, j: usize, fpj: Self) -> bool;
}

/// Defualt fingerprint for lms-substring comparisons.
pub type DefaultLmsFingerprint = UintFingerprint<u128>;

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

/// Permutate lms-suffixes in place, using the suffix array of subproblem in head of workspace.
#[inline]
pub fn permut_lmssufs<C, I>(text: &[C], suf: &mut [I], n: usize)
where
    C: SacaChar,
    I: SacaIndex,
{
    // get the original problem in tail of workspace.
    let mut p = suf.len();
    foreach_lmschars(text, |i, _| {
        p -= 1;
        suf[p] = I::from_index(i);
    });

    // permutate lms-substrings in place.
    for i in 0..n {
        let j = suf[i].as_index();
        suf[i] = suf[p + j];
    }
}
