use rayon::prelude::*;

use super::common::*;
use super::types::*;

/// Threshold to name lms-substrings in parallel.
const PARALLEL_NAME_THRESHOLD: usize = 256;

/// Fingerprints for lms-substrings comparison.
pub trait Fingerprint<C: SacaChar>: Copy + Eq {
    fn get(text: &[C], i: usize) -> Self;
    fn equals(text: &[C], i: usize, fpi: Self, j: usize, fpj: Self) -> bool;
}

/// Defualt fingerprint type for lms-substring comparisons.
pub type DefaultFingerprint = UintFingerprint<u128>;

/// Name sorted lms-substrings in head of workspace, then produce subproblem in tail.
#[inline]
pub fn name_lmssubstrings<C, I>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
{
    name_lmssubstrings_using::<C, I, DefaultFingerprint>(text, suf, n)
}

/// Name sorted lms-substrings in head of workspace, then produce subproblem in tail,
/// using customed fingerprint type.
#[inline]
pub fn name_lmssubstrings_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
    FP: Fingerprint<C>,
{
    if n < PARALLEL_NAME_THRESHOLD || n <= rayon::current_num_threads() {
        // serial version, for small amount of lms-substrings.
        nonpar_name_lmssubstrings_using::<C, I, FP>(text, suf, n)
    } else if text.len() <= I::LOWER_BITS.as_index() {
        // faster version, only works when lms-substrings <= I::MAX/2.
        fastpar_name_lmssubstrings_using::<C, I, FP>(text, suf, n)
    } else {
        // fallback version, runs slower to avoid data race.
        par_name_lmssubstrings_using::<C, I, FP>(text, suf, n)
    }
}

/// Name sorted lms-substrings in serial, using customed fingerprint type.
#[inline]
fn nonpar_name_lmssubstrings_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
    FP: Fingerprint<C>,
{
    if n == 0 {
        return 0;
    }
    let (lmssubs, work) = suf.split_at_mut(n);

    // compare lms-substrings, make the interim subproblem, and get buckets.
    let mut k = 0;
    let mut h = 0;
    let mut m = 1;
    let mut i0 = lmssubs[0].as_index();
    let mut fp0 = FP::get(text, i0);
    reset_slice(work, I::MAX);
    work[i0 / 2] = I::ZERO;
    for i in 1..n {
        let i1 = lmssubs[i].as_index();
        let fp1 = FP::get(text, i1);
        if !FP::equals(text, i0, fp0, i1, fp1) {
            // store bucket tails in suf[..k].
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

    // compact the interim subproblem, and translate characters.
    let p = work.len() - n;
    compact_right(work, I::MAX);
    translate(&mut work[p..], &lmssubs[..k]);
    k
}

/// Name sorted lms-substrings partially in parallel, using customed fingerprint type.
#[inline]
fn par_name_lmssubstrings_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
    FP: Fingerprint<C>,
{
    let jobs = rayon::current_num_threads();

    // compare lms-substrings in parallel.
    let (lmssubs, work) = suf.split_at_mut(n);
    reset_slice(work, I::LOWER_BITS);
    let chunk_size = ceil_divide(n - 1, jobs);
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
                    // store difference bit in suf[n..2*n].
                    work_chunk[j] |= I::HIGHEST_BIT;
                }
                i0 = i1;
                fp0 = fp1;
            }
        });

    // make the interim subproblem, and get buckets.
    let mut k = 0;
    let mut h = 0;
    let mut m = 0;
    for i in 0..n {
        let x = lmssubs[i].as_index();
        if work[i] & I::HIGHEST_BIT != I::ZERO {
            // store bucket tails in suf[..k].
            lmssubs[k] = I::from_index(h + m);
            k += 1;
            h += m;
            m = 0;
        }
        m += 1;
        let j = x / 2;
        if j < n {
            // preserve difference bits.
            work[j] &= I::HIGHEST_BIT;
            work[j] |= I::from_index(k) & I::LOWER_BITS;
        } else {
            // k <= I::LOWER_BITS.
            work[j] = I::from_index(k);
        }
    }
    lmssubs[k] = I::from_index(h + m);
    k += 1;

    // clean up the difference bits in suf[n..2*n].
    work[..n]
        .par_chunks_mut(ceil_divide(n, jobs))
        .for_each(|chunk| chunk.iter_mut().for_each(|p| *p &= I::LOWER_BITS));

    // compact the interim subproblem, and translate characters.
    let p = work.len() - n;
    compact_right(work, I::LOWER_BITS);
    translate(&mut work[p..], &lmssubs[..k]);
    k
}

/// Name lms-substrings in parallel, using customed fingerprint type.
#[inline]
fn fastpar_name_lmssubstrings_using<C, I, FP>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex + SacaChar,
    FP: Fingerprint<C>,
{
    debug_assert!(text.len() <= I::LOWER_BITS.as_index());

    let jobs = rayon::current_num_threads();

    // compare lms-substrings in parallel, allocating `2 * jobs * I::SIZE` extra space.
    let (lmssubs, work) = suf.split_at_mut(n);
    let chunk_size = ceil_divide(n - 1, jobs);
    let i0s = (0..ceil_divide(n - 1, chunk_size))
        .into_iter()
        .map(|i| lmssubs[i * chunk_size])
        .collect::<Vec<_>>();
    let mut k0s = lmssubs[1..]
        .par_chunks_mut(chunk_size)
        .zip(i0s.into_par_iter())
        .map(|(chunk, i0)| {
            let mut k0 = I::ZERO;
            let mut i0 = i0.as_index();
            let mut fp0 = FP::get(text, i0);
            for i in 0..chunk.len() {
                let i1 = chunk[i].as_index();
                let fp1 = FP::get(text, i1);
                if !FP::equals(text, i0, fp0, i1, fp1) {
                    // directly store difference bit in suf[..n].
                    chunk[i] |= I::HIGHEST_BIT;
                    k0 += I::ONE;
                }
                i0 = i1;
                fp0 = fp1;
            }
            k0
        })
        .collect::<Vec<_>>();
    k0s.iter_mut().fold(I::ZERO, |k0, k| {
        let k1 = k0 + *k;
        *k = k0;
        k1
    });

    // compute the interim subproblem in parallel.
    reset_slice(work, I::MAX);
    work[lmssubs[0].as_index() / 2] = I::ZERO;
    {
        let work = AtomicSlice::new(work);
        lmssubs[1..]
            .par_chunks(chunk_size)
            .zip(k0s.into_par_iter())
            .for_each(|(chunk, k0)| {
                let mut k = k0.as_index();
                for i in 0..chunk.len() {
                    if chunk[i] & I::HIGHEST_BIT != I::ZERO {
                        k += 1;
                    }
                    // provable: no element would be written twice.
                    unsafe {
                        let pos = (chunk[i] & I::LOWER_BITS).as_index() / 2;
                        work.set(pos, I::from_index(k));
                    }
                }
            });
    }

    // get buckets.
    let mut k = 0;
    let mut h = 0;
    let mut m = 0;
    for i in 0..n {
        if lmssubs[i] & I::HIGHEST_BIT != I::ZERO {
            // store bucket tails in suf[..k].
            lmssubs[k] = I::from_index(h + m);
            k += 1;
            h += m;
            m = 0;
        }
        m += 1;
    }
    lmssubs[k] = I::from_index(h + m);
    k += 1;

    // compact the interim subproblem, and translate characters.
    let p = work.len() - n;
    compact_right(work, I::MAX);
    translate(&mut work[p..], &lmssubs[..k]);
    k
}

/// Preprocess the interim subproblem by translating characters to bucket pointers.
#[inline]
fn translate<I: SacaIndex + SacaChar>(text: &mut [I], bkt_tails: &[I]) {
    foreach_typedchars_mut(text, |_, t, p| {
        let c = p.as_index();
        if !t.stype {
            // l-type => bucket head.
            *p = if c == 0 { I::ZERO } else { bkt_tails[c - 1] };
        } else if t.stype {
            // s-type => bucket tail - 1.
            *p = bkt_tails[c] - I::ONE;
        }
    });
}

// Lengths of lms-substrings as legacy fingerprints.
impl<C: SacaChar> Fingerprint<C> for usize {
    /// Calculate the length of lms-substring (probably contains sentinel).
    #[inline(always)]
    fn get(text: &[C], i: usize) -> usize {
        if i == text.len() {
            return 1;
        }

        let mut n = 1;

        while i + n < text.len() && text[i + n] >= text[i + n - 1] {
            n += 1;
        }

        while i + n < text.len() && text[i + n] <= text[i + n - 1] {
            n += 1;
        }

        // include the sentinel.
        if i + n == text.len() {
            return n + 1;
        }

        // exclude the redundant trailing s-characters.
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

        text[i..p] == text[j..q]
    }
}

/// Fingerprint that stores a short prefix of the lms-substring in a big integer.
///
/// In most cases, lms-substrings are very short.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct UintFingerprint<X: Uint> {
    /// Length of lms-substring, including the sentinel.
    len: usize,

    /// An encoded short prefix of lms-substring as integer value.
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

        while i + n < text.len() && text[i + n] >= text[i + n - 1] {
            if n < X::SIZE / C::SIZE {
                fp <<= C::BIT_WIDTH;
                fp |= text[i + n].r#as();
            }
            n += 1;
        }

        while i + n < text.len() && text[i + n] <= text[i + n - 1] {
            if n < X::SIZE / C::SIZE {
                fp <<= C::BIT_WIDTH;
                fp |= text[i + n].r#as();
            }
            n += 1;
        }

        // include the sentinel.
        if i + n == text.len() {
            n += 1;
            if n < X::SIZE / C::SIZE {
                fp <<= C::BIT_WIDTH;
                fp |= C::MAX.r#as(); // C::MAX cannot be s-type.
            }
        } else {
            // exclude the redundant trailing s-characters.
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
        text[i..p] == text[j..q]
    }
}

/// Permutate lms-suffixes in place, using the suffix array of subproblem in head of workspace.
#[inline]
pub fn permutate_lmssuffixes<C, I>(text: &[C], suf: &mut [I], n: usize)
where
    C: SacaChar,
    I: SacaIndex,
{
    // get lms-substrings ordered by the position in text.
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

#[cfg(test)]
mod tests {
    use super::super::common::*;
    use super::super::types::*;
    use super::*;

    macro_rules! quickcheck_naming {
        ($($par_name:ident $fastpar_name:ident : $chr:ty, $idx:ty, $fp:ty;)*) => {
            $(
                #[quickcheck]
                fn $par_name(text: Vec<$chr>) -> bool {
                    let (n, mut suf0) = make_workspace::<$chr, $idx>(&text[..]);
                    if n < rayon::current_num_threads() + 1 {
                        return true;
                    }
                    let mut suf1 = suf0.clone();
                    nonpar_name_lmssubstrings_using::<$chr, $idx, $fp>(&text[..], &mut suf0[..], n);
                    par_name_lmssubstrings_using::<$chr, $idx, $fp>(&text[..], &mut suf1[..], n);
                    suf0[suf0.len()-n..] == suf1[suf1.len()-n..]
                }

                #[quickcheck]
                fn $fastpar_name(text: Vec<$chr>) -> bool {
                    let (n, mut suf0) = make_workspace::<$chr, $idx>(&text[..]);
                    if n <= rayon::current_num_threads() || n > <$idx>::LOWER_BITS as usize {
                        return true;
                    }
                    let mut suf1 = suf0.clone();
                    nonpar_name_lmssubstrings_using::<$chr, $idx, $fp>(&text[..], &mut suf0[..], n);
                    fastpar_name_lmssubstrings_using::<$chr, $idx, $fp>(&text[..], &mut suf1[..], n);
                    suf0[suf0.len()-n..] == suf1[suf1.len()-n..]
                }
            )*
        };
    }

    quickcheck_naming! {
        quickcheck_par_naming_lenfp    quickcheck_fastpar_naming_lenfp:    u8,  u32, usize;
        quickcheck_par_naming8_u128fp  quickcheck_fastpar_naming8_u128fp:  u8,  u32, UintFingerprint<u128>;
        quickcheck_par_naming32_u128fp quickcheck_fastpar_naming32_u128fp: u32, u32, UintFingerprint<u128>;
    }

    // helper functions.

    fn make_workspace<C, I>(text: &[C]) -> (usize, Vec<I>)
    where
        C: SacaChar,
        I: SacaIndex,
    {
        let mut suf = Vec::with_capacity(text.len());
        foreach_lmschars(text, |i, _| {
            suf.push(I::from_index(i));
        });
        let n = suf.len();
        suf.sort_by(|&i, &j| Ord::cmp(&text[i.as_index()..], &text[j.as_index()..]));
        suf.resize(text.len(), I::ZERO);
        (n, suf)
    }
}
