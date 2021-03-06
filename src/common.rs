use std::mem::align_of;
use std::ops::{Bound, RangeBounds};
use std::sync::atomic::Ordering;

use rayon::prelude::*;

use super::types::*;

/// Threshold to reset slice of integers in parallel.
const THRESHOLD_PARALLEL_RESET: usize = 1024 * 1024;

/// A stupid suffix array construction algorithm that handles tiny input.
pub fn saca_tiny<C, I>(text: &[C], mut suf: &mut [I])
where
    C: SacaChar,
    I: SacaIndex,
{
    debug_assert!(text.len() < u32::MAX as usize);

    suf = &mut suf[..text.len()];
    for (i, p) in suf.iter_mut().enumerate() {
        *p = I::from_index(i);
    }
    suf[..text.len()].sort_by(|&i, &j| Ord::cmp(&text[i.as_index()..], &text[j.as_index()..]));
}

/// Character type.
#[derive(Copy, Clone, Default, Debug, Eq, PartialEq)]
pub struct SacaCharType {
    pub stype: bool,
    pub left_most: bool,
}

impl SacaCharType {
    #[inline(always)]
    pub fn is_lms(self) -> bool {
        self.stype && self.left_most
    }
}

/// Enumerate characters and types, traversal in reversed order.
#[inline(always)]
pub fn foreach_typedchars<C, F>(text: &[C], mut f: F)
where
    C: SacaChar,
    F: FnMut(usize, SacaCharType, C),
{
    let n = text.len();
    if n == 0 {
        return;
    }

    let mut stype = false;
    for i in (1..n).rev() {
        let next_stype = text[i - 1] < text[i] || (text[i - 1] == text[i] && stype);
        let t = SacaCharType {
            stype,
            left_most: stype != next_stype,
        };
        f(i, t, text[i]);
        stype = next_stype;
    }
    let t = SacaCharType {
        stype,
        left_most: false,
    };
    f(0, t, text[0]);
}

/// Enumerate characters and types, traversal in reversed order.
#[inline(always)]
pub fn foreach_typedchars_mut<C, F>(text: &mut [C], mut f: F)
where
    C: SacaChar,
    F: FnMut(usize, SacaCharType, &mut C),
{
    let n = text.len();
    if n == 0 {
        return;
    }

    let mut stype = false;
    for i in (1..n).rev() {
        let next_stype = text[i - 1] < text[i] || (text[i - 1] == text[i] && stype);
        let t = SacaCharType {
            stype,
            left_most: stype != next_stype,
        };
        f(i, t, &mut text[i]);
        stype = next_stype;
    }
    let t = SacaCharType {
        stype,
        left_most: false,
    };
    f(0, t, &mut text[0]);
}

/// Specialized lms-characters enumerator, traversal in reversed order.
#[inline(always)]
pub fn foreach_lmschars<C, F>(text: &[C], mut f: F)
where
    C: SacaChar,
    F: FnMut(usize, C),
{
    if text.len() < 3 {
        return;
    }

    let mut stype = false;
    for i in (1..text.len() - 1).rev() {
        stype = text[i] < text[i + 1] || (text[i] == text[i + 1] && stype);
        if stype && text[i - 1] > text[i] {
            f(i, text[i]);
        }
    }
}

/// Compact all the elements that not eqauls to `except` to the left side of array.
///
/// Returns count of collected elements.
#[inline]
pub fn compact_left<T: Uint>(data: &mut [T], except: T) -> usize {
    let mut n = 0;
    for i in 0..data.len() {
        let x = data[i];
        if x != except {
            data[n] = x;
            n += 1;
        }
    }
    n
}

/// Compact all the elements that not eqauls to `except` to the right side of array.
///
/// Returns count of collected elements.
#[inline]
pub fn compact_right<T: Uint>(data: &mut [T], except: T) -> usize {
    let mut p = data.len();
    for i in (0..data.len()).rev() {
        let x = data[i];
        if x != except {
            p -= 1;
            data[p] = x;
        }
    }
    data.len() - p
}

/// Compact all the elements in given range of value to the left side of array.
///
/// Returns count of collected elements.
#[inline]
pub fn compact_left_range<T, R>(data: &mut [T], range: R) -> usize
where
    T: Uint,
    R: RangeBounds<T>,
{
    let (low, high) = range_to_bounds(range);

    let mut n = 0;
    for i in 0..data.len() {
        let x = data[i];
        if x >= low && x <= high {
            data[n] = x;
            n += 1;
        }
    }
    n
}

/// Reset slice of integers.
#[inline]
pub fn reset_slice<T: Uint>(slice: &mut [T], reset: T) {
    if slice.len() * T::SIZE <= THRESHOLD_PARALLEL_RESET {
        slice.iter_mut().for_each(|p| *p = reset);
    } else {
        // simply use the thread count to divide huge data for now.
        let jobs = rayon::current_num_threads();
        let chunk_size = ceil_divide(slice.len(), jobs);
        slice
            .par_chunks_mut(chunk_size)
            .for_each(|chunk| chunk.iter_mut().for_each(|p| *p = reset));
    }
}

/// Convert generic ranges to range bounds.
#[inline(always)]
fn range_to_bounds<T, R>(range: R) -> (T, T)
where
    T: Uint,
    R: RangeBounds<T>,
{
    let low = match range.start_bound() {
        Bound::Included(&x) => x,
        Bound::Excluded(&x) => x.saturating_add(T::ONE),
        Bound::Unbounded => T::ZERO,
    };
    let high = match range.end_bound() {
        Bound::Included(&x) => x,
        Bound::Excluded(&x) => x.saturating_sub(T::ONE),
        Bound::Unbounded => T::MAX,
    };
    (low, high)
}

/// Calculate `ceil(x/y)`.
#[inline(always)]
pub fn ceil_divide(x: usize, y: usize) -> usize {
    if x != 0 {
        1 + ((x - 1) / y)
    } else {
        0
    }
}

/// Mutable atomic unsigned integer slice adaptor.
#[derive(Debug)]
pub struct AtomicSlice<'a, T: Uint + HasAtomic> {
    slice: &'a [T],
}

unsafe impl<'a, T: Uint + HasAtomic> Sync for AtomicSlice<'a, T> {}

impl<'a, T: Uint + HasAtomic> AtomicSlice<'a, T> {
    /// Create new mutable atomic unsigned integer slice adaptor.
    #[inline]
    pub fn new(slice: &'a mut [T]) -> Self {
        // assert array is well-aligned, thus we could cast the pointer to perform atomic operations,
        // because atomic integers are guaranteed to have the same layout as the plain integers.
        assert_eq!(0, slice.as_ptr().align_offset(align_of::<T>()));
        AtomicSlice { slice }
    }

    /// Get slice length.
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.slice.len()
    }

    /// Make subslice.
    #[inline(always)]
    pub fn slice<R: RangeBounds<usize>>(&self, range: R) -> Self {
        let start = match range.start_bound() {
            Bound::Included(&x) => x,
            Bound::Excluded(&x) => x.saturating_add(1),
            Bound::Unbounded => 0,
        };
        match range.end_bound() {
            Bound::Included(&end) => AtomicSlice {
                slice: &self.slice[start..=end],
            },
            Bound::Excluded(&end) => AtomicSlice {
                slice: &self.slice[start..end],
            },
            Bound::Unbounded => AtomicSlice {
                slice: &self.slice[start..],
            },
        }
    }

    /// Get element atomically.
    #[inline(always)]
    pub unsafe fn get(&self, i: usize) -> T {
        let p = (&self.slice[i] as *const T as *const T::Atomic).as_ref().unwrap();
        p.load(Ordering::Relaxed)
    }

    /// Set element atomically.
    #[inline(always)]
    pub unsafe fn set(&self, i: usize, x: T) {
        let p = (&self.slice[i] as *const T as *const T::Atomic).as_ref().unwrap();
        p.store(x, Ordering::Relaxed);
    }

    /// Decrease element atomically, and get the original value.
    #[inline(always)]
    pub unsafe fn decrease(&self, i: usize) -> T {
        let p = (&self.slice[i] as *const T as *const T::Atomic).as_ref().unwrap();
        p.fetch_sub(T::ONE, Ordering::Relaxed)
    }

    // Create a parallel atomic iterator.
    #[inline]
    pub unsafe fn par_iter(&self) -> impl IndexedParallelIterator<Item = T> + 'a {
        self.slice.par_iter().map(|p| {
            let p = (p as *const T as *const T::Atomic).as_ref().unwrap();
            p.load(Ordering::Relaxed)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[quickcheck]
    fn quickcheck_foreach_typedchars(text: Vec<u8>) -> bool {
        let mut success = true;
        let mut prev_i = text.len();
        foreach_typedchars(&text[..], |i, t, c| {
            if prev_i != i + 1 {
                success = false;
            }
            if t.stype != is_stype(&text[..], i) {
                success = false;
            }
            if t.left_most != is_left_most(&text[..], i) {
                success = false;
            }
            if c != text[i] {
                success = false;
            }
            prev_i = i;
        });
        success
    }

    #[quickcheck]
    fn quickcheck_foreach_typedchars_mut(mut text: Vec<u8>) -> bool {
        let mut types = vec![SacaCharType::default(); text.len()];
        foreach_typedchars(&text[..], |i, t, _| {
            types[i] = t;
        });

        let mut success = true;
        let mut prev_i = text.len();
        foreach_typedchars_mut(&mut text[..], |i, t, _| {
            if prev_i != i + 1 {
                success = false;
            }
            if types[i] != t {
                success = false;
            }
            prev_i = i;
        });
        success
    }

    #[quickcheck]
    fn quickcheck_foreach_lmschars(text: Vec<u8>) -> bool {
        let mut lmsmap = vec![false; text.len()];
        foreach_typedchars(&text[..], |i, t, _| {
            if t.is_lms() {
                lmsmap[i] = true;
            }
        });

        let mut success = true;
        let mut prev_i = text.len();
        foreach_lmschars(&text[..], |i, _| {
            if prev_i <= i {
                success = false;
            }
            if !lmsmap[i] {
                success = false;
            }
            lmsmap[i] = false;
            prev_i = i;
        });
        if lmsmap.into_iter().any(|lms| lms) {
            success = false;
        }
        success
    }

    // helper functions.

    fn is_stype(text: &[u8], i: usize) -> bool {
        let c0 = text[i];
        for j in i + 1..text.len() {
            let c1 = text[j];
            if c0 > c1 {
                return false;
            } else if c0 < c1 {
                return true;
            }
        }
        false
    }

    fn is_left_most(text: &[u8], i: usize) -> bool {
        if i == 0 {
            return false;
        }
        let c0 = text[i - 1];
        let c1 = text[i];
        if c0 == c1 {
            false
        } else if is_stype(text, i) {
            c0 > c1
        } else {
            c0 < c1
        }
    }
}
