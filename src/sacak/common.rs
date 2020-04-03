use std::mem::{align_of, size_of, transmute};
use std::ops::{Bound, RangeBounds};
use std::sync::atomic::{fence, AtomicU16, AtomicU32, AtomicU64, AtomicU8, AtomicUsize, Ordering};

use super::types::*;

/// A stupid suffix array construction algorithm that handles tiny input.
pub fn saca_tiny<C, I>(text: &[C], suf: &mut [I])
where
    C: SacaChar,
    I: SacaIndex,
{
    for i in 0..text.len() {
        suf[i] = I::from_index(i);
    }
    suf[..text.len()].sort_by(|&i, &j| Ord::cmp(&text[i.as_index()..], &text[j.as_index()..]));
}

/// Character type.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct SacaCharType {
    pub stype: bool,
    pub left_most: bool,
}

impl SacaCharType {
    #[inline(always)]
    pub fn is_lms(self) -> bool {
        self.stype && self.left_most
    }

    #[inline(always)]
    pub fn is_lml(self) -> bool {
        !self.stype && self.left_most
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

/// Compact all the elements that not eqauls to `exclude` to the left (or right) side of array.
///
/// Returns count of collected elements.
#[inline]
pub fn compact_exclude<T: Uint>(data: &mut [T], exclude: T, right_side: bool) -> usize {
    if !right_side {
        let mut n = 0;
        for i in 0..data.len() {
            let x = data[i];
            if x != exclude {
                data[n] = x;
                n += 1;
            }
        }
        n
    } else {
        let mut p = data.len();
        for i in (0..data.len()).rev() {
            let x = data[i];
            if x != exclude {
                p -= 1;
                data[p] = x;
            }
        }
        data.len() - p
    }
}

/// Compact all the elements in given range of value to the left (or right) side of array.
///
/// Returns count of collected elements.
#[inline]
pub fn compact_include<T, R>(data: &mut [T], include: R, right_side: bool) -> usize
where
    T: Uint,
    R: RangeBounds<T>,
{
    let ge = match include.start_bound() {
        Bound::Included(&x) => x,
        Bound::Excluded(&x) => x.saturating_add(T::ONE),
        Bound::Unbounded => T::ZERO,
    };
    let le = match include.end_bound() {
        Bound::Included(&x) => x,
        Bound::Excluded(&x) => x.saturating_sub(T::ONE),
        Bound::Unbounded => T::MAX,
    };
    compact_between(data, ge, le, right_side)
}

/// Compact all the elements that in range `ge..=le` to the left (or right) side of array.
///
/// Returns count of collected elements.
#[inline(always)]
fn compact_between<T: Uint>(data: &mut [T], low: T, high: T, right_side: bool) -> usize {
    if !right_side {
        let mut n = 0;
        for i in 0..data.len() {
            let x = data[i];
            if x >= low && x <= high {
                data[n] = x;
                n += 1;
            }
        }
        n
    } else {
        let mut p = data.len();
        for i in (0..data.len()).rev() {
            let x = data[i];
            if x >= low && x <= high {
                p -= 1;
                data[p] = x;
            }
        }
        data.len() - p
    }
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

/// Atomic unsigned integer array reader/writer.
#[derive(Debug)]
pub struct AtomicSlice<'a, T: Uint + HasAtomic> {
    slice: &'a [T],
}

unsafe impl<'a, T: Uint + HasAtomic> Sync for AtomicSlice<'a, T> {}

impl<'a, T: Uint + HasAtomic> AtomicSlice<'a, T> {
    /// Create an atomic unsigned integer array reader/writer from a mutable slice.
    #[inline]
    pub fn new(slice: &'a mut [T]) -> Self {
        assert_eq!(size_of::<T>(), size_of::<T::Atomic>());
        assert_eq!(0, (&slice[0] as *const T).align_offset(align_of::<T>()));
        AtomicSlice { slice }
    }

    /// Overwrite destination vector with this slice.
    ///
    /// Element copying is not guaranteed to be atomic.
    #[inline]
    pub unsafe fn copy_exclude(&self, dest: &mut Vec<T>, exclude: T) {
        dest.resize(self.len(), T::ZERO);
        dest.copy_from_slice(self.slice);
        let n = compact_exclude(&mut dest[..], exclude, false);
        dest.truncate(n);
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
        let end = match range.end_bound() {
            Bound::Included(&x) => x.saturating_add(1),
            Bound::Excluded(&x) => x,
            Bound::Unbounded => self.len(),
        };
        AtomicSlice {
            slice: &self.slice[start..end],
        }
    }

    /// Get element atomically.
    #[inline(always)]
    pub unsafe fn get(&self, i: usize) -> T {
        T::Atomic::load(unsafe { transmute(&self.slice[i]) }, Ordering::Relaxed)
    }

    /// Set element atomically.
    #[inline(always)]
    pub unsafe fn set(&self, i: usize, x: T) {
        T::Atomic::store(unsafe { transmute(&self.slice[i]) }, x, Ordering::Relaxed)
    }

    /// Explicit memory fence.
    #[inline(always)]
    pub fn fence(&self, order: Ordering) {
        std::sync::atomic::fence(order)
    }
}
