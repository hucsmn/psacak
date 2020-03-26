use std::ops::{Bound, RangeBounds};

use super::types::*;

/* TODO: rewrite using simd, when std::simd and specialization is stable. */

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
