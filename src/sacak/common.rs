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

    /// Overwrite destination vector with the whole slice.
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

/// Helper macro for the debug inspector.
#[macro_export]
macro_rules! inspect_here {
    ($( $fmt:expr $( , $args:expr )* );* => text: $text:expr, suf: $suf:expr $(, $ptr:expr)*) => {
        if cfg!(debug_assertions) {
            $( eprintln!($fmt $(, $args)*); )*
            inspect($text, $suf, &[$($ptr ,)*]);
        }
    };
    ($( $fmt:expr $( , $args:expr )* );* => text: $text:expr $(, $ptr:expr)*) => {
        if cfg!(debug_assertions) {
            $( eprintln!($fmt $(, $args)*); )*
            let suf: &[u32] = &[];
            inspect($text, suf, &[$($ptr ,)*]);
        }
    };
    ($( $fmt:expr $( , $args:expr )* );* => suf: $suf:expr $(, $ptr:expr)*) => {
        if cfg!(debug_assertions) {
            $( eprintln!($fmt $(, $args)*); )*
            inspect(b"", $suf, &[$($ptr ,)*]);
        }
    };
    ($( $fmt:expr $( , $args:expr )* );*) => {
        if cfg!(debug_assertions) {
            $( eprintln!($fmt $(, $args)*); )*
        }
    };
}

/// Debug inspector routine.
#[allow(unused)]
pub fn inspect<C, I>(text: &[C], suf: &[I], cursors: &[usize])
where
    C: SacaChar,
    I: SacaIndex,
{
    use term_table::row::Row;
    use term_table::table_cell::{Alignment, TableCell};
    use term_table::{Table, TableStyle};

    // table style.
    let mut table = Table::new();
    table.max_column_width = 150;
    table.style = TableStyle::blank();
    table.style.horizontal = '-';
    table.style.outer_top_horizontal = '-';
    table.style.outer_bottom_horizontal = '-';
    table.separate_rows = false;
    let align = Alignment::Center;

    // indeices.
    let n = Ord::max(text.len(), suf.len());
    let k = text.iter().map(|&c| c.as_index()).max().unwrap_or(0) + 1;
    let mut num_row = vec![TableCell::new("")];
    (0..n).for_each(|i| num_row.push(TableCell::new_with_alignment(format!("[{}]", i), 1, align)));
    table.add_row(Row::new(num_row));

    if text.len() > 0 {
        let n = text.len();

        // characters.
        let mut txt_row = vec![TableCell::new("text")];
        text.iter()
            .for_each(|&c| txt_row.push(TableCell::new_with_alignment(c, 1, align)));
        table.add_row(Row::new(txt_row));

        // types.
        let mut typ = vec![""; n];
        foreach_typedchars(text, |i, t, _| typ[i] = if t.stype { "S" } else { "L" });
        let mut typ_row = vec![TableCell::new("type")];
        typ.iter()
            .for_each(|&s| typ_row.push(TableCell::new_with_alignment(s, 1, align)));
        table.add_row(Row::new(typ_row));

        // lms marks.
        let mut lms = vec![""; n];
        foreach_typedchars(text, |i, t, _| {
            if t.is_lms() {
                lms[i] = "*"
            }
        });
        let mut lms_row = vec![TableCell::new("lms")];
        lms.iter()
            .for_each(|&s| lms_row.push(TableCell::new_with_alignment(s, 1, align)));
        table.add_row(Row::new(lms_row));
    }

    if suf.len() > 0 {
        let n = suf.len();

        // workspace.
        let mut suf_row = vec![TableCell::new("suf")];
        let empty = I::HIGHEST_BIT;
        suf.iter().for_each(|&i| {
            let val;
            if i == empty {
                val = String::from("E");
            } else if i > empty {
                val = format!("-{}", I::ONE + !i);
            } else {
                val = format!("{}", i);
            }
            suf_row.push(TableCell::new_with_alignment(val, 1, align))
        });
        table.add_row(Row::new(suf_row));

        if text.len() == text.len() {
            // bucket indicators.
            let mut bkt = vec![String::new(); n];
            let mut bkt_head = vec![0; k + 1];
            text.iter().for_each(|&c| bkt_head[c.as_index() + 1] += 1);
            bkt_head[1..].iter_mut().fold(0, |sum, p| {
                *p += sum;
                *p
            });
            let mut bkt_tail = Vec::from(&bkt_head[1..]);
            foreach_typedchars(text, |_, t, c| {
                let i = c.as_index();
                if t.stype {
                    bkt_tail[i] -= 1;
                    bkt[bkt_tail[i]] = format!("{}S", c);
                } else {
                    bkt[bkt_head[i]] = format!("{}L", c);
                    bkt_head[i] += 1;
                }
            });
            let mut bkt_row = vec![TableCell::new("bkt")];
            bkt.iter()
                .for_each(|s| bkt_row.push(TableCell::new_with_alignment(s, 1, align)));
            table.add_row(Row::new(bkt_row));
        }
    }

    // additional cursors.
    if cursors.len() > 0 {
        let mut ptr = vec![String::new(); n];
        cursors.iter().enumerate().for_each(|(i, &p)| {
            let c = std::char::from_u32('i' as u32 + i as u32).unwrap_or('!');
            ptr[p] = format!("↑\n{}", c);
        });
        let mut ptr_row = vec![TableCell::new("ptr")];
        ptr.iter()
            .for_each(|p| ptr_row.push(TableCell::new_with_alignment(p, 1, align)));
        table.add_row(Row::new(ptr_row));
    }

    eprintln!("{}", table.render());
}
