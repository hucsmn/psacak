use super::types::*;

/// Switch: lms substring compare by fingerprint.
pub const LMS_FINGERPRINT: bool = true;

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
            inspect($text, &[], &[$($ptr ,)*]);
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

/// Naive suffix array construction that sorts tiny input.
pub fn saca_tiny<C, I>(text: &[C], suf: &mut [I])
where
    C: SacaChar,
    I: SacaIndex,
{
    suf[..text.len()]
        .iter_mut()
        .enumerate()
        .for_each(|(i, p)| *p = I::from_index(i));
    suf[..text.len()].sort_by(|&i, &j| {
        let i = i.as_index();
        let j = j.as_index();
        Ord::cmp(&text[i..], &text[j..])
    });
}

/// Test if character is s-typed.
#[inline(always)]
pub fn is_schar<C: SacaChar>(text: &[C], i: usize) -> bool {
    let c = text[i];
    for &next in text[i + 1..].iter() {
        if c < next {
            return true;
        }
        if c > next {
            return false;
        }
    }
    false
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

/// Enumerate characters and types, in reversed order.
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
        let left_most = stype != next_stype;
        f(i, SacaCharType { stype, left_most }, text[i]);
        stype = next_stype;
    }
    f(
        0,
        SacaCharType {
            stype,
            left_most: false,
        },
        text[0],
    );
}

/// Enumerate characters and types, in reversed order.
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
        let left_most = stype != next_stype;
        f(i, SacaCharType { stype, left_most }, &mut text[i]);
        stype = next_stype;
    }
    f(
        0,
        SacaCharType {
            stype,
            left_most: false,
        },
        &mut text[0],
    );
}

/// Specialized lms-characters enumerator, in reversed order.
#[inline(always)]
pub fn foreach_lmschars<C, F>(text: &[C], mut f: F)
where
    C: SacaChar,
    F: FnMut(usize, C),
{
    foreach_typedchars(text, |i, t, c| {
        if t.is_lms() {
            f(i, c)
        }
    });
}

/// Calculate the length of lms-substring (probably contains sentinel).
#[inline(always)]
pub fn lmssubs_getlen<C: SacaChar>(text: &[C], i: usize) -> usize {
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
pub fn lmssubs_equal<C: SacaChar>(text: &[C], i: usize, m: usize, j: usize, n: usize) -> bool {
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

/// Finger print of lms-substring.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub struct FingerPrint<X: Uint> {
    // length including sentinel.
    len: usize,
    // prefix of the lms-substring.
    fp: X,
}

impl<X: Uint> FingerPrint<X> {
    fn new(len: usize, fp: X) -> Self {
        FingerPrint { len, fp }
    }
}

/// Calculate the fingerprint of lms-substring.
#[inline(always)]
pub fn lmssubs_getfp<C, X>(text: &[C], i: usize) -> FingerPrint<X>
where
    C: SacaChar + As<X>,
    X: Uint,
{
    if i == text.len() {
        return FingerPrint::new(1, C::MAX.r#as());
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

    FingerPrint::new(n, fp)
}

/// Test if two lms-substrings of known fingerprints are equal.
#[inline(always)]
pub fn lmssubs_equalfp<C, X>(text: &[C], mut i: usize, fpi: FingerPrint<X>, mut j: usize, fpj: FingerPrint<X>) -> bool
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

/// Get ranks of the sorted lms-substrings in the head, to the tail of workspace.
#[inline]
pub fn rank_sorted_lmssubs<C, I, FP, GET, CMP>(text: &[C], suf: &mut [I], n: usize, getfp: GET, cmpfp: CMP) -> usize
where
    C: SacaChar,
    I: SacaIndex,
    FP: Copy + Eq,
    GET: Fn(&[C], usize) -> FP,
    CMP: Fn(&[C], usize, FP, usize, FP) -> bool,
{
    if n == 0 {
        return 0;
    }

    // insert ranks of lms-substrings by their occurrence order in text.
    let mut k = 1;
    let mut fp0 = getfp(text, suf[0].as_index());
    let (lmssubs, work) = suf.split_at_mut(n);
    work.iter_mut().for_each(|p| *p = I::MAX);
    work[lmssubs[0].as_index() / 2] = I::ZERO;
    for i in 1..n {
        let x1 = lmssubs[i].as_index();
        let x0 = lmssubs[i - 1].as_index();
        let fp1 = getfp(text, x1);
        if !cmpfp(text, x0, fp0, x1, fp1) {
            k += 1;
        }
        work[x1 / 2] = I::from_index(k - 1);
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

/// Debug inspector.
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
        let empty = I::ONE << (I::BIT_WIDTH - 1);
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
