use super::types::*;

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
#[inline]
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
#[inline]
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
#[inline]
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

/// Specialized rml-characters enumerator, in reversed order.
#[inline]
pub fn foreach_rmlchars<C, F>(text: &[C], mut f: F)
where
    C: SacaChar,
    F: FnMut(usize, C),
{
    if text.len() > 0 {
        f(text.len() - 1, text[text.len() - 1])
    }

    let mut prev_lms = false;
    foreach_typedchars(text, |i, t, c| {
        if t.is_lms() {
            prev_lms = true;
        } else if prev_lms {
            f(i, c);
            prev_lms = false;
        } else {
            prev_lms = false;
        }
    });
}

/// Test if two lms-substrings of known fingerprint are equal.
#[inline]
fn lmssubs_equalfp<C: SacaChar>(text: &[C], i: usize, fpi: u128, j: usize, fpj: u128) -> bool {
    if i == j {
        return true;
    }

    if fpi != fpj {
        return false;
    }

    let fptype = (fpi >> 120) as u8;
    if fptype < 0x80 {
        return true;
    }

    let width = fptype - 0x80;
    let fpvalue = (fpi & ((1 << 120) - 1));
    let n = (fpvalue >> (120 - width)) as usize;
    let m = ((120 - width) / C::BIT_WIDTH) as usize;
    lmssubs_equal(text, i + m, j + m, n - m)
}

/// Test if two lms-substrings of known length are equal.
#[inline]
fn lmssubs_equal<C: SacaChar>(text: &[C], i: usize, j: usize, n: usize) -> bool {
    let p = i + n;
    let q = j + n;
    if i == j {
        return true;
    }

    if p <= text.len() && q <= text.len() {
        &text[i..p] == &text[j..q]
    } else if p > text.len() {
        false
    } else {
        false
    }
}

/// Calculate the length of lms-substring (probably contains sentinel).
#[inline]
fn lmssubs_getlen<C: SacaChar>(text: &[C], i: usize) -> usize {
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

    // exclude trailing part of valley.
    while n > 0 && text[i + n - 1] == text[i + n - 2] {
        n -= 1;
    }

    n
}

/// Calculate the fingerprint of lms-substring (probably contains sentinel).
///
/// The MSB is the fingerprint type:
///     0x00-0x1f: following a short string (length is type),
///     0x20-0x30: following a sentinel terminated short string (length is `type - 0x10`, include sentinel),
///     0x80-0xff: following the length (bit width is `type - 0x80`), and a short prefix of the lms-substring,
///
/// Therefore, `fp0 == fp1 && fp_type < 0x80 && lmssubs_equal(text, i0, i1, fp_len)`
/// is equivalent to `the two lms-substrings are equal`.
#[inline]
fn lmssubs_getfp<C, I>(text: &[C], i: usize) -> u128
where
    C: SacaChar,
    I: SacaIndex,
{
    if i == text.len() {
        return 0x21 << 120;
    }

    let mut n = 1;
    let mut fp = text[i].as_u64() as u128;
    let mut sentinel = false;

    // upslope and plateau.
    while i + n < text.len() && text[i + n] >= text[i + n - 1] {
        if n <= 15 / C::SIZE {
            fp <<= C::BIT_WIDTH;
            fp |= text[i + n].as_u64() as u128;
        }
        n += 1;
    }

    // downslope and valley.
    while i + n < text.len() && text[i + n] <= text[i + n - 1] {
        if n <= 15 / C::SIZE {
            fp <<= C::BIT_WIDTH;
            fp |= text[i + n].as_u64() as u128;
        }
        n += 1;
    }

    if i + n == text.len() {
        // include sentinel.
        n += 1;
        sentinel = true;
    } else {
        // exclude trailing part of valley.
        while n > 0 && text[i + n - 1] == text[i + n - 2] {
            if n <= 15 / C::SIZE {
                fp >>= C::BIT_WIDTH;
            }
            n -= 1;
        }
    }

    // pack fingerprint.
    if sentinel && n <= 15 / C::SIZE + 1 {
        fp |= ((0x20 + n) as u128) << 120;
    } else if n > 15 / C::SIZE {
        fp >>= (15 / C::SIZE - (15 - I::SIZE) / C::SIZE) as u8 * C::BIT_WIDTH;
        fp |= (n as u128) << (120 - I::BIT_WIDTH);
        fp |= ((0x80 + I::BIT_WIDTH) as u128) << 120;
    } else {
        fp |= (n as u128) << 120;
    }
    fp
}

/// Debug only utility to show lms-substring fingerprint.
fn fmtfp<C: SacaChar>(fp: u128) -> String {
    let fptype = (fp >> 120) as u8;
    let fpvalue = (fp & ((1 << 120) - 1));
    if fptype < 0x20 {
        let n = fptype as usize;
        let mut x = fp & ((1 << 120) - 1);
        let mut s = vec![C::ZERO; n];
        for i in (0..n).rev() {
            s[i] = C::from_u64(x as u64);
            x >>= C::BIT_WIDTH;
        }
        format!("{:02x}:{:030x} <short {}+{:?}>", fptype, fpvalue, n, s)
    } else if fptype < 0x80 {
        let n = (fptype - 0x20) as usize;
        let mut x = fpvalue;
        let mut s = vec![C::ZERO; n - 1];
        for i in (0..n - 1).rev() {
            s[i] = C::from_u64(x as u64);
            x >>= C::BIT_WIDTH;
        }
        format!("{:02x}:{:030x} <short {}+{:?}$>", fptype, fpvalue, n, s)
    } else {
        let w = (fptype - 0x80) as u8;
        let n = (fpvalue >> (120 - w)) as usize;
        let m = ((120 - w) / C::BIT_WIDTH) as usize;
        let pre = fpvalue & ((1 << (120 - w)) - 1);
        let mut x = pre;
        let mut s = vec![C::ZERO; m];
        for i in (0..m).rev() {
            s[i] = C::from_u64(x as u64);
            x >>= C::BIT_WIDTH;
        }
        format!(
            "{:02x}:{:08x}:{:022x} <long {}+{:?}+@{}...>",
            fptype, n, pre, n, s, m
        )
    }
}

/// Get ranks of the sorted lms-substrings in the head, to the tail of workspace.
pub fn rank_sorted_lmssubs<C, I>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex,
{
    if n == 0 {
        return 0;
    }

    //inspect_here!("rank_sorted_lmssubs" => text: text, suf: suf);

    // insert ranks of lms-substrings by their occurrence order in text.
    let mut k = 1;
    let mut fp0 = lmssubs_getfp::<C, I>(text, suf[0].as_index());
    //eprintln!("  fp[{:02}]: {}", suf[0], fmtfp::<C>(fp0));
    let (lmssubs, work) = suf.split_at_mut(n);
    work.iter_mut().for_each(|p| *p = I::MAX);
    work[lmssubs[0].as_index() / 2] = I::ZERO;
    for i in 1..n {
        let x1 = lmssubs[i].as_index();
        let x0 = lmssubs[i - 1].as_index();
        let fp1 = lmssubs_getfp::<C, I>(text, x1);
        //eprintln!("  fp[{:02}]: {}", x1, fmtfp::<C>(fp1));
        if !lmssubs_equalfp(text, x0, fp0, x1, fp1) {
            //eprintln!("    -> not equal");
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
    //inspect_here!("done" => text: text, suf: suf);
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
            let mut bkt_head = vec![0; n + 1];
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
