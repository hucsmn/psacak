use super::types::*;

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

/// A stupid suffix array construction algorithm that handles tiny input.
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
