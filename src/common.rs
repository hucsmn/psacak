use super::types::*;

/// Naive suffix array construction algorithm that works on tiny input.
pub fn saca_tiny<C, I>(text: &[C], suf: &mut [I])
where
    C: SacaChar,
    I: SacaIndex,
{
    for i in 0..text.len() {
        suf[i] = I::from_index(i);
    }
    suf[..text.len()].sort_by(|&i, &j| {
        let i = i.as_index();
        let j = j.as_index();
        Ord::cmp(&text[i..], &text[j..])
    });
}

/// Enumerate characters and types, in reversed order.
#[inline]
pub fn foreach_typedchars<C, F>(text: &[C], mut f: F)
where
    C: SacaChar,
    F: FnMut(usize, bool, C),
{
    if text.len() < 1 {
        return;
    }

    let mut prev = text[text.len() - 1];
    let mut stype = false;
    f(text.len() - 1, false, prev);

    for (i, &c) in text[..text.len() - 1].iter().enumerate().rev() {
        stype = c < prev || (c == prev && stype);
        f(i, false, prev);
    }
}

/// Enumerate lms-characters, in reversed order.
#[inline]
pub fn foreach_lmschars<C, F>(text: &[C], mut f: F)
where
    C: SacaChar,
    F: FnMut(usize, C),
{
    if text.len() < 3 {
        return;
    }

    let mut stype = text[text.len() - 2] < text[text.len() - 1];
    for i in (1..text.len() - 1).rev() {
        let c = text[i];
        let next = text[i - 1];
        let next_stype = next < c || (next == c && stype);
        if stype && !next_stype {
            f(i, c);
        }
        stype = next_stype;
    }
}

/// Make subproblem in the tail of workspace with alphabet size `k`, from the sorted lms-substrings in the head.
///
/// Returns the alphabet size `k`.
pub fn make_subproblem<C, I>(text: &[C], suf: &mut [I], n: usize) -> usize
where
    C: SacaChar,
    I: SacaIndex,
{
    let (subs, work) = suf.split_at_mut(n);
    work.iter_mut().for_each(|p| *p = I::MAX);

    // rename lms-substrings, and do bucket sort by index value.
    let mut k = 0;
    let mut prev = text.len();
    for &i in subs.iter() {
        let x = i.as_index();
        if !lmssubstrs_equal(text, x, prev) {
            k += 1;
        }
        work[x / 2] = I::from_index(k - 1);
        prev = x;
    }

    // no need to gather subproblem if k==n.
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

/// Test if lms-substrings are equal (support sentinel).
#[inline]
fn lmssubstrs_equal<C: SacaChar>(text: &[C], mut i: usize, mut j: usize) -> bool {
    if i > j {
        std::mem::swap(&mut i, &mut j);
    }

    // typically, lms-substrings are very short
    let m = lmssubstrs_getlen(text, i);
    let n = lmssubstrs_getlen(text, j);
    m == n && j <= text.len() - n && text[i..i + m] == text[j..j + n]
}

#[inline]
fn lmssubstrs_getlen<C: SacaChar>(text: &[C], i: usize) -> usize {
    if i == text.len() {
        return 1;
    }

    let mut n = 1;

    // upslope and plateau.
    while text[i + n] >= text[i + n - 1] {
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

    assert_eq!(text.len(), suf.len());
    let n = text.len();

    let mut table = Table::new();
    let align = Alignment::Right;
    table.style = TableStyle::blank();

    let mut num_row = vec![TableCell::new("")];
    (0..n).for_each(|i| num_row.push(TableCell::new_with_alignment(i, 1, align)));
    table.add_row(Row::new(num_row));

    let mut txt_row = vec![TableCell::new("text")];
    text.iter()
        .for_each(|&c| txt_row.push(TableCell::new_with_alignment(c, 1, align)));
    table.add_row(Row::new(txt_row));

    let mut suf_row = vec![TableCell::new("suf")];
    suf.iter()
        .for_each(|&i| suf_row.push(TableCell::new_with_alignment(i, 1, align)));
    table.add_row(Row::new(suf_row));

    let mut lms = vec![""; n];
    foreach_lmschars(text, |i, _| lms[i] = "*");
    let mut lms_row = vec![TableCell::new("lms")];
    lms.iter()
        .for_each(|&s| lms_row.push(TableCell::new_with_alignment(s, 1, align)));
    table.add_row(Row::new(lms_row));

    if cursors.len() > 0 {
        let mut ptr = vec![""; n];
        cursors.iter().for_each(|&p| ptr[p] = "^");
        let mut ptr_row = vec![TableCell::new("ptr")];
        ptr.iter()
            .for_each(|&p| ptr_row.push(TableCell::new_with_alignment(p, 1, align)));
        table.add_row(Row::new(ptr_row));
    }

    eprintln!("{}", table.render());
}
