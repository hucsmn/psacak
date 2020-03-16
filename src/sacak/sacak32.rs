use super::common::*;
use super::types::*;

/// Empty mark in the workspace.
const EMPTY: u32 = 1 << 31;

/// Sort suffix array for integer string (alphabet size k).
pub fn sacak32(text: &mut [u32], suf: &mut [u32], k: usize) {
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // make bucket pointers by rewriting text.
    transform_text(text, suf, k);

    // induce sort lms-substrings.
    put_lmschars(text, suf);
    induce_lchars(text, suf, true);
    induce_schars(text, suf, true);

    // collect sorted lms-substrings into the head of workspace.
    let mut n = 0;
    for i in 0..suf.len() {
        if suf[i] < EMPTY {
            suf[n] = suf[i];
            n += 1;
        }
    }

    // get ranks of lms-substrings into the tail of workspace.
    let k = if LMS_FINGERPRINT {
        rank_sorted_lmssubs(text, suf, n, lmssubs_getfp::<u32, u128>, lmssubs_equalfp::<u32, u128>)
    } else {
        rank_sorted_lmssubs(text, suf, n, lmssubs_getlen, lmssubs_equal)
    };

    if k < n {
        // order of lms-suffixes != order of lms-substrings
        {
            let (subsuf, subtext) = suf.split_at_mut(suf.len() - n);
            sacak32(subtext, subsuf, k);
        }

        // get the original problem.
        let mut p = suf.len();
        foreach_lmschars(text, |i, _| {
            p -= 1;
            suf[p] = u32::from_index(i);
        });

        // permutate lms-substrings by the suffix array of subproblem.
        for i in 0..n {
            let j = suf[i].as_index();
            suf[i] = suf[suf.len() - n + j];
        }
    }

    // induce sort the suffix array from sorted lms-suffixes.
    put_lmssufs(text, suf, n);
    induce_lchars(text, suf, false);
    induce_schars(text, suf, false);
}

/// Map integer characters to bucket pointers, keeping character types unchanged.
#[inline(never)]
fn transform_text(text: &mut [u32], suf: &mut [u32], k: usize) {
    // make bucket.
    suf[..k].iter_mut().for_each(|p| *p = 0);
    text.iter().for_each(|&c| suf[c.as_index()] += 1);
    suf[..k].iter_mut().fold(0, |sum, p| {
        *p += sum;
        *p
    });

    // l-type => bucket_head, s-type => bucket_tail - 1.
    foreach_typedchars_mut(text, |i, t, p| {
        let c = p.as_index();
        if !t.stype {
            *p = if c > 0 { suf[c - 1] } else { 0 };
        } else {
            *p = suf[c] - 1;
        }
    })
}

/// Put lms-characters.
#[inline(never)]
fn put_lmschars(text: &[u32], suf: &mut [u32]) {
    suf.iter_mut().for_each(|p| *p = EMPTY);

    foreach_lmschars(text, |i, c| {
        let p = c.as_index();
        if suf[p] < EMPTY {
            // right shift the borrowed chunk.
            let n = suf[p..].iter().take_while(|&&x| x < EMPTY).count();
            suf.copy_within(p..p + n, p + 1);
            suf[p] = EMPTY;
        }

        if suf[p] == EMPTY {
            if p > 0 && suf[p - 1] == EMPTY {
                suf[p] = from_counter(1);
                suf[p - 1] = u32::from_index(i);
            } else {
                suf[p] = u32::from_index(i);
            }
        } else {
            let q = p - as_counter(suf[p]);
            if q > 0 && suf[q - 1] == EMPTY {
                suf[p] -= 1;
                suf[q - 1] = u32::from_index(i);
            } else {
                suf.copy_within(q..p, q + 1);
                suf[q] = u32::from_index(i);
            }
        }
    });

    // clean up counters.
    for i in (1..suf.len()).rev() {
        if suf[i] > EMPTY {
            let n = as_counter(suf[i]);
            suf.copy_within(i - n..i, i - n + 1);
            suf[i - n] = EMPTY;
        }
    }
}

/// Put the sorted lms-suffixes in head of workspace to the right place.
#[inline(never)]
fn put_lmssufs(text: &[u32], suf: &mut [u32], n: usize) {
    suf[n..].iter_mut().for_each(|p| *p = EMPTY);

    // copy sorted lms-suffixes bucket by bucket.
    let mut prev = text.len();
    let mut m = 0;
    for i in (0..n).rev() {
        let x = suf[i];
        suf[i] = EMPTY;

        let p = text[x.as_index()].as_index();
        if p != prev {
            m = 0;
        }
        suf[p - m] = x;
        prev = p;
        m += 1;
    }
}

/// Induce (left most) l-typed characters from left most s-type characters.
#[inline(never)]
fn induce_lchars(text: &[u32], suf: &mut [u32], left_most: bool) {
    // sentinel.
    let p = text[text.len() - 1].as_index();
    if p + 1 < suf.len() && suf[p + 1] == EMPTY {
        suf[p] = from_counter(1);
        suf[p + 1] = u32::from_index(text.len() - 1);
    } else {
        suf[p] = u32::from_index(text.len() - 1);
    }

    let mut i = 0;
    while i < suf.len() {
        if suf[i] > 0 && suf[i] < EMPTY {
            // text[suf[i]] is non-empty, and has preceding character.
            let j = (suf[i] - 1).as_index();
            let clean = is_schar(text, suf[i].as_index()); // we need manually clean up lms.
            if text[j] >= text[j + 1] {
                // preceding character text[j] is l-type.
                let p = text[j].as_index();
                if suf[p] < EMPTY {
                    // left shift the borrowed chunk.
                    let n = suf[..=p].iter().rev().take_while(|&&x| x <= EMPTY).count();
                    suf.copy_within(p + 1 - n..p + 1, p - n);
                    suf[p] = EMPTY;
                    if i > p - n {
                        // shift cursor.
                        i -= 1;
                    }
                }

                if suf[p] == EMPTY {
                    if p + 1 < suf.len() && suf[p + 1] == EMPTY {
                        suf[p] = from_counter(1);
                        suf[p + 1] = u32::from_index(j);
                    } else {
                        suf[p] = u32::from_index(j);
                    }
                } else {
                    let q = p + 1 + as_counter(suf[p]);
                    if q < suf.len() && suf[q] == EMPTY {
                        suf[p] -= 1;
                        suf[q] = u32::from_index(j);
                    } else {
                        suf.copy_within(p + 1..q, p);
                        suf[q - 1] = u32::from_index(j);
                        if i > p {
                            // shift cursor.
                            i -= 1;
                        }
                    }
                }

                if left_most || clean {
                    // leave lml-characters only if left_most is toggled.
                    suf[i] = EMPTY;
                }
            }
        }
        i += 1;
    }

    // clean up counters.
    for i in 0..suf.len() {
        if suf[i] > EMPTY {
            let n = as_counter(suf[i]);
            suf.copy_within(i + 1..i + 1 + n, i);
            suf[i + n] = EMPTY;
        }
    }
    if left_most {
        // text[0] is not lml-character.
        suf.iter_mut().filter(|p| **p == 0).for_each(|p| *p = EMPTY);
    }
}

/// Induce (left most) s-typed characters from (left most) l-type characters.
#[inline(never)]
fn induce_schars(text: &[u32], suf: &mut [u32], left_most: bool) {
    let mut i = text.len() - 1;
    loop {
        if suf[i] > 0 && suf[i] < EMPTY {
            // text[suf[i]] is non-empty, and has preceding character.
            let j = (suf[i] - 1).as_index();
            if text[j] < text[j + 1] || (text[j] == text[j + 1] && text[j].as_index() > i) {
                // preceding character text[j] is s-type.
                let p = text[j].as_index();
                if suf[p] < EMPTY {
                    // right shift the borrowed chunk.
                    let n = suf[p..].iter().take_while(|&&x| x < EMPTY).count();
                    suf.copy_within(p..p + n, p + 1);
                    suf[p] = EMPTY;
                    if i < p + n {
                        // shift cursor.
                        i += 1;
                    }
                }

                if suf[p] == EMPTY {
                    if p > 0 && suf[p - 1] == EMPTY {
                        suf[p] = from_counter(1);
                        suf[p - 1] = u32::from_index(j);
                    } else {
                        suf[p] = u32::from_index(j);
                    }
                } else {
                    let q = p - as_counter(suf[p]);
                    if q > 0 && suf[q - 1] == EMPTY {
                        suf[p] -= 1;
                        suf[q - 1] = u32::from_index(j);
                    } else {
                        suf.copy_within(q..p, q + 1);
                        suf[q] = u32::from_index(j);
                        if i < p {
                            // shift cursor.
                            i += 1;
                        }
                    }
                }

                if left_most {
                    // leave lms-characters only.
                    suf[i] = EMPTY;
                }
            }
        }

        // avoid underflow.
        if let Some(next_i) = i.checked_sub(1) {
            i = next_i;
        } else {
            break;
        }
    }

    // clean up counters.
    for i in (0..suf.len()).rev() {
        if suf[i] > EMPTY {
            let n = as_counter(suf[i]);
            suf.copy_within(i - n..i, i - n + 1);
            suf[i - n] = EMPTY;
        }
    }
    if left_most {
        // text[0] is not lms-character.
        suf.iter_mut().filter(|p| **p == 0).for_each(|p| *p = EMPTY);
    }
}

#[inline(always)]
fn as_counter(i: u32) -> usize {
    (-(i as i32)) as u32 as usize
}

#[inline(always)]
fn from_counter(n: usize) -> u32 {
    (-(n as u32 as i32)) as u32
}

// tests for sacak32.
#[cfg(test)]
mod tests {
    use super::super::common::saca_tiny;
    use super::super::types::*;
    use super::sacak32;
    use std::collections::BTreeMap;

    #[test]
    fn tablecheck_sacak32() {
        let texts: &[&[u32]] = &[
            &[],
            &[0],
            &[0, 0, 0, 0, 0, 0],
            &[0, 0, 0, 0, 0, 1],
            &[5, 4, 3, 2, 1, 0],
            &[3, 4, 5, 2, 0, 1],
            &[2, 0, 2, 0, 2, 1, 4, 3],
            &[3, 2, 1, 3, 2, 3, 2, 1, 0, 1],
            &[2, 1, 4, 1, 1, 4, 1, 3, 1],
            &[2, 1, 1, 3, 3, 1, 1, 3, 3, 1, 2, 1],
            &[2, 2, 1, 4, 4, 1, 4, 4, 1, 3, 3, 1, 1],
            &[6, 8, 9, 5, 2, 4, 3, 0, 0, 7, 1, 2],
            &[
                1, 2, 2, 1, 1, 0, 0, 1, 1, 2, 2, 0, 0, 2, 2, 0, 1, 0, 2, 0, 1, 1, 1, 1, 2, 2, 0, 0, 2, 1, 2, 1, 1, 0,
                2, 1, 2, 2, 0, 2, 1, 1, 2, 2, 2, 1, 2, 0, 0, 1, 2, 0, 0, 0, 1, 2, 2, 2, 1, 1, 1, 1, 2, 0, 2, 1, 1, 1,
                2, 1, 0, 1,
            ],
        ];

        for &text in texts.iter() {
            assert_eq!(calc_sacak32(text), calc_naive(text));
        }
    }

    #[quickcheck]
    fn quickcheck_sacak32(text: Vec<u32>) -> bool {
        calc_sacak32(&text[..]) == calc_naive(&text[..])
    }

    // helpers.

    fn calc_sacak32(text: &[u32]) -> Vec<u32> {
        let (mut text, k) = reduce_alphabet(text);
        let mut suf = vec![0u32; text.len()];
        sacak32(&mut text[..], &mut suf[..], k);
        suf
    }

    fn calc_naive(text: &[u32]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        saca_tiny(text, &mut suf[..]);
        suf
    }

    fn reduce_alphabet(text: &[u32]) -> (Vec<u32>, usize) {
        let mut dic = BTreeMap::new();
        text.iter().for_each(|&c| {
            dic.insert(c, 0);
        });
        let mut k = 0;
        for (_, v) in &mut dic {
            *v = k;
            k += 1;
        }
        (text.iter().map(|c| *dic.get(c).unwrap() as u32).collect(), k)
    }
}
