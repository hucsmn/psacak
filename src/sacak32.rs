use super::common::*;
use super::types::*;

/// Empty mark in the workspace.
const EMPTY: u32 = 1 << 31;

/// Sort suffix array for integer string (alphabet size k).
pub fn sacak32(text: &mut [u32], suf: &mut [u32], k: usize) {
    debug_assert!(text.len() <= suf.len());
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // make bucket pointers bt rewriting text.
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

    // sort lms-substrings lexicographically in the head of workspace.
    let k = make_subproblem(text, suf, n);
    if k < n {
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
fn transform_text(text: &mut [u32], suf: &mut [u32], k: usize) {
    // make bucket.
    suf[..k].iter_mut().for_each(|p| *p = 0);
    text.iter().for_each(|&c| suf[c.as_index()] += 1);
    suf[..k].iter_mut().fold(0, |sum, p| {
        *p += sum;
        *p
    });

    // l-type => bucket_head, s-type => bucket_tail - 1.
    foreach_typedchars_mut(text, |i, stype, p| {
        let c = p.as_index();
        if !stype {
            *p = if c > 0 { suf[c - 1] } else { 0 };
        } else {
            *p = suf[c] - 1;
        }
    })
}

/// Put lms-characters.
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
fn put_lmssufs(text: &[u32], suf: &mut [u32], n: usize) {
    suf[n..].iter_mut().for_each(|p| *p = EMPTY);

    let mut prev = text.len(); // any integer not in text.
    let mut m = 0; // destination bucket offset counter.
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
            // non-empty, and has preceding character.
            let j = (suf[i] - 1).as_index();
            if text[j] >= text[j + 1] {
                // preceding character is l-type.
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

                if left_most || is_schar(text, suf[i].as_index()) {
                    // clean up lms-characters.
                    // when left_most is toggled, leave lml-characters only.
                    suf[i] = EMPTY;
                }
            }
        }
        i += 1;
    }

    // clean up counters.
    for i in 0..suf.len() {
        if left_most && suf[i] == 0 {
            // text[0] is not lms-character.
            suf[i] = EMPTY;
        }
        if suf[i] > EMPTY {
            let n = as_counter(suf[i]);
            suf.copy_within(i + 1..i + 1 + n, i);
            suf[i + n] = EMPTY;
        }
    }
}

/// Induce (left most) s-typed characters from (left most) l-type characters.
fn induce_schars(text: &[u32], suf: &mut [u32], left_most: bool) {
    let mut i = text.len() - 1;
    loop {
        if suf[i] > 0 && suf[i] < EMPTY {
            // non-empty, and has preceding character.
            let j = (suf[i] - 1).as_index();
            if text[j] <= text[j + 1] {
                // preceding character is s-type.
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
    for i in (1..suf.len()).rev() {
        if left_most && suf[i] == 0 {
            // text[0] is not lms-character.
            suf[i] = EMPTY;
        }
        if suf[i] > EMPTY {
            let n = as_counter(suf[i]);
            suf.copy_within(i - n..i, i - n + 1);
            suf[i - n] = EMPTY;
        }
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
#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use super::super::types::*;
    use super::super::common::saca_tiny;
    use super::sacak32;

    #[test]
    fn tablecheck_sacak32() {
        let texts: &[&[u32]] = &[
            &[],
            &[0],
            &[0, 0, 0, 0, 0, 0],
            &[2, 0, 2, 0, 2, 1, 4, 3],
            &[3, 2, 1, 3, 2, 3, 2, 1, 0, 1],
            &[2, 1, 4, 1, 1, 4, 1, 3, 1],
            &[2, 1, 1, 3, 3, 1, 1, 3, 3, 1, 2, 1],
            &[2, 2, 1, 4, 4, 1, 4, 4, 1, 3, 3, 1, 1],
            &[
                1, 2, 2, 1, 1, 0, 0, 1, 1, 2, 2, 0, 0, 2, 2, 0, 1, 0, 2, 0, 1, 1, 1, 1, 2, 2, 0, 0, 2,
                1, 2, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 1, 2, 2, 2, 1, 2, 0, 0, 1, 2, 0, 0, 0, 1, 2, 2, 2,
                1, 1, 1, 1, 2, 0, 2, 1, 1, 1, 2, 1, 0, 1,
            ],
        ];

        for &text in texts.iter() {
            assert_eq!(sacak(text), naive(text));
        }
    }

    quickcheck! {
        fn quickcheck_sacak8(text: Vec<u32>) -> bool {
            naive(&text[..]) == sacak(&text[..])
        }
    }

    fn sacak(text: &[u32]) -> Vec<u32> {
        let mut dic = BTreeMap::new();
        text.iter().for_each(|&c| { dic.insert(c, 0); });

        let mut k = 1;
        for (_, v) in &mut dic {
            *v = k - 1;
            k += 1;
        }
        let mut text: Vec<_> = text.iter().map(|c| *dic.get(c).unwrap() as u32).collect();

        let mut suf = vec![0u32; text.len()];
        sacak32(&mut text[..], &mut suf[..], k);
        suf
    }

    fn naive(text: &[u32]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        saca_tiny(text, &mut suf[..]);
        suf
    }
}
