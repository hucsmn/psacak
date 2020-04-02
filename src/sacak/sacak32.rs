use super::common::*;
use super::naming::*;
use super::types::*;

/// Empty symbol in workspace.
const EMPTY: u32 = 1 << 31;

/// Inner level of SACA-K for integer strings.
///
/// Assumes that characters are correctly translated to bucket pointers.
#[inline]
pub fn sacak32(text: &[u32], suf: &mut [u32]) {
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // induce sort lms-substrings.
    put_lmscharacters(text, suf);
    induce_sort(text, suf, true);

    // construct subproblem, compute its suffix array, and get sorted lms-suffixes.
    let n = compact_include(suf, 1..EMPTY, false);
    let k = name_lmssubstrings(text, suf, n);
    if k < n {
        // need to solve the subproblem recursively.
        let (suf1, text1) = suf.split_at_mut(suf.len() - n);
        sacak32(text1, suf1);
    } else {
        // the subproblem itself is the inversed suffix array.
        let (suf1, text1) = suf.split_at_mut(suf.len() - n);
        for i in 0..n {
            suf1[text1[i].as_index()] = i as u32;
        }
    }
    permutate_lmssuffixes(text, suf, n);

    // induce sort the suffix array from sorted lms-suffixes.
    put_lmssuffixes(text, suf, n);
    induce_sort(text, suf, false);
}

/// Put lms-characters to their corresponding bucket tails, in arbitary order.
#[inline]
fn put_lmscharacters(text: &[u32], suf: &mut [u32]) {
    suf.iter_mut().for_each(|p| *p = EMPTY);

    foreach_lmschars(text, |i, c| {
        let p = c as usize;
        if suf[p] < EMPTY {
            // right shift the borrowed chunk.
            let n = suf[p..].iter().take_while(|&&x| x < EMPTY).count();
            suf.copy_within(p..p + n, p + 1);
            suf[p] = EMPTY;
        }

        if suf[p] == EMPTY {
            if p > 0 && suf[p - 1] == EMPTY {
                suf[p] = to_counter(1);
                suf[p - 1] = i as u32;
            } else {
                suf[p] = i as u32;
            }
        } else {
            let q = p - get_counter(suf[p]);
            if q > 0 && suf[q - 1] == EMPTY {
                suf[p] -= 1;
                suf[q - 1] = i as u32;
            } else {
                suf.copy_within(q..p, q + 1);
                suf[q] = i as u32;
            }
        }
    });

    // clean up counters.
    for i in (1..suf.len()).rev() {
        if suf[i] > EMPTY {
            let n = get_counter(suf[i]);
            suf.copy_within(i - n..i, i - n + 1);
            suf[i - n] = EMPTY;
        }
    }
}

/// Put the sorted lms-suffixes, originally located in head of workspace, to their corresponding bucket tails.
#[inline]
fn put_lmssuffixes(text: &[u32], suf: &mut [u32], n: usize) {
    suf[n..].iter_mut().for_each(|p| *p = EMPTY);

    // copy sorted lms-suffixes bucket by bucket.
    let mut prev = text.len();
    let mut m = 0;
    for i in (0..n).rev() {
        let x = suf[i] as usize;
        suf[i] = EMPTY;

        let p = text[x] as usize;
        if p != prev {
            m = 0;
        }
        suf[p - m] = x as u32;
        prev = p;
        m += 1;
    }
}

/// Induce sort all the suffixes (or lms-substrings) from the sorted lms-suffixes (or lms-characters).
#[inline]
fn induce_sort(text: &[u32], suf: &mut [u32], left_most: bool) {
    // stage 1. induce l (or lml) from lms.

    // the sentinel.
    let p = text[text.len() - 1] as usize;
    if p + 1 < suf.len() && suf[p + 1] == EMPTY {
        suf[p] = to_counter(1);
        suf[p + 1] = (text.len() - 1) as u32;
    } else {
        suf[p] = (text.len() - 1) as u32;
    }

    let mut i = 0;
    while i < suf.len() {
        if suf[i] > 0 && suf[i] < EMPTY {
            // c1 is non-empty, and has a preceding character c0.
            let j = (suf[i] - 1) as usize;
            let c0 = text[j];
            let c1 = text[j + 1];
            let p = c0 as usize;
            if c0 >= c1 {
                // c0 is l-type.
                if suf[p] < EMPTY {
                    // left shift the borrowed chunk.
                    let n = suf[..p + 1].iter().rev().take_while(|&&x| x <= EMPTY).count();
                    suf.copy_within(p + 1 - n..p + 1, p - n);
                    suf[p] = EMPTY;
                    if i > p - n {
                        // shift cursor.
                        i -= 1;
                    }
                }

                if suf[p] == EMPTY {
                    if p + 1 < suf.len() && suf[p + 1] == EMPTY {
                        suf[p] = to_counter(1);
                        suf[p + 1] = j as u32;
                    } else {
                        suf[p] = j as u32;
                    }
                } else {
                    let q = p + 1 + get_counter(suf[p]);
                    if q < suf.len() && suf[q] == EMPTY {
                        suf[p] -= 1;
                        suf[q] = j as u32;
                    } else {
                        suf.copy_within(p + 1..q, p);
                        suf[q - 1] = j as u32;
                        if i > p {
                            // shift cursor.
                            i -= 1;
                        }
                    }
                }

                // manually clean up lms-suffixes.
                let mut ltype = true;
                if j + 2 < text.len() {
                    let c2 = text[j + 2];
                    ltype = c1 > c2 || (c1 == c2 && i > c1 as usize);
                }
                if left_most || !ltype {
                    // only keep lml-suffixes if toggled left_most.
                    suf[i] = EMPTY;
                }
            }
        }
        i += 1;
    }

    // clean up bucket counters.
    for i in 0..suf.len() {
        if suf[i] > EMPTY {
            let n = get_counter(suf[i]);
            suf.copy_within(i + 1..i + 1 + n, i);
            suf[i + n] = EMPTY;
        }
    }

    // stage 2. induce s (or lms) from l (or lml).

    let mut i = text.len() - 1;
    loop {
        if suf[i] > 0 && suf[i] < EMPTY {
            // c1 is non-empty, and has a preceding character c0.
            let j = (suf[i] - 1) as usize;
            let c0 = text[j];
            let c1 = text[j + 1];
            let p = c0 as usize;
            if c0 < c1 || (c0 == c1 && i < p) {
                // c0 is s-type.
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
                        suf[p] = to_counter(1);
                        suf[p - 1] = j as u32;
                    } else {
                        suf[p] = j as u32;
                    }
                } else {
                    let q = p - get_counter(suf[p]);
                    if q > 0 && suf[q - 1] == EMPTY {
                        suf[p] -= 1;
                        suf[q - 1] = j as u32;
                    } else {
                        suf.copy_within(q..p, q + 1);
                        suf[q] = j as u32;
                        if i < p {
                            // shift cursor.
                            i += 1;
                        }
                    }
                }

                if left_most {
                    // only keep lms-suffixes if toggled left_most.
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
}

#[inline(always)]
fn get_counter(i: u32) -> usize {
    (-(i as i32)) as u32 as usize
}

#[inline(always)]
fn to_counter(n: usize) -> u32 {
    (-(n as u32 as i32)) as u32
}

// Simple sacak32 tests.
#[cfg(test)]
mod tests {
    use super::super::common::*;
    use super::super::types::*;
    use super::sacak32;
    use std::collections::BTreeMap;

    #[test]
    fn tablecheck_sacak32() {
        let texts: &[&[u32]] = &[
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

    // helper functions.

    fn calc_sacak32(text: &[u32]) -> Vec<u32> {
        let (mut text, k) = reduce_alphabet(text);
        let mut suf = vec![0u32; text.len()];
        if k >= text.len() {
            // sacak32 cannot sort permutations, simply skip those test cases.
            saca_tiny(&text[..], &mut suf[..]);
        } else {
            translate_text(&mut text[..], &mut suf[..], k);
            sacak32(&text[..], &mut suf[..]);
        }
        suf
    }

    fn calc_naive(text: &[u32]) -> Vec<u32> {
        let mut suf = vec![0u32; text.len()];
        saca_tiny(text, &mut suf[..]);
        suf
    }

    fn reduce_alphabet(text: &[u32]) -> (Vec<u32>, usize) {
        let mut ranking = BTreeMap::new();
        for c in text.iter().cloned() {
            ranking.insert(c, 0);
        }
        let mut k = 0;
        for (_, rank) in &mut ranking {
            *rank = k;
            k += 1;
        }
        (text.iter().map(|c| *ranking.get(c).unwrap() as u32).collect(), k)
    }

    fn translate_text(text: &mut [u32], suf: &mut [u32], k: usize) {
        suf[..k + 1].iter_mut().for_each(|p| *p = 0);
        for c in text.iter().cloned() {
            suf[c as usize + 1] += 1;
        }
        let mut p = 0;
        for i in 1..k + 1 {
            let cnt = suf[i];
            suf[i] += p;
            p += cnt;
        }

        foreach_typedchars_mut(text, |i, t, p| {
            let c = *p as usize;
            if !t.stype {
                *p = suf[c];
            } else if t.stype {
                *p = suf[c + 1] - 1;
            }
        })
    }
}
