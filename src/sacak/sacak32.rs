use super::common::*;
use super::ranking::*;
use super::types::*;

/// Empty symbol in workspace.
const EMPTY: u32 = 1 << 31;

/// SACA-K inner level sort algorithm for integer strings (alphabet size k).
///
/// Assumes that for each c in text, c < k.
#[inline]
pub fn sacak32(text: &mut [u32], suf: &mut [u32], k: usize) {
    let suf = &mut suf[..text.len()];

    if text.len() <= 3 {
        saca_tiny(text, suf);
        return;
    }

    // make buckets by rewriting text.
    make_buckets(text, suf, k);

    // induce sort lms-substrings.
    put_lmschars(text, suf);
    induce_lchars(text, suf, true);
    induce_schars(text, suf, true);
    let n = compact_lmssubs(text, suf);

    // get ranks of lms-substrings into the tail of workspace.
    let k = rank_lmssubs(text, suf, n);

    if k < n {
        // order of lms-suffixes != order of lms-substrings.
        // need to further sort the lms-suffixes.
        {
            let (subsuf, subtext) = suf.split_at_mut(suf.len() - n);
            sacak32(subtext, subsuf, k);
        }
        unrank_lmssufs(text, suf, n);
    }

    // induce sort the suffix array from sorted lms-suffixes.
    put_lmssufs(text, suf, n);
    induce_lchars(text, suf, false);
    induce_schars(text, suf, false);
}

/// Rewrite characters to their corresponding bucket pointers.
///
/// The lexicographical order would keep unchanged.
#[inline]
fn make_buckets(text: &mut [u32], suf: &mut [u32], k: usize) {
    // calculate bucket pointers.
    suf[..k + 1].iter_mut().for_each(|p| *p = 0);
    text.iter().for_each(|&c| suf[c as usize + 1] += 1);
    suf[1..k + 1].iter_mut().fold(0, |sum, p| {
        *p += sum;
        *p
    });

    // translate characters to bucket pointers.
    foreach_typedchars_mut(text, |i, t, p| {
        let c = *p as usize;
        if !t.stype {
            // l-type => bucket head.
            *p = suf[c];
        } else if t.stype {
            // s-type => bucket tail - 1.
            *p = suf[c + 1] - 1;
        }
    })
}

/// Put lms-characters to their corresponding bucket tails, in arbitary order.
#[inline]
fn put_lmschars(text: &[u32], suf: &mut [u32]) {
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

/// Compact the sorted lms-substrings into the head of workspace.
#[inline]
fn compact_lmssubs(text: &[u32], suf: &mut [u32]) -> usize {
    let mut n = 0;
    for i in 0..suf.len() {
        if suf[i] < EMPTY {
            suf[n] = suf[i];
            n += 1;
        }
    }
    n
}

/// Put the sorted lms-suffixes, originally located in the head of workspace,
/// to their corresponding bucket tails.
#[inline]
fn put_lmssufs(text: &[u32], suf: &mut [u32], n: usize) {
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

/// Induce l-suffixes (or lml-suffixes) from sorted lms-suffixes.
///
/// Assumes that non lms-suffixes among the input `suf` have been reset as `EMPTY`.
///
/// Outputs the induced l-suffixes with the input lms-suffixes reset to `EMPTY`,
/// or the induced lml-suffixes with all other suffixes reset to `EMPTY`.
#[inline]
fn induce_lchars(text: &[u32], suf: &mut [u32], left_most: bool) {
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

                // manually clean up lms.
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
    if left_most {
        // text[0..] is not a lml-suffix.
        suf.iter_mut().filter(|p| **p == 0).for_each(|p| *p = EMPTY);
    }
}

/// Induce s-suffixes (or lms-suffixes) from sorted l-suffixes (or lml-suffixes).
///
/// Assumes that non l-suffixes (or non lml-suffixes) among the input `suf`
/// have been reset as `EMPTY`.
///
/// Outputs the induced s-suffixes together with the input l-suffixes,
/// or the induced lms-suffixes with all other suffixes reset to `EMPTY`.
#[inline]
fn induce_schars(text: &[u32], suf: &mut [u32], left_most: bool) {
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

    // clean up bucket counters.
    for i in (0..suf.len()).rev() {
        if suf[i] > EMPTY {
            let n = get_counter(suf[i]);
            suf.copy_within(i - n..i, i - n + 1);
            suf[i - n] = EMPTY;
        }
    }
    if left_most {
        // text[0..] is not a lms-suffix.
        suf.iter_mut().filter(|p| **p == 0).for_each(|p| *p = EMPTY);
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

    // helper funtions.

    fn calc_sacak32(text: &[u32]) -> Vec<u32> {
        let (mut text, k) = reduce_alphabet(text);
        let mut suf = vec![0u32; text.len()];
        if k >= text.len() {
            saca_tiny(&text[..], &mut suf[..]);
        } else {
            sacak32(&mut text[..], &mut suf[..], k);
        }
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
