use super::common::*;
use super::types::*;

pub fn sacak32(text: &mut [u32], suf: &mut [u32]) {
    debug_assert!(text.len() <= suf.len());
    let suf = &mut suf[..text.len()];

    saca_tiny(text, suf)
}
