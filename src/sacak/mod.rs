mod common;
mod naming;
mod pipeline;
mod sacak32;
mod sacak8;
mod types;

/// The SACA-K algorithm.
pub fn sacak(text: &[u8], suf: &mut [u32]) {
    assert!(text.len() <= suf.len());
    let suf = &mut suf[..text.len()];
    sacak8::sacak8(text, suf);
}
