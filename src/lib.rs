#[cfg(test)]
#[macro_use]
extern crate quickcheck_macros;

#[macro_use]
extern crate cfg_if;

mod common;
mod naming;
mod pipeline;
mod sacak32;
mod sacak8;
mod types;

/// The pSACAK algorithm.
pub fn psacak(text: &[u8], suf: &mut [u32]) {
    assert!(text.len() <= suf.len());
    let suf = &mut suf[..text.len()];
    sacak8::sacak8(text, suf);
}
