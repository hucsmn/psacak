//! The pSACAK suffix sorting algorithm for huge in-memory data on multicore machine.
//!
//! The algorithm is described in this [paper](https://doi.org/10.1109/TC.2018.2842050).
//!
//! # Examples
//!
//! Create new suffix array.
//!
//! ```rust
//! use psacak::psacak;
//!
//! let text = b"mississippi";
//! let suf = psacak(text);
//!
//! assert_eq!(suf, vec![10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]);
//! ```
//!
//! Compute the suffix array in a slice.
//!
//! ```rust
//! use psacak::psacak_inplace;
//!
//! let text = b"mississippi";
//! let mut suf = vec![0; text.len()];
//! psacak_inplace(text, &mut suf[..]);
//!
//! assert_eq!(suf, vec![10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]);
//! ```
//!
//! Note that the position of sentinel is not presented in the suffix array calculated by pSACAK.
//!
//! To obtain a suffix array which included the position of sentinel, you can try this:
//!
//! ```rust
//! use psacak::psacak_inplace;
//!
//! let text = b"mississippi";
//! let mut suf = vec![0; text.len() + 1];
//! suf[0] = text.len() as u32;
//! psacak_inplace(text, &mut suf[1..]);
//!
//! assert_eq!(suf, vec![11, 10, 7, 4, 1, 0, 9, 8, 6, 3, 5, 2]);
//! ```

#[cfg(test)]
#[macro_use]
extern crate quickcheck_macros;

#[macro_use]
extern crate cfg_if;

mod common;
mod naming;
mod pipeline;
mod psacak32;
mod psacak8;
mod types;

use psacak8::psacak8;

/// The pSACAK suffix sorting algorithm.
///
/// It would panic if the length of text has exceeded the limitation of 32-bit integer.
pub fn psacak(text: &[u8]) -> Vec<u32> {
    assert!(text.len() < std::u32::MAX as usize);
    let mut suf = vec![0; text.len()];
    psacak8(text, &mut suf[..]);
    suf
}

/// In-place version of the pSACAK suffix sorting algorithm.
///
/// If any of the following conditions is met, it would panic:
/// * The length of text has exceeded the limitation of 32-bit integer.
/// * The length of text is greater than the length of suffix array.
/// * Memory layout of the suffix array slice is not aligned.
pub fn psacak_inplace(text: &[u8], suf: &mut [u32]) {
    assert!(text.len() <= suf.len());
    assert!(text.len() < std::u32::MAX as usize);
    let suf = &mut suf[..text.len()];
    psacak8(text, suf);
}
