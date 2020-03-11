#![allow(unused)]

#[cfg(test)]
#[macro_use]
extern crate quickcheck;

mod common;
mod pipeline;
mod sacak32;
mod sacak8;
mod types;

use std::env::args;
use std::fs;
use std::iter::FromIterator;
use std::time;

fn timeit<F, T>(f: F) -> (T, time::Duration)
where
    F: FnOnce() -> T,
{
    let start = time::Instant::now();
    let ret = f();
    let dur = start.elapsed();
    (ret, dur)
}

fn main() {
    for path in args().skip(1) {
        let data = fs::read(&path).unwrap();
        println!("* file `{}` ({} bytes) *", &path, data.len());

        {
            let mut suf = vec![0; data.len()];
            let ((), t) = timeit(|| sacak8::sacak8(&data[..], &mut suf[..]));
            let sum = suf.iter().fold(0u32, |sum, &x| sum.wrapping_add(x));
            println!("sacak {:.3}s (checksum: 0x{:x})", t.as_secs_f64(), sum);
        }
    }
}
