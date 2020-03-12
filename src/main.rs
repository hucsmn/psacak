#![allow(unused)]

#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;

mod sacak;

use std::env;
use std::fs;
use std::io;
use std::time;

fn main() {
    for name in env::args().skip(1) {
        eprintln!("* {} *", name);
        let text = fs::read(&name).unwrap();
        eprintln!("load: {} bytes", text.len());
        let mut suf = vec![0; text.len()];
        let ((), dur) = timeit(|| sacak::sacak(&text[..], &mut suf[..]));
        eprintln!("saca-k: runs in {:.3}s", dur.as_secs_f64());
        eprintln!("validate: {}", validate(&text[..], &suf[..]));
        eprintln!("");
    }
}

fn timeit<F, T>(f: F) -> (T, time::Duration)
where
    F: FnOnce() -> T,
{
    let start = time::Instant::now();
    let ret = f();
    let dur = start.elapsed();
    (ret, dur)
}

fn validate(text: &[u8], suf: &[u32]) -> bool {
    for i in 1..suf.len() {
        if text[suf[i - 1] as usize..] >= text[suf[i] as usize..] {
            return false;
        }
    }
    true
}
