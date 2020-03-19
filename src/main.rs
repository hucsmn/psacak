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
use std::mem;
use std::process;
use std::time;

use libc::{getrusage, rusage, RUSAGE_SELF};

const TIMES: usize = 1;

fn main() {
    let text;
    if let Some(name) = env::args().skip(1).next() {
        text = fs::read(&name).unwrap();
        eprintln!("load {} bytes from `{}`", text.len(), &name);
    } else {
        eprintln!("error: require a input file");
        process::exit(1);
    }

    let mut suf = vec![0; text.len()];
    let mut times = String::new();
    for _ in 0..TIMES {
        let ((), dur) = timeit(|| sacak::sacak(&text[..], &mut suf[..]));
        times.push_str(format!("{:.3}s ", dur.as_secs_f64()).as_str());
    }
    eprintln!(" time: {}", times);

    eprintln!("  rss: {:.3}MiB", peak_rss_kib() as f64 / 1024.0);
    eprintln!("check: {}", validate(&text[..], &suf[..]));
}

fn peak_rss_kib() -> u64 {
    let mut ru;
    unsafe {
        ru = mem::zeroed::<rusage>();
        getrusage(RUSAGE_SELF, &mut ru as *mut rusage);
    }
    ru.ru_maxrss as u64
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
