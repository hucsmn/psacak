#[macro_use]
extern crate clap;

use std::fs;
use std::mem;
use std::process;
use std::time;

use libc::{getrusage, rusage, RUSAGE_SELF};
use psacak::psacak_inplace;

fn main() {
    let matches = clap_app!(psacak =>
        (about: "pSACAK suffix sorting test")
        (@arg TIMES: -t --times +takes_value "repeat multiple times")
        (@arg CHECK: -c --check "check the suffix array")
        (@arg INPUT: +required "the data to sort")
    )
    .get_matches();

    let input_file = matches.value_of("INPUT").unwrap();
    let check_suffix_array = matches.is_present("CHECK");
    let repeat_times = Ord::max(
        matches
            .value_of("TIMES")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(1),
        1,
    );

    let text = match fs::read(&input_file) {
        Ok(text) => text,
        Err(err) => {
            eprintln!("error: {:?}", err);
            process::exit(1);
        }
    };
    eprintln!("load {} bytes from `{}`", text.len(), &input_file);

    let mut suf = vec![0; text.len()];
    let mut times = String::new();
    for _ in 0..repeat_times {
        let ((), dur) = timeit(|| psacak_inplace(&text[..], &mut suf[..]));
        times.push_str(format!("{:.3}s ", dur.as_secs_f64()).as_str());
    }
    eprintln!(" time: {}", times);

    eprintln!("  rss: {:.3}MiB", get_peak_rss_kib() as f64 / 1024.0);
    if check_suffix_array {
        eprintln!("check: {}", check(&text[..], &suf[..]));
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

fn get_peak_rss_kib() -> u64 {
    let mut ru;
    unsafe {
        ru = mem::zeroed::<rusage>();
        getrusage(RUSAGE_SELF, &mut ru as *mut rusage);
    }
    ru.ru_maxrss as u64
}

fn check(text: &[u8], suf: &[u32]) -> bool {
    for i in 1..suf.len() {
        if text[suf[i - 1] as usize..] >= text[suf[i] as usize..] {
            return false;
        }
    }
    true
}
