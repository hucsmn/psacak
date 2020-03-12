#![allow(unused)]

#[cfg(test)]
extern crate quickcheck;
#[cfg(test)]
#[macro_use(quickcheck)]
extern crate quickcheck_macros;

mod common;
mod pipeline;
mod sacak32;
mod sacak8;
mod types;

use std::env;
use std::fs;
use std::time;
use std::path;
use std::io;

use byteorder::{NativeEndian, WriteBytesExt};

fn main() {
    for inname in env::args().skip(1) {
        eprintln!("* {} *", inname);
        let data = fs::read(&inname).unwrap();
        eprintln!("load file `{}` of {} bytes", inname, data.len());

        let mut suf = vec![0u32; data.len()];
        let ((), t) = timeit(|| sacak8::sacak8(&data[..], &mut suf[..]));
        eprintln!("construct suffix array in {:.3}s", t.as_secs_f64());

        let mut outname = path::PathBuf::from(inname);
        outname.set_extension("suffix");
        eprintln!("store file `{:?}` of {} bytes", outname.as_path(), 4 * suf.len());
        let mut file = io::BufWriter::new(fs::File::create(outname.as_path()).unwrap());
        suf.iter().cloned().for_each(|x| file.write_u32::<NativeEndian>(x).unwrap());
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
