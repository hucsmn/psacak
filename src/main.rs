#![allow(unused)]

mod common;
mod pipeline;
mod sacak32;
mod sacak8;
mod types;

fn sacak(text: &[u8]) -> Vec<u32> {
    let mut suf = vec![0u32; text.len() + 1];
    suf[0] = text.len() as u32;
    sacak8::sacak8(text, &mut suf[1..]);
    suf
}

fn naive(text: &[u8]) -> Vec<u32> {
    let mut suf = vec![0u32; text.len() + 1];
    suf[0] = text.len() as u32;
    common::saca_tiny(text, &mut suf[1..]);
    suf
}

fn main() {
    let text = &[2, 0, 2, 0, 2, 1, 4, 3];
    println!("sacak: {:?}", sacak(text));
    println!("naive: {:?}", naive(text));
}
