#![allow(unused)]

mod common;
mod pipeline;
mod sacak32;
mod sacak8;
mod types;

fn sacak_bytes(text: &[u8]) -> Vec<u32> {
    let mut suf = vec![0u32; text.len() + 1];
    suf[0] = text.len() as u32;
    sacak8::sacak8(text, &mut suf[1..]);
    suf
}

fn sacak_uints(text: &[u32]) -> Vec<u32> {
    let k = text.iter().map(|&x| x as usize).max().unwrap_or(0) + 1;
    let mut text = Vec::from(text);
    let mut suf = vec![0u32; text.len() + 1];
    suf[0] = text.len() as u32;
    sacak32::sacak32(&mut text[..], &mut suf[1..], k);
    suf
}

fn naive_saca<C: types::SacaChar>(text: &[C]) -> Vec<u32> {
    let mut suf = vec![0u32; text.len() + 1];
    suf[0] = text.len() as u32;
    common::saca_tiny(text, &mut suf[1..]);
    suf
}

fn main() {
    /*
    let texts: &[&[u8]] = &[
        &[],
        &[0],
        &[0, 0, 0, 0, 0, 0],
        &[2, 0, 2, 0, 2, 1, 4, 3],
        &[3, 2, 1, 3, 2, 3, 2, 1, 0, 1],
        &[2, 1, 4, 1, 1, 4, 1, 3, 1],
        &[2, 1, 1, 3, 3, 1, 1, 3, 3, 1, 2, 1],
        &[2, 2, 1, 4, 4, 1, 4, 4, 1, 3, 3, 1, 1],
        &[
            1, 2, 2, 1, 1, 0, 0, 1, 1, 2, 2, 0, 0, 2, 2, 0, 1, 0, 2, 0, 1, 1, 1, 1, 2, 2, 0, 0, 2,
            1, 2, 1, 1, 0, 2, 1, 2, 2, 0, 2, 1, 1, 2, 2, 2, 1, 2, 0, 0, 1, 2, 0, 0, 0, 1, 2, 2, 2,
            1, 1, 1, 1, 2, 0, 2, 1, 1, 1, 2, 1, 0, 1,
        ],
    ];
    for &text in texts.iter() {
        eprintln!("text: {:?}", text);
        assert_eq!(sacak_bytes(text), naive_saca(text));
    }
    */

    sacak_bytes(&[2, 1, 1, 3, 3, 1, 1, 3, 3, 1, 2, 1]);
}
