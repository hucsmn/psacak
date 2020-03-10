use super::types::*;

pub struct InducePipeline<'a, C: SacaChar, I: SacaIndex> {
    forward: bool,
    chunk: usize,
    text: &'a [C],
    suf: &'a mut [I],
}

pub struct InduceContext<'a, C: SacaChar, I: SacaIndex> {
    pipeline: InducePipeline<'a, C, I>,
    rcache: Vec<(usize, C)>,
    wcache: Vec<(usize, I)>,
}
