use indexmap::IndexMap;

use super::recmap::RecMap;

pub type GenomeCoord = i32;

pub struct Genome {
    pub seqlens: IndexMap<String, GenomeCoord>,
    pub recmap: RecMap,
}

