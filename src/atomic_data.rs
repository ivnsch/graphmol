use std::collections::HashMap;

pub struct AtomicData {
    pub anum: AtomicNumber,
    pub symb: String,
    pub name: String,
    pub valence: Vec<i32>,
    pub mass: f64,
    pub n_val: u8,
    pub common_isotope: u8,
    pub common_isotope_mass: f64,
    pub row: u32,
    // port: TODO populate
    pub isotope_info_map: HashMap<u32, (f64, f64)>,
}

impl AtomicData {
    pub fn default_valence(&self) -> i32 {
        // port: TODO review unwrap: consider using data structure that has at least one element
        *self.valence.first().unwrap()
    }
}

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct AtomicNumber(pub usize);
