use crate::atomic_data::{AtomicData, AtomicNumber};

pub struct PeriodicTable {
    by_anum: Vec<AtomicData>,
}

impl PeriodicTable {
    // port: TODO singleton
    pub fn new() -> PeriodicTable {
        // port: TODO fill with data
        PeriodicTable { by_anum: vec![] }
    }

    pub fn max_atom_num() -> AtomicNumber {
        AtomicNumber(103)
    }
}

impl PeriodicTable {
    pub fn valence_list(&self, atomic_number: AtomicNumber) -> Vec<i32> {
        // port: TODO clone - performance
        self.by_anum[atomic_number.0].valence.clone()
    }

    pub fn default_valence(&self, atomic_number: AtomicNumber) -> i32 {
        self.by_anum[atomic_number.0].default_valence()
    }

    pub fn element_symbol(&self, atomic_number: AtomicNumber) -> String {
        self.by_anum[atomic_number.0].symb.clone()
    }

    pub fn mass_for_isotope(&self, atomic_number: AtomicNumber, isotope: u32) -> f64 {
        let info_map = &self.by_anum[atomic_number.0].isotope_info_map;
        info_map
            .get(&isotope)
            .map(|value| value.0)
            .unwrap_or_else(|| 0.0)
    }

    pub fn atomic_weight(&self, atomic_number: AtomicNumber) -> f64 {
        self.by_anum[atomic_number.0].mass
    }

    pub fn get_row(&self, atomic_number: AtomicNumber) -> u32 {
        self.by_anum[atomic_number.0].row
    }
}
