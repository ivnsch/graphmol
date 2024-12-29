use crate::{atom::Atom, bond::BondWithMol};

#[derive(Debug)]
struct QueryBond {
    bond: BondWithMol,
}

impl QueryBond {
    fn get_valence_contrib(&self, atom: &Atom) -> f64 {
        // port: for now ignoring rest of implementation
        self.bond.get_valence_contrib(atom)
    }
}
