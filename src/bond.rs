use crate::{atom::Atom, mol::Mol};

#[derive(PartialEq, Eq, Debug)]
pub enum BondType {
    Unspecified,
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Hextuple,
    OneAndAHalf,
    TwoAndAHalf,
    ThreeAndAHalf,
    FourAndAHalf,
    FiveAndAHalf,
    Aromatic,
    Ionic,
    Hydrogen,
    ThreeCenter,
    DativeOne, // one-electron dative (e.g. from a c in a cp ring to a metal)
    Dative,    // standard two-electron dative
    DativeL,   // standard two-electron dative
    DativeR,   // standard two-electron dative
    Other,
    Zero, // Zero-order bond (from
          // http://pubs.acs.org/doi/abs/10.1021/ci200488k)
}

#[derive(PartialEq, Eq, Debug)]
/// the bond's direction (for chirality)
pub enum BondDir {
    None = 0,   // no special style
    BeginWedge, // wedged: narrow at begin
    BeginDash,  // dashed: narrow at begin
    // fix: this may not really be adequate
    EndDownRight, // for cis/trans
    EndUpRight,   //  ditto
    EitherDouble, // a "crossed" double bond
    Unknown,      // intentionally unspecified stereochemistry
}

/// the nature of the bond's stereochem (for cis/trans)
#[derive(PartialEq, Eq, Debug)]
pub enum BondStereo {
    // stereochemistry of double bonds
    StereOnOne = 0, // no special style
    StereoAny,      // intentionally unspecified
    // -- put any true specifications about this point so
    // that we can do comparisons like if(bond->getstereo()>bond::stereoany)
    StereoZ,        // z double bond
    StereoE,        // e double bond
    StereoCis,      // cis double bond
    StereoTrans,    // trans double bond
    StereoAtropCW,  //  atropisomer clockwise rotation
    StereoAtropCCW, //  atropisomer counter clockwise rotation
}

#[derive(Debug, PartialEq, Eq)]
pub struct Bond {
    pub stereo_atoms: Vec<i32>,
    pub index: usize,
    pub begin_atom_idx: usize,
    pub end_atom_idx: usize,
    pub is_aromatic: bool,
    pub is_conjugated: bool,
    pub bond_type: BondType,
    pub dir_tag: BondDir,
    pub stereo: BondStereo,
}

#[derive(Debug)]
pub struct BondWithMol {
    pub bond: Bond,
    pub mol: Option<Mol>,
}

impl Bond {
    fn is_dative(&self) -> bool {
        let t = &self.bond_type;
        *t == BondType::Dative
            || *t == BondType::DativeL
            || *t == BondType::DativeR
            || *t == BondType::DativeOne
    }

    fn can_set_double_bond_stereo(bond: Bond) -> bool {
        let t = &bond.bond_type;
        *t == BondType::Single || *t == BondType::Aromatic || bond.is_dative()
    }

    fn can_have_direction(bond: Bond) -> bool {
        let t = &bond.bond_type;
        *t == BondType::Single || *t == BondType::Aromatic
    }

    fn get_bond_type_as_float(&self) -> f32 {
        match self.bond_type {
            BondType::Unspecified | BondType::Ionic | BondType::Zero => 0.0,
            BondType::Single => 1.0,
            BondType::Double => 2.0,
            BondType::Triple => 3.0,
            BondType::Quadruple => 4.0,
            BondType::Quintuple => 5.0,
            BondType::Hextuple => 6.0,
            BondType::OneAndAHalf => 1.5,
            BondType::TwoAndAHalf => 2.5,
            BondType::ThreeAndAHalf => 3.5,
            BondType::FourAndAHalf => 4.5,
            BondType::FiveAndAHalf => 5.5,
            BondType::Aromatic => 1.5,
            BondType::Hydrogen => 0.0,
            _ => panic!("Bad bond type"),
        }
    }
    // This method can be used to distinguish query bonds from standard bonds
    pub fn has_query(&self) -> bool {
        false
    }
}

impl BondWithMol {
    fn get_begin_atom(&self) -> Option<&Atom> {
        self.mol.as_ref()?.atom_with_idx(self.bond.begin_atom_idx)
    }

    fn get_end_atom(&self) -> Option<&Atom> {
        self.mol.as_ref()?.atom_with_idx(self.bond.end_atom_idx)
    }

    fn get_other_atom(&self) -> Option<&Atom> {
        self.mol.as_ref()?.atom_with_idx(self.bond.index)
    }

    pub fn get_valence_contrib(&self, atom: &Atom) -> f64 {
        if Some(atom) != self.get_begin_atom() && Some(atom) != self.get_end_atom() {
            0.
        } else if (self.bond.bond_type == BondType::Dative
            || self.bond.bond_type == BondType::DativeOne)
            && self.bond.index != self.bond.end_atom_idx
        {
            0.
        } else {
            self.bond.get_bond_type_as_float() as f64
        }
    }
}
