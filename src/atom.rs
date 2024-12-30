use crate::{
    atomic_data::AtomicNumber,
    bond::BondType,
    common_props::{common_props, CommonPropsKey, CommonPropsValue},
    mol::Mol,
    periodic_table::PeriodicTable,
};
use anyhow::{anyhow, Result};

enum HybridizationType {
    Unspecified, // hybridization that hasn't been specified
    S,
    Sp,
    Sp2,
    Sp3,
    Sp2d,
    Sp3d,
    Sp3d2,
    Other,
}

#[derive(Debug, PartialEq, Eq)]
enum ChiralType {
    ChiUnspecified = 0, // chirality that hasn't been specified
    ChiTetrahedralCw,   // tetrahedral: clockwise rotation (smiles \@\@)
    ChiTetrahedralCcw,  // tetrahedral: counter-clockwise rotation (smiles

    ChiOther,               // some unrecognized type of chirality
    ChiTetrahedral,         // tetrahedral, use permutation flag
    ChiAllene,              // allene, use permutation flag
    ChiSquareplanar,        // square planar, use permutation flag
    ChiTrigonalbipyramidal, // trigonal bipyramidal, use permutation flag
    ChiOctahedral,          // octahedral, use permutation flag
}

#[derive(Debug)]
pub struct AtomWithMol {
    pub atom: Atom,
    // Note: while informally expected, this molecule is not guaranteed to include the atom
    pub atom_mol: Option<AtomMol>,
}

impl AtomWithMol {
    fn calc_implicit_valence(&self, strict: bool) -> Result<i32> {
        calc_implicit_valence(self, strict, false)
    }

    fn calc_explicit_valence(&self, strict: bool) -> Result<i32> {
        calc_explicit_valence(self, strict, false)
    }

    fn degree(&self) -> Option<i32> {
        self.atom_mol
            .as_ref()
            .map(|a| a.mol.atom_degree(&self.atom))
    }

    fn total_degree(&self) -> Result<Option<i32>> {
        if let Some(degree) = self.degree() {
            Ok(Some(self.total_num_hs(false)? + degree))
        } else {
            Ok(None)
        }
    }

    fn total_num_hs(&self, include_neighbors: bool) -> Result<i32> {
        let mut res = self.atom.num_explicit_hs + self.num_implicit_hs()?;

        if let Some(mol) = &self.atom_mol {
            if include_neighbors {
                for nbr in mol.mol.atom_neighbors(&self.atom) {
                    if nbr.atomic_num.0 == 1 {
                        res += 1;
                    }
                }
            }
        }

        Ok(res)
    }

    fn num_implicit_hs(&self) -> Result<i32> {
        if self.atom.no_implicit {
            return Ok(0);
        }
        self.implicit_valence()
    }

    fn implicit_valence(&self) -> Result<i32> {
        if let Some(atom_mol) = &self.atom_mol {
            Ok(atom_mol.implicit_valence)
        } else {
            Err(anyhow!(
                "valence not defined for atoms not associated with molecules"
            ))
        }
    }

    // port: it's weird that explicit valence checks the field and implicit_valence the mol optional
    // probably it's equivalent in practice, but they should do the same check, for consistency
    // (or this logic should be changed more thoroughly)
    fn explicit_valence(&self) -> Result<i32> {
        if self.atom.explicit_valence < 0 {
            return Err(anyhow!(
                "getExplicitValence() called without call to calcExplicitValence()"
            ));
        }
        Ok(self.atom.explicit_valence)
    }

    fn total_valence(&self) -> Result<i32> {
        Ok(self.implicit_valence()? + self.explicit_valence()?)
    }

    fn has_valence_violation(&self) -> Result<bool> {
        // Ignore dummy atoms, query atoms, or atoms attached to query bonds
        if let Some(mol) = &self.atom_mol {
            let mut bonds = mol.mol.atom_bonds(&self.atom);
            // port: TODO query things
            // auto is_query = [](auto b) { return b->hasQuery(); };
            if self.atom.atomic_num.0 == 0
                || self.atom.has_query()
                || bonds.any(|b| b.bond.has_query())
            {
                return Ok(false);
            }
        }

        let check_it = true;
        let effective_atomic_num_res = get_effective_atom_num(&self.atom, check_it);
        if effective_atomic_num_res.is_err() {
            return Ok(true);
        }
        let effective_atomic_num = effective_atomic_num_res?;

        // special case for H:
        if self.atom.atomic_num.0 == 1 {
            if self.atom.formal_charge > 1 || self.atom.formal_charge < -1 {
                return Ok(true);
            }
        } else {
            // Non-H checks for absurd charge values:
            //   1. the formal charge is larger than the atomic number
            //   2. the formal charge moves us to a different row of the periodic table
            let pd = PeriodicTable::new();

            if self.atom.formal_charge as usize > self.atom.atomic_num.0
                || pd.get_row(self.atom.atomic_num) != pd.get_row(effective_atomic_num)
            {
                return Ok(true);
            }
        }

        let strict = false;
        let check_it = true;
        if calc_explicit_valence(self, strict, check_it)? == -1
            || calc_implicit_valence(self, strict, check_it)? == -1
        {
            return Ok(true);
        }

        Ok(false)
    }
}

// port: TODO review number types, at least use u* for positive
#[derive(Debug, PartialEq, Eq)]
pub struct Atom {
    pub is_aromatic: bool,
    pub no_implicit: bool,
    pub num_explicit_hs: i32,
    pub formal_charge: i32,
    pub atomic_num: AtomicNumber,
    pub implicit_valence: i32,
    pub explicit_valence: i32,
    pub num_radical_electrons: i32,
    pub chiral_tag: ChiralType,
    pub hybrid: u8,
    pub isotope: u32,
    pub index: usize,
}

/// Molecule associated with an atom and related data
#[derive(Debug)]
struct AtomMol {
    mol: Mol,
    implicit_valence: i32,
}

impl AtomMol {}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            is_aromatic: false,
            no_implicit: false,
            num_explicit_hs: 0,
            formal_charge: 0,
            atomic_num: AtomicNumber(0),
            implicit_valence: 0,
            explicit_valence: 0,
            num_radical_electrons: 0,
            chiral_tag: ChiralType::ChiUnspecified,
            hybrid: 0,
            isotope: 0,
            index: 0,
        }
    }
}

impl Atom {
    fn with_atomic_number(num: AtomicNumber) -> Atom {
        Atom {
            atomic_num: num,
            ..Default::default()
        }
    }

    fn get_mass(&self) -> f64 {
        // if self.iso
        let pd = PeriodicTable::new();
        if self.isotope != 0 {
            let mut res = pd.mass_for_isotope(self.atomic_num, self.isotope);
            if self.atomic_num.0 != 0 && res == 0.0 {
                res = self.isotope as f64;
            }
            res
        } else {
            pd.atomic_weight(self.atomic_num)
        }
    }

    // This method can be used to distinguish query atoms from standard atoms:
    pub fn has_query(&self) -> bool {
        false
    }

    pub fn invert_chirality(&mut self) -> Result<bool> {
        match self.chiral_tag {
            ChiralType::ChiTetrahedralCw => {
                self.chiral_tag = ChiralType::ChiTetrahedralCcw;
                return Ok(true);
            }
            ChiralType::ChiTetrahedralCcw => {
                self.chiral_tag = ChiralType::ChiTetrahedralCw;
                return Ok(true);
            }
            ChiralType::ChiTetrahedral => {
                let mut common_props = common_props()?;
                let perm = common_props.get_mut(&CommonPropsKey::ChiralPermutation);
                if let Some(CommonPropsValue::U32(perm)) = perm {
                    match perm {
                        1 => *perm = 2,
                        2 => *perm = 1,
                        _ => *perm = 0,
                    }
                    return Ok(*perm != 0);
                } else {
                    return Ok(false);
                }
            }
            ChiralType::ChiTrigonalbipyramidal => {
                let mut common_props = common_props()?;
                let perm = common_props.get_mut(&CommonPropsKey::ChiralPermutation);
                if let Some(CommonPropsValue::U32(perm)) = perm {
                    let perm_res = if *perm <= 20 {
                        trigonalbipyramidal_invert()[*perm as usize]
                    } else {
                        0
                    };
                    *perm = perm_res as u32;
                    return Ok(perm_res != 0);
                } else {
                    return Ok(false);
                }
            }
            ChiralType::ChiOctahedral => {
                let mut common_props = common_props()?;
                let perm = common_props.get_mut(&CommonPropsKey::ChiralPermutation);
                if let Some(CommonPropsValue::U32(perm)) = perm {
                    let perm_res = if *perm <= 30 {
                        trigonalbipyramidal_invert()[*perm as usize]
                    } else {
                        0
                    };
                    *perm = perm_res as u32;
                    return Ok(perm_res != 0);
                } else {
                    return Ok(false);
                }
            }
            _ => {
                return Ok(false);
            }
        }
    }
}

// NOTE: this uses the explicitValence, so it will call
// calc_explicit_valence if it is not set on the given atom
fn calc_implicit_valence(atom: &AtomWithMol, strict: bool, check_it: bool) -> Result<i32> {
    if atom.atom.no_implicit {
        return Ok(0);
    }

    let mut explicit_valence = atom.atom.explicit_valence;
    if explicit_valence == -1 {
        explicit_valence = calc_explicit_valence(atom, strict, check_it)?;
    }
    let atomic_num = atom.atom.atomic_num;
    if atomic_num.0 == 0 {
        return Ok(0);
    }

    let formal_charge = atom.atom.formal_charge;
    let num_radical_electrons = atom.atom.num_radical_electrons;

    if explicit_valence == 0 && num_radical_electrons == 0 && atomic_num.0 == 1 {
        if formal_charge == 1 || formal_charge == -1 {
            return Ok(0);
        } else if formal_charge == 0 {
            return Ok(1);
        } else {
            if strict {
                return Err(anyhow!(
                    "Unreasonable formal charge on atom # {}",
                    atom.atom.index
                ));
            } else if check_it {
                return Ok(-1);
            } else {
                return Ok(0);
            }
        }
    }

    let mut explicit_plus_rad_v = atom.atom.explicit_valence + atom.atom.num_radical_electrons;

    let periodic_table = PeriodicTable::new();
    let ovalens = periodic_table.valence_list(atom.atom.atomic_num);

    // if we start with an atom that doesn't have specified valences, we stick
    // with that. otherwise we will use the effective valence for the rest of
    // this.

    let mut effective_atomic_num = atom.atom.atomic_num;
    if ovalens.len() > 1 || ovalens.first() != Some(&1) {
        effective_atomic_num = get_effective_atom_num(&atom.atom, check_it)?;
    }
    if effective_atomic_num.0 == 0 {
        return Ok(0);
    }

    // this is basically the difference between the allowed valence of
    // the atom and the explicit valence already specified - tells how
    // many Hs to add
    //

    // The d-block and f-block of the periodic table (i.e. transition metals,
    // lanthanoids and actinoids) have no default valence.
    let dv = periodic_table.default_valence(atom.atom.atomic_num);
    if dv == -1 {
        return Ok(0);
    }

    // here is how we are going to deal with the possibility of
    // multiple valences
    // - check the explicit valence "ev"
    // - if it is already equal to one of the allowed valences for the
    //    atom return 0
    // - otherwise take return difference between next larger allowed
    //   valence and "ev"
    // if "ev" is greater than all allowed valences for the atom raise an
    // exception
    // finally aromatic cases are dealt with differently - these atoms are allowed
    // only default valences

    // we have to include a special case here for negatively charged P, S, As,
    // and Se, which all support "hypervalent" forms, but which can be
    // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.

    if can_be_hypervalent(&atom.atom, effective_atomic_num) {
        effective_atomic_num = atomic_num;
        explicit_plus_rad_v -= atom.atom.formal_charge;
    }

    let valens = periodic_table.valence_list(effective_atomic_num);

    let mut res;

    if is_aromatic_atom(&atom) {
        if explicit_plus_rad_v <= dv {
            res = dv - explicit_plus_rad_v;
        } else {
            // As we assume when finding the explicitPlusRadValence if we are
            // aromatic we should not be adding any hydrogen and already
            // be at an accepted valence state,

            // FIX: this is just ERROR checking and probably moot - the
            // explicitPlusRadValence function called above should assure us that
            // we satisfy one of the accepted valence states for the
            // atom. The only diff I can think of is in the way we handle
            // formal charge here vs the explicit valence function.

            let mut satis = false;
            for vi in valens {
                if explicit_plus_rad_v == vi {
                    satis = true;
                    break;
                }
            }

            if !satis && (strict || check_it) {
                if strict {
                    return Err(anyhow!("Explicit valence for aromatic atom # {} not equal to any accepted valence\n", atom.atom.index));
                } else {
                    return Ok(-1);
                }
            }
            res = 0;
        }
    } else {
        // non-aromatic case we are allowed to have non default valences
        // and be able to add hydrogens

        res = -1;
        for vi in &valens {
            let tot = vi;
            if explicit_plus_rad_v <= *tot {
                res = tot - explicit_plus_rad_v;
                break;
            }
        }

        if res < 0 {
            if (strict || check_it) && valens.last() != Some(&-1) && ovalens.last() > Some(&0) {
                // this means that the explicit valence is greater than any
                // allowed valence for the atoms
                if strict {
                    // raise an error
                    return Err(anyhow!(
                        "Explicit valence for atom # {} {} greater than permitted",
                        atom.atom.index,
                        periodic_table.element_symbol(atom.atom.atomic_num)
                    ));
                }
            } else {
                res = 0
            }
        }
    }

    Ok(res)
}

fn calc_explicit_valence(atom: &AtomWithMol, strict: bool, check_it: bool) -> Result<i32> {
    // FIX: contributions of bonds to valence are being done at best
    // approximately

    let mut accum: f64 = 0.;
    if let Some(mol) = &atom.atom_mol {
        for bond in mol.mol.atom_bonds(&atom.atom) {
            accum += bond.get_valence_contrib(&atom.atom);
        }
    }
    accum += atom.atom.num_explicit_hs as f64;

    let periodic_table = PeriodicTable::new();
    let ovalens = periodic_table.valence_list(atom.atom.atomic_num);

    // if we start with an atom that doesn't have specified valences, we stick
    // with that. otherwise we will use the effective valence
    let mut effective_atomic_num = atom.atom.atomic_num;
    if ovalens.len() > 1
    // port: TODO: why/when can this be -1?
    // || ovalens[0] != -1
    {
        effective_atomic_num = get_effective_atom_num(&atom.atom, check_it)?;
    }

    let dv = periodic_table.default_valence(atom.atom.atomic_num);
    let valens = periodic_table.valence_list(atom.atom.atomic_num);

    if accum > dv as f64 && is_aromatic_atom(&atom) {
        // this needs some explanation : if the atom is aromatic and
        // accum > dv we assume that no hydrogen can be added
        // to this atom.  We set x = (v + chr) such that x is the
        // closest possible integer to "accum" but less than
        // "accum".
        //
        // "v" here is one of the allowed valences. For example:
        //    sulfur here : O=c1ccs(=O)cc1
        //    nitrogen here : c1cccn1C

        let mut pval = dv;
        for val in &valens {
            // port: TODO review negative valence
            if *val == -1 {
                break;
            }
            if *val as f64 > accum {
                break;
            } else {
                pval = *val;
            }
        }
        // if we're within 1.5 of the allowed valence, go ahead and take it.
        // this reflects things like the N in c1cccn1C, which starts with
        // accum of 4, but which can be kekulized to C1=CC=CN1C, where
        // the valence is 3 or the bridging N in c1ccn2cncc2c1, which starts
        // with a valence of 4.5, but can be happily kekulized down to a valence
        // of 3
        if accum - pval as f64 <= 1.5 {
            accum = pval as f64;
        }
    }

    // despite promising to not to blame it on him - this a trick Greg
    // came up with: if we have a bond order sum of x.5 (i.e. 1.5, 2.5
    // etc) we would like it to round to the higher integer value --
    // 2.5 to 3 instead of 2 -- so we will add 0.1 to accum.
    // this plays a role in the number of hydrogen that are implicitly
    // added. This will only happen when the accum is a non-integer
    // value and less than the default valence (otherwise the above if
    // statement should have caught it). An example of where this can
    // happen is the following smiles:
    //    C1ccccC1
    // Daylight accepts this smiles and we should be able to Kekulize
    // correctly.
    accum += 0.1;

    let res = accum.round() as i32;

    if strict || check_it {
        // port: TODO review unwrap
        let mut max_valence = valens.last().unwrap();
        let mut offset = 0;

        // we have to include a special case here for negatively charged P, S, As,
        // and Se, which all support "hypervalent" forms, but which can be
        // isoelectronic to Cl/Ar or Br/Kr, which do not support hypervalent forms.
        if can_be_hypervalent(&atom.atom, effective_atomic_num) {
            // port: TODO review unwrap
            max_valence = valens.last().unwrap();
            offset -= atom.atom.formal_charge;
        }
        // we have historically accepted two-coordinate [H-] as a valid atom. This
        // is highly questionable, but changing it requires some thought. For now we
        // will just keep accepting it
        if atom.atom.atomic_num.0 == 1 && atom.atom.formal_charge == -1 {
            max_valence = &2;
        }

        // maxValence == -1 signifies that we'll take anything at the high end
        // port: TODO review unwrap
        if *max_valence >= 0
            && *ovalens.last().unwrap() >= 0
            && (res + offset) > *max_valence as i32
        {
            if strict {
                return Err(anyhow!(
                    "Explicit valence for atom # {} {} is greater than permitted",
                    atom.atom.index,
                    periodic_table.element_symbol(atom.atom.atomic_num)
                ));
            }
        }
    }

    Ok(res)
}

fn can_be_hypervalent(atom: &Atom, effective_atomic_num: AtomicNumber) -> bool {
    (effective_atomic_num.0 > 16 && (atom.atomic_num.0 == 15 || atom.atomic_num.0 == 16))
        || (effective_atomic_num.0 > 34 && (atom.atomic_num.0 == 33 || atom.atomic_num.0 == 34))
}

fn get_effective_atom_num(atom: &Atom, check_value: bool) -> Result<AtomicNumber> {
    let effective_atom_num = atom.atomic_num.0 - atom.formal_charge as usize;

    if check_value {
        if effective_atom_num < 0 || effective_atom_num > PeriodicTable::max_atom_num().0 {
            return Err(anyhow!(
                "Effective atomic number out of range: {:?}",
                atom.atomic_num
            ));
        }
    }

    Ok(AtomicNumber(
        atom.atomic_num.0.clamp(0, PeriodicTable::max_atom_num().0),
    ))
}

// port: I don't like this function, if one calls "getIsAromatic" (c++) it can give a different result than "isAromaticAtom"
// the atom being aromatic or not should always mean the same
// also TODO fields isAromatic and type Bond::BondType::AROMATIC is redundant, former should be derived
pub fn is_aromatic_atom(atom: &AtomWithMol) -> bool {
    if atom.atom.is_aromatic {
        return true;
    }

    if let Some(mol) = &atom.atom_mol {
        for bond in mol.mol.atom_bonds(&atom.atom) {
            // port: also not good, should check/use only one of these fields (bond_type probably)
            if bond.bond.is_aromatic || bond.bond.bond_type == BondType::Aromatic {
                return true;
            }
        }
    }
    return false;
}

// Determine whether or not an element is to the left of carbon.
fn is_early_atom(atomic_num: AtomicNumber) -> bool {
    let table = [
        false, // #0 *
        false, // #1 H
        false, // #2 He
        true,  // #3 Li
        true,  // #4 Be
        true,  // #5 B
        false, // #6 C
        false, // #7 N
        false, // #8 O
        false, // #9 F
        false, // #10 Ne
        true,  // #11 Na
        true,  // #12 Mg
        true,  // #13 Al
        false, // #14 Si
        false, // #15 P
        false, // #16 S
        false, // #17 Cl
        false, // #18 Ar
        true,  // #19 K
        true,  // #20 Ca
        true,  // #21 Sc
        true,  // #22 Ti
        false, // #23 V
        false, // #24 Cr
        false, // #25 Mn
        false, // #26 Fe
        false, // #27 Co
        false, // #28 Ni
        false, // #29 Cu
        true,  // #30 Zn
        true,  // #31 Ga
        true,  // #32 Ge  see github #2606
        false, // #33 As
        false, // #34 Se
        false, // #35 Br
        false, // #36 Kr
        true,  // #37 Rb
        true,  // #38 Sr
        true,  // #39 Y
        true,  // #40 Zr
        true,  // #41 Nb
        false, // #42 Mo
        false, // #43 Tc
        false, // #44 Ru
        false, // #45 Rh
        false, // #46 Pd
        false, // #47 Ag
        true,  // #48 Cd
        true,  // #49 In
        true,  // #50 Sn  see github #2606
        true,  // #51 Sb  see github #2775
        false, // #52 Te
        false, // #53 I
        false, // #54 Xe
        true,  // #55 Cs
        true,  // #56 Ba
        true,  // #57 La
        true,  // #58 Ce
        true,  // #59 Pr
        true,  // #60 Nd
        true,  // #61 Pm
        false, // #62 Sm
        false, // #63 Eu
        false, // #64 Gd
        false, // #65 Tb
        false, // #66 Dy
        false, // #67 Ho
        false, // #68 Er
        false, // #69 Tm
        false, // #70 Yb
        false, // #71 Lu
        true,  // #72 Hf
        true,  // #73 Ta
        false, // #74 W
        false, // #75 Re
        false, // #76 Os
        false, // #77 Ir
        false, // #78 Pt
        false, // #79 Au
        true,  // #80 Hg
        true,  // #81 Tl
        true,  // #82 Pb  see github #2606
        true,  // #83 Bi  see github #2775
        false, // #84 Po
        false, // #85 At
        false, // #86 Rn
        true,  // #87 Fr
        true,  // #88 Ra
        true,  // #89 Ac
        true,  // #90 Th
        true,  // #91 Pa
        true,  // #92 U
        true,  // #93 Np
        false, // #94 Pu
        false, // #95 Am
        false, // #96 Cm
        false, // #97 Bk
        false, // #98 Cf
        false, // #99 Es
        false, // #100 Fm
        false, // #101 Md
        false, // #102 No
        false, // #103 Lr
        true,  // #104 Rf
        true,  // #105 Db
        true,  // #106 Sg
        true,  // #107 Bh
        true,  // #108 Hs
        true,  // #109 Mt
        true,  // #110 Ds
        true,  // #111 Rg
        true,  // #112 Cn
        true,  // #113 Nh
        true,  // #114 Fl
        true,  // #115 Mc
        true,  // #116 Lv
        true,  // #117 Ts
        true,  // #118 Og
    ];
    if atomic_num.0 < 119 {
        table[atomic_num.0]
    } else {
        false
    }
}

fn octahedral_invert() -> [u8; 31] {
    [
        0,  //  0 -> 0
        2,  //  1 -> 2
        1,  //  2 -> 1
        16, //  3 -> 16
        14, //  4 -> 14
        15, //  5 -> 15
        18, //  6 -> 18
        17, //  7 -> 17
        10, //  8 -> 10
        11, //  9 -> 11
        8,  // 10 -> 8
        9,  // 11 -> 9
        13, // 12 -> 13
        12, // 13 -> 12
        4,  // 14 -> 4
        5,  // 15 -> 5
        3,  // 16 -> 3
        7,  // 17 -> 7
        6,  // 18 -> 6
        24, // 19 -> 24
        23, // 20 -> 23
        22, // 21 -> 22
        21, // 22 -> 21
        20, // 23 -> 20
        19, // 24 -> 19
        30, // 25 -> 30
        29, // 26 -> 29
        28, // 27 -> 28
        27, // 28 -> 27
        26, // 29 -> 26
        25, // 30 -> 25
    ]
}

fn trigonalbipyramidal_invert() -> [u8; 21] {
    [
        0,  //  0 -> 0
        2,  //  1 -> 2
        1,  //  2 -> 1
        4,  //  3 -> 4
        3,  //  4 -> 3
        6,  //  5 -> 6
        5,  //  6 -> 5
        8,  //  7 -> 8
        7,  //  8 -> 7
        11, //  9 -> 11
        12, // 10 -> 12
        9,  // 11 -> 9
        10, // 12 -> 10
        14, // 13 -> 14
        13, // 14 -> 13
        20, // 15 -> 20
        19, // 16 -> 19
        18, // 17 -> 28
        17, // 18 -> 17
        16, // 19 -> 16
        15, // 20 -> 15
    ]
}
