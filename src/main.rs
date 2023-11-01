use std::collections::HashMap;
use std::str::FromStr;

use molrs::{
    atom::{Atom, PointChirality},
    bond::{Bond, BondType},
    molecule::{Molecule, MoleculeError},
};
use pertable::Element;

enum AtomAttribute {
    Isotope,
    NImplicitHydrogens,
    Charge,
}

fn smiles_parser_v1(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];

    for c in smi.chars() {
        if c == 'B'
            || c == 'C'
            || c == 'N'
            || c == 'O'
            || c == 'P'
            || c == 'S'
            || c == 'F'
            || c == 'I'
            || c == '*'
        {
            atoms.push(Atom {
                element: Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap(),
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: molrs::atom::PointChirality::Undefined,
            });
        } else if c == 'l' && atoms.last().unwrap().element == Element::C {
            atoms.last_mut().unwrap().element = Element::Cl;
        } else if c == 'r' && atoms.last().unwrap().element == Element::B {
            atoms.last_mut().unwrap().element = Element::Br;
        } else {
            return Err(MoleculeError::SmilesParseError(format!(
                "{smi} | invalid char {c}"
            )));
        };
    }

    Ok(Molecule {
        atoms,
        bonds: vec![],
        rings: None,
    })
}

fn smiles_parser_v2(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];

    // *** NEW CODE BEGINS ***
    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();
    // *** NEW CODE ENDS ***

    for c in smi.chars() {
        // *** NEW CODE BEGINS ***
        if c_is_in_bracket {
            let atom: &mut Atom = atoms.last_mut().unwrap();
            if c == ']' {
                c_is_in_bracket = false;
                atom.element = match element_str.parse() {
                    Ok(element) => element,
                    Err(_) => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | invalid element {element_str}"
                        )));
                    }
                };
                if element_str.chars().next().unwrap().is_lowercase() {
                    atom.delocalized = true;
                };
                element_str = String::new();
            } else if c == 'H' {
                if element_str.is_empty() {
                    element_str.push('H');
                } else {
                    atom.n_implicit_hydrogens = Some(1);
                    atom_attribute = AtomAttribute::NImplicitHydrogens;
                };
            } else if c == '+' {
                atom.charge = 1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '-' {
                atom.charge = -1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '@' {
                match atom.point_chirality {
                    PointChirality::Undefined => {
                        atom.point_chirality = PointChirality::CounterClockwise;
                    }
                    PointChirality::CounterClockwise => {
                        atom.point_chirality = PointChirality::Clockwise;
                    }
                    _ => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | chirality error"
                        )));
                    }
                }
            } else if c.is_alphabetic() {
                element_str.push(c);
            } else if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap();
                match atom_attribute {
                    AtomAttribute::Isotope => match atom.isotope {
                        None => {
                            atom.isotope = Some(c_as_digit as u16);
                        }
                        Some(isotope) => atom.isotope = Some(isotope * 10 + c_as_digit as u16),
                    },
                    AtomAttribute::Charge => {
                        atom.charge *= c_as_digit as i8;
                    }
                    AtomAttribute::NImplicitHydrogens => {
                        atom.n_implicit_hydrogens = Some(c_as_digit as u8);
                    }
                };
            };
        } else if c == '[' {
            c_is_in_bracket = true;
            atom_attribute = AtomAttribute::Isotope;
            atoms.push(Atom::default());
        // *** NEW CODE ENDS ***
        } else if c == 'B'
            || c == 'C'
            || c == 'N'
            || c == 'O'
            || c == 'P'
            || c == 'S'
            || c == 'F'
            || c == 'I'
            || c == '*'
        {
            atoms.push(Atom {
                element: Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap(),
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: molrs::atom::PointChirality::Undefined,
            });
        } else if c == 'l' && atoms.last().unwrap().element == Element::C {
            atoms.last_mut().unwrap().element = Element::Cl;
        } else if c == 'r' && atoms.last().unwrap().element == Element::B {
            atoms.last_mut().unwrap().element = Element::Br;
        } else {
            return Err(MoleculeError::SmilesParseError(format!(
                "{smi} | invalid char {c}"
            )));
        }
    }

    Ok(Molecule {
        atoms,
        bonds: vec![],
        rings: None,
    })
}

fn smiles_parser_v3(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];
    // *** NEW CODE BEGINS ***
    let mut bond: Bond = Bond::default();
    let mut bonds: Vec<Bond> = vec![];
    // *** NEW CODE ENDS ***

    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();

    // *** NEW CODE BEGINS ***
    for (i, c) in smi.chars().enumerate() {
        // *** NEW CODE ENDS ***
        if c_is_in_bracket {
            let atom: &mut Atom = atoms.last_mut().unwrap();
            if c == ']' {
                c_is_in_bracket = false;
                atom.element = match element_str.parse() {
                    Ok(element) => element,
                    Err(_) => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | invalid element {element_str}"
                        )));
                    }
                };
                if element_str.chars().next().unwrap().is_lowercase() {
                    atom.delocalized = true;
                };
                element_str = String::new();
            } else if c == 'H' {
                if element_str.is_empty() {
                    element_str.push('H');
                } else {
                    atom.n_implicit_hydrogens = Some(1);
                    atom_attribute = AtomAttribute::NImplicitHydrogens;
                };
            } else if c == '+' {
                atom.charge = 1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '-' {
                atom.charge = -1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '@' {
                match atom.point_chirality {
                    PointChirality::Undefined => {
                        atom.point_chirality = PointChirality::CounterClockwise;
                    }
                    PointChirality::CounterClockwise => {
                        atom.point_chirality = PointChirality::Clockwise;
                    }
                    _ => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | chirality error"
                        )));
                    }
                }
            } else if c.is_alphabetic() {
                element_str.push(c);
            } else if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap();
                match atom_attribute {
                    AtomAttribute::Isotope => match atom.isotope {
                        None => {
                            atom.isotope = Some(c_as_digit as u16);
                        }
                        Some(isotope) => atom.isotope = Some(isotope * 10 + c_as_digit as u16),
                    },
                    AtomAttribute::Charge => {
                        atom.charge *= c_as_digit as i8;
                    }
                    AtomAttribute::NImplicitHydrogens => {
                        atom.n_implicit_hydrogens = Some(c_as_digit as u8);
                    }
                };
            };
        } else if c == '[' {
            c_is_in_bracket = true;
            atom_attribute = AtomAttribute::Isotope;
            atoms.push(Atom::default());
            // *** NEW CODE BEGINS ***
            if i > 0 {
                bond.i = i - 1;
                bond.j = i;
                bonds.push(bond);
                bond = Bond::default();
            }
            // *** NEW CODE ENDS ***
            // *** NEW CODE BEGINS ***
        } else if c == '-' || c == '/' || c == '\\' || c == ':' || c == '=' || c == '#' || c == '$'
        {
            bond.bond_type = BondType::try_from(c).unwrap();
        // *** NEW CODE ENDS ***
        } else if c == 'B'
            || c == 'C'
            || c == 'N'
            || c == 'O'
            || c == 'P'
            || c == 'S'
            || c == 'F'
            || c == 'I'
            || c == '*'
        {
            atoms.push(Atom {
                element: Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap(),
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: molrs::atom::PointChirality::Undefined,
            });
            // *** NEW CODE BEGINS ***
            if i > 0 {
                bond.i = i - 1;
                bond.j = i;
                bonds.push(bond);
                bond = Bond::default();
            }
            // *** NEW CODE ENDS ***
        } else if c == 'l' && atoms.last().unwrap().element == Element::C {
            atoms.last_mut().unwrap().element = Element::Cl;
        } else if c == 'r' && atoms.last().unwrap().element == Element::B {
            atoms.last_mut().unwrap().element = Element::Br;
        } else {
            return Err(MoleculeError::SmilesParseError(format!(
                "{smi} | invalid char {c}"
            )));
        }
    }

    Ok(Molecule {
        atoms,
        bonds,
        rings: None,
    })
}

fn smiles_parser_v4(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];
    let mut bond: Bond = Bond::default();
    let mut bonds: Vec<Bond> = vec![];

    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();
    // *** NEW CODE BEGINS ***
    let mut root_atom: Vec<usize> = vec![];
    // *** NEW CODE ENDS ***

    // *** NEW CODE BEGINS ***
    for c in smi.chars() {
        // *** NEW CODE ENDS ***
        if c_is_in_bracket {
            let atom: &mut Atom = atoms.last_mut().unwrap();
            if c == ']' {
                c_is_in_bracket = false;
                atom.element = match element_str.parse() {
                    Ok(element) => element,
                    Err(_) => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | invalid element {element_str}"
                        )));
                    }
                };
                if element_str.chars().next().unwrap().is_lowercase() {
                    atom.delocalized = true;
                };
                element_str = String::new();
            } else if c == 'H' {
                if element_str.is_empty() {
                    element_str.push('H');
                } else {
                    atom.n_implicit_hydrogens = Some(1);
                    atom_attribute = AtomAttribute::NImplicitHydrogens;
                };
            } else if c == '+' {
                atom.charge = 1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '-' {
                atom.charge = -1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '@' {
                match atom.point_chirality {
                    PointChirality::Undefined => {
                        atom.point_chirality = PointChirality::CounterClockwise;
                    }
                    PointChirality::CounterClockwise => {
                        atom.point_chirality = PointChirality::Clockwise;
                    }
                    _ => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | chirality error"
                        )));
                    }
                }
            } else if c.is_alphabetic() {
                element_str.push(c);
            } else if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap();
                match atom_attribute {
                    AtomAttribute::Isotope => match atom.isotope {
                        None => {
                            atom.isotope = Some(c_as_digit as u16);
                        }
                        Some(isotope) => atom.isotope = Some(isotope * 10 + c_as_digit as u16),
                    },
                    AtomAttribute::Charge => {
                        atom.charge *= c_as_digit as i8;
                    }
                    AtomAttribute::NImplicitHydrogens => {
                        atom.n_implicit_hydrogens = Some(c_as_digit as u8);
                    }
                };
            };
        } else if c == '[' {
            c_is_in_bracket = true;
            atom_attribute = AtomAttribute::Isotope;
            atoms.push(Atom::default());
            // *** NEW CODE BEGINS ***
            if !root_atom.is_empty() {
                bond.i = root_atom.pop().unwrap();
                bond.j = atoms.len() - 1;
                bonds.push(bond);
                bond = Bond::default();
            }
            root_atom.push(atoms.len() - 1);
            // *** NEW CODE ENDS ***
        } else if c == '-' || c == '/' || c == '\\' || c == ':' || c == '=' || c == '#' || c == '$'
        {
            bond.bond_type = BondType::try_from(c).unwrap();
        // *** NEW CODE BEGINS ***
        } else if c == '(' {
            root_atom.push(*root_atom.last().unwrap());
        } else if c == ')' {
            root_atom.pop();
        // *** NEW CODE ENDS ***
        } else if c == 'B'
            || c == 'C'
            || c == 'N'
            || c == 'O'
            || c == 'P'
            || c == 'S'
            || c == 'F'
            || c == 'I'
            || c == '*'
        {
            atoms.push(Atom {
                element: Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap(),
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: molrs::atom::PointChirality::Undefined,
            });
            // *** NEW CODE BEGINS ***
            if !root_atom.is_empty() {
                bond.i = root_atom.pop().unwrap();
                bond.j = atoms.len() - 1;
                bonds.push(bond);
                bond = Bond::default();
            }
            root_atom.push(atoms.len() - 1);
            // *** NEW CODE ENDS ***
        } else if c == 'l' && atoms.last().unwrap().element == Element::C {
            atoms.last_mut().unwrap().element = Element::Cl;
        } else if c == 'r' && atoms.last().unwrap().element == Element::B {
            atoms.last_mut().unwrap().element = Element::Br;
        } else {
            return Err(MoleculeError::SmilesParseError(format!(
                "{smi} | invalid char {c}"
            )));
        }
    }

    Ok(Molecule {
        atoms,
        bonds,
        rings: None,
    })
}

fn smiles_parser_v5(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];
    let mut bond: Bond = Bond::default();
    let mut bonds: Vec<Bond> = vec![];
    // *** NEW CODE BEGINS ***
    let mut ring_closures: HashMap<usize, usize> = HashMap::new();
    // *** NEW CODE ENDS ***

    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();
    let mut root_atom: Vec<usize> = vec![];

    for c in smi.chars() {
        if c_is_in_bracket {
            let atom: &mut Atom = atoms.last_mut().unwrap();
            if c == ']' {
                c_is_in_bracket = false;
                atom.element = match element_str.parse() {
                    Ok(element) => element,
                    Err(_) => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | invalid element {element_str}"
                        )));
                    }
                };
                if element_str.chars().next().unwrap().is_lowercase() {
                    atom.delocalized = true;
                };
                element_str = String::new();
            } else if c == 'H' {
                if element_str.is_empty() {
                    element_str.push('H');
                } else {
                    atom.n_implicit_hydrogens = Some(1);
                    atom_attribute = AtomAttribute::NImplicitHydrogens;
                };
            } else if c == '+' {
                atom.charge = 1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '-' {
                atom.charge = -1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '@' {
                match atom.point_chirality {
                    PointChirality::Undefined => {
                        atom.point_chirality = PointChirality::CounterClockwise;
                    }
                    PointChirality::CounterClockwise => {
                        atom.point_chirality = PointChirality::Clockwise;
                    }
                    _ => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | chirality error"
                        )));
                    }
                }
            } else if c.is_alphabetic() {
                element_str.push(c);
            } else if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap();
                match atom_attribute {
                    AtomAttribute::Isotope => match atom.isotope {
                        None => {
                            atom.isotope = Some(c_as_digit as u16);
                        }
                        Some(isotope) => atom.isotope = Some(isotope * 10 + c_as_digit as u16),
                    },
                    AtomAttribute::Charge => {
                        atom.charge *= c_as_digit as i8;
                    }
                    AtomAttribute::NImplicitHydrogens => {
                        atom.n_implicit_hydrogens = Some(c_as_digit as u8);
                    }
                };
            };
        } else if c == '[' {
            c_is_in_bracket = true;
            atom_attribute = AtomAttribute::Isotope;
            atoms.push(Atom::default());
            if !root_atom.is_empty() {
                bond.i = root_atom.pop().unwrap();
                bond.j = atoms.len() - 1;
                bonds.push(bond);
                bond = Bond::default();
            }
            root_atom.push(atoms.len() - 1);
        } else if c == '-' || c == '/' || c == '\\' || c == ':' || c == '=' || c == '#' || c == '$'
        {
            bond.bond_type = BondType::try_from(c).unwrap();
        // *** NEW CODE BEGINS ***
        } else if c.is_numeric() {
            let c_as_digit: usize = c.to_digit(10).unwrap() as usize;
            if let std::collections::hash_map::Entry::Vacant(e) = ring_closures.entry(c_as_digit) {
                e.insert(atoms.len() - 1);
            } else {
                let ring_closure_bond = Bond {
                    i: *ring_closures.get(&c_as_digit).unwrap(),
                    j: atoms.len() - 1,
                    bond_type: bond.bond_type,
                };
                bond.bond_type = BondType::Default;
                bonds.push(ring_closure_bond);
                ring_closures.remove(&c_as_digit);
            }
        // *** NEW CODE ENDS ***
        } else if c == '(' {
            root_atom.push(*root_atom.last().unwrap());
        } else if c == ')' {
            root_atom.pop();
        } else if c == 'B'
            || c == 'C'
            || c == 'N'
            || c == 'O'
            || c == 'P'
            || c == 'S'
            || c == 'F'
            || c == 'I'
            || c == '*'
        {
            atoms.push(Atom {
                element: Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap(),
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: molrs::atom::PointChirality::Undefined,
            });
            if !root_atom.is_empty() {
                bond.i = root_atom.pop().unwrap();
                bond.j = atoms.len() - 1;
                bonds.push(bond);
                bond = Bond::default();
            }
            root_atom.push(atoms.len() - 1);
        } else if c == 'l' && atoms.last().unwrap().element == Element::C {
            atoms.last_mut().unwrap().element = Element::Cl;
        } else if c == 'r' && atoms.last().unwrap().element == Element::B {
            atoms.last_mut().unwrap().element = Element::Br;
        } else {
            return Err(MoleculeError::SmilesParseError(format!(
                "{smi} | invalid char {c}"
            )));
        }
    }

    Ok(Molecule {
        atoms,
        bonds,
        rings: None,
    })
}

fn main() {
    let smi = "C";
    dbg!(smi);
    dbg!(smiles_parser_v1(smi).unwrap());

    let smi = "[18OH-]";
    dbg!(smi);
    dbg!(smiles_parser_v2(smi).unwrap());

    let smi = "CC=C";
    dbg!(smi);
    dbg!(smiles_parser_v3(smi).unwrap());

    let smiles = ["CC(=O)C", "CS(=O)(=O)C", "CC(C(F)F)C"];
    for smi in smiles {
        dbg!(smi);
        dbg!(smiles_parser_v4(smi).unwrap());
    }

    let smiles = ["C1CC1"];
    for smi in smiles {
        dbg!(smi);
        dbg!(smiles_parser_v5(smi).unwrap());
    }
}
