use std::str::FromStr;

use molrs::{
    atom::{Atom, PointChirality},
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

    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();

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

fn main() {
    let smiles = ["C"];
    for smi in smiles {
        dbg!(smi);
        dbg!(smiles_parser_v1(smi).unwrap());
    }

    let smiles = ["C", "[18OH-]"];
    for smi in smiles {
        dbg!(smi);
        dbg!(smiles_parser_v2(smi).unwrap());
    }
}
