pub mod alphabets {
    
use std::hash::{Hash, Hasher};
use std::fmt::{Debug};

#[derive(Debug)]
pub enum Dna4 {
   A,
   C,
   G,
   T,
}

impl Clone for Dna4 {
    fn clone(&self) -> Dna4 {
        match self {
            Dna4::A => Dna4::A,
            Dna4::C => Dna4::C,
            Dna4::G => Dna4::G,
            Dna4::T => Dna4::T,
        }
    }
}

impl Hash for Dna4 {
    fn hash<H : Hasher>(&self, state: &mut H) {
        let c = match self {
            Dna4::A => 0,
            Dna4::C => 1,
            Dna4::G => 2,
            Dna4::T => 3,
        };
        c.hash(state);
    }
} 

impl PartialEq for Dna4 {
    fn eq(&self, other: &Self) -> bool {
        //return self.hash() == other.hash()
        match (self, other) {
            (Dna4::A, Dna4::A) => true,
            (Dna4::C, Dna4::C) => true,
            (Dna4::G, Dna4::G) => true,
            (Dna4::T, Dna4::T) => true,
            _ => false,
        }
    }
}

// full equality requires reflexivity which can not be checked by the compiler
impl Eq for Dna4 {}

}

