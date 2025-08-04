mod alphabets;

fn main() {
    use markov::Chain;
    use crate::alphabets::alphabets::*;

    let mut chain = Chain::of_order(3);

    // train the chain on some vectors
    chain.feed(vec![Dna4::C, Dna4::C, Dna4::G, Dna4::C])
         .feed(vec![Dna4::A, Dna4::C, Dna4::G, Dna4::T]);
    
    // constants
    let seq_len : usize = 10;

    // generate sequences
    let gen = chain.iter_for(seq_len);
    for c in gen {
        print!("{:?}", &c[0]);
    }
    println!();
}
