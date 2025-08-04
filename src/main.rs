mod alphabets;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[arg(short, long)]
    input: String,

    #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    output: String,
    
    #[arg(short, long, default_value_t = 100)]
    seq_len: usize,
}

fn main() {
    use crate::alphabets::alphabets::*;
    
    let args = Args::parse();

    use markov::Chain;
    let mut chain = Chain::of_order(3);

    // train the chain on some vectors
    chain.feed(vec![Dna4::C, Dna4::C, Dna4::G, Dna4::C])
         .feed(vec![Dna4::A, Dna4::C, Dna4::G, Dna4::T]);
    
    // generate sequences
    let gen = chain.iter_for(args.seq_len);
    for c in gen {
        print!("{:?}", &c[0]);
    }
    println!();
}

