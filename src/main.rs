use clap::Parser;
use bio::io::fasta::{Reader, Writer};
//use seq_io::fasta::{Reader,Record,Writer};
use std::collections::HashMap;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[arg(short, long)]
    input: String,

    #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    output: String,
    
    #[arg(short, long)]
    lens: Vec<usize>,

    #[arg(long, default_value_t = 3)]
    order: usize,

    #[arg(long, default_value_t = 42)]
    seed: usize,

    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

pub fn print_dna4_vec(v : &Vec<Vec<char>>) {
    for i in v {
        for j in i {
            print!("{}", j);
        }
        println!();
    }
}

pub fn char_to_u8(c : &char) -> u8 {
    u8::try_from(c.clone()).expect("Char out of range")
}

pub fn u8_to_char(i : &u8) -> char {
    char::from_u32(*i as u32).expect("can not convert to char").to_ascii_uppercase()
}

/*
pub fn print_dna4_vec(v : &Vec<Vec<u8>>) {
    for i in v {
        for j in i {
            let c = match j {
                b'A' | b'a' => 'A',
                b'C' | b'c' => 'C',
                b'G' | b'g' => 'G',
                b'T' | b't' => 'T',
                _ => panic!("Unknown char\t{:b}", j),
            };
            print!("{}", j);
        }
        println!();
    }
}*/

fn main() {
    let mut args = Args::parse();

    if args.lens.len() == 0 {
        args.lens.push(1000);
    }
    let reader = Reader::from_file(args.input);
    let records = reader.expect("fasta reader: got an io::Error or could not read_line()").records();

    // hashmap of k-mer frequencies
    let mut kmer_counts = HashMap::new();

    let mut rec_in_vec : Vec<Vec<char>> = Vec::new(); 
    for result in records {
        let record = result.expect("Error during fasta record parsing");
        
        let mut char_seq : Vec<char> = Vec::new();
        for i in record.seq() {
            char_seq.push(u8_to_char(i));
        }        
        rec_in_vec.push(char_seq.clone());
        // count k-mers
        
        for i in args.order..char_seq.len()-3 {
            let kmer = char_seq[i..i+3].to_vec();
            if let Some(x) = kmer_counts.get_mut(&kmer) {
                *x = *x + 1;
            }
            else {
                kmer_counts.insert(kmer, 1);
            }
        }
    } 

    for (k, c) in &kmer_counts {
        for i in k {
            print!("{i}");
        }
        print!(":{c}\n");
    }

    if args.verbose {
        println!("Input FASTA");
        print_dna4_vec(&rec_in_vec);
    }
    
    // markov library
    use markov::Chain;
    let mut chain = Chain::of_order(args.order);
    // train the chain on some vectors
    for i in 0..rec_in_vec.len() {
        chain.feed(rec_in_vec[i].clone());
    }

    // generate sequences
    let mut rec_out_vec : Vec<Vec<char>> = Vec::new();

    for l in args.lens {
        //TODO: make deterministic using a seed
        let mut iter = chain.iter().scan(args.seed, |_state, x| {return Some(x)});
        let mut rec_out : Vec<char> = Vec::new();
        for _ in 0..l {
            let r = iter.next();
            match r {
                Some(c) => rec_out.push(c[0]),
                None    => panic!("Reached iterator end"),
            }
        }
        rec_out_vec.push(rec_out.clone());
    }
   
    let mut writer = Writer::to_file(args.output);
    let mut id = 0;
    for char_seq in &rec_out_vec {
        let mut int_seq : Vec<u8> = Vec::new();
        
        for c in char_seq {
            int_seq.push(char_to_u8(c));
        }
        
        let _ = writer.as_mut().expect("Error writing record").write(&id.to_string(), None, int_seq.as_slice());
        id += 1;
    }
    
    if args.verbose {
        println!("Output FASTA");
        print_dna4_vec(&rec_out_vec);
    }
}

