use clap::Parser;
use bio::io::fasta::{Reader, Writer};
//use seq_io::fasta::{Reader,Record,Writer};
use std::collections::HashMap;
use std::fmt::Display;
use rand::prelude::*;

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
    seed: u64,

    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

pub fn char_to_int(c : &char) -> u8 {
    u8::try_from(c.clone()).expect("Char out of range")
}

pub fn int_to_char(i : &u8) -> char {
    char::from_u32(*i as u32).expect("can not convert to char").to_ascii_uppercase()
}

pub fn print_record<I: ?Sized>(container : &I, id : impl Display) 
where
    for<'a> &'a I: IntoIterator<Item = &'a u8>,
{
    println!(">{}", id);
    for i in container {
        print!("{}", int_to_char(i));
    }
    println!();
}

fn main() {
    let mut args = Args::parse();

    if args.lens.len() == 0 {
        args.lens.push(1000);
    }

    let reader = Reader::from_file(args.input);
    let records = reader.expect("fasta reader: got an io::Error or could not read_line()").records();

    // hashmap of k-mer and nucleotide requencies
    let mut kmer_counts = HashMap::new();
    let mut char_counts = HashMap::new();

    // markov library
    use markov::Chain;
    let mut chain = Chain::of_order(args.order);

    // manual Markov model
    let mut rng = StdRng::seed_from_u64(args.seed);

    if args.verbose {
        println!("Input FASTA");
    }

    for result in records {
        let record = result.as_ref().expect("Error during fasta record parsing");
        
        if args.verbose {
            print_record(record.seq(), record.id());
        }
        
        // count k-mers
        for i in 0..record.seq().len()-args.order {
            let mut kmer = record.seq()[i..i+args.order].to_vec();
            for i in 0..kmer.len() {
                //print!("{}\t", kmer[i]);
                let mut c = &mut kmer[i];
                //print!("{}\t", int_to_char(c));
                kmer[i] = char_to_int(&mut int_to_char(c)); // ignore case
                //print!("{}\n", kmer[i]);
            }
            
            let c = kmer[0];
            if let Some(count) = char_counts.get_mut(&c) {
                *count = *count + 1;
            }
            else {
                char_counts.insert(c, 1);
            }
            
            if let Some(count) = kmer_counts.get_mut(&kmer) {
                *count = *count + 1;
            }
            else {
                kmer_counts.insert(kmer, 1);
            }
        }
        
        // train the chain on some vectors
        chain.feed(record.seq());
        //TODO: learn Markov chain Übergangswahrscheinlichkeiten from kmer_counts
    }
    
    let mut ref_len : usize = 0;
    
    if args.verbose {
        for (k, n) in &kmer_counts {
            for i in k {
                let c = int_to_char(i);
                print!("{c}");
            }
            print!(":{n}\n");
        }
   
        for (i, n) in &char_counts {
            let c = int_to_char(i);
            //print!("{i}\t");
            print!("{c}");
            print!(":{n}\n");
            ref_len += n;
        }
    }

    // generate sequences
    let mut writer = Writer::to_file(args.output);
    let mut id = 0;

    if args.verbose {
        println!("Output FASTA");
    }
    
    for l in args.lens {
        //TODO: make deterministic using a seed
        let mut iter = chain.iter().scan(args.seed, |_state, x| {return Some(x)});
        let mut rec_out: Vec<u8> = Vec::new();

        // initialize sequence by sampling from char probability distribution
        for _ in 0..args.order {
            let i = rng.random_range(0..ref_len);
            //TODO: make cum probability distribution of chars
            let mut cum_sum : usize = 0;
            for (c, n) in &char_counts {
                cum_sum += n;
                if cum_sum >= i {
                    rec_out.push(*c);
                }
            }
        }
        for _ in args.order..l {
            // Markov library
            let r = iter.next();
            match r {
                Some(c) => rec_out.push(c[0]),
                None    => panic!("Reached iterator end"),
            }

            // std
            let p = rng.random_range(0.0..1.0);
            //TODO: get k-mer that corresponds to probability
            /*
            let mut cum_sum : usize = 0;
            for (k, n) in &kmer_counts {
                cum_sum += n; 
                if cum_sum >= 1 {
                    rec_out.push(*k);
                }
            }*
        }
        
        if args.verbose {
            print_record(&rec_out, id);
        }
        let _ = writer.as_mut().expect("Error writing record").write(&id.to_string(), None, rec_out.as_slice());
        id += 1;
    }
}

