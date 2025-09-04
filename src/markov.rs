use bio::io::fasta::{Reader, Writer};
use std::collections::HashMap;
use std::fmt::Display;
use rand::prelude::*;
use std::hash::Hash;

use crate::args::SimulateArgs;

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

pub fn update_count_map<K>(map : &mut HashMap<K, usize>, key : K)
where K: Eq, K: Hash
{ 
    // count k-mers
    if let Some(count) = map.get_mut(&key) {
        *count = *count + 1;
    }
    else {
        map.insert(key, 1);
    }
}

pub fn run_markov_simulation(args : &SimulateArgs) {
    // hashmaps of k-mer and nucleotide frequencies
    let mut kmer_counts = HashMap::new();
    let mut char_counts = HashMap::new();

    // seed for reproducible results
    let mut rng = StdRng::seed_from_u64(args.seed);

    if args.verbose {
        println!("Input FASTA");
    }

    let reader = Reader::from_file(args.input.clone());
    let records = reader.expect("fasta reader: got an io::Error or could not read_line()").records();

    // learn Markov probabilities
    let mut ref_len : usize = 0;
    for result in records {
        let record = result.as_ref().expect("Error during fasta record parsing");
        
        if args.verbose {
            print_record(record.seq(), record.id());
        }
        
        for i in 0..record.seq().len()-args.order {
            let mut kmer = record.seq()[i..i+args.order].to_vec();
            for i in 0..kmer.len() {
                let c = &mut kmer[i];
                kmer[i] = char_to_int(&mut int_to_char(c)); // ignore case
                /*
                print!("{}\t", kmer[i]);
                print!("{}\t", int_to_char(c));
                print!("{}\n", kmer[i]);
                */
            }
            
            let c = kmer[0];
            // count chars and k-mers
            update_count_map(&mut char_counts, c);
            update_count_map(&mut kmer_counts, kmer);
            ref_len += 1;
        }
    }
    
    // test
    if args.verbose {
        let mut kmer_count_total = 0;
        for (k, n) in &kmer_counts {
            kmer_count_total += n;
            for i in k {
                let c = int_to_char(i);
                print!("{c}");
            }
            print!(":{n}\n");
        }
        assert_eq!(ref_len, kmer_count_total);
   
        let mut char_count_total = 0;
        for (i, n) in &char_counts {
            char_count_total += n;
            let c = int_to_char(i);
            //print!("{i}\t");
            print!("{c}");
            print!(":{n}\n");
        }
        assert_eq!(ref_len, char_count_total);
    }

    // generate sequences
    let mut writer = Writer::to_file(args.output.clone());
    let mut id = 0;

    if args.verbose {
        println!("Output FASTA");
    }
    
    let mut alphabet: Vec<u8> = char_counts.clone().into_keys().collect();
    alphabet.sort_unstable(); // make deterministic
    
    for l in &args.lens {
        let mut rec_out: Vec<u8> = Vec::new();

        // initialize sequence by sampling from char probability distribution
        for _ in 0..args.order-1 {
            let i = rng.random_range(0..ref_len);
            //TODO: make cum probability distribution of chars
            // hold in memory instead of recalculating each time
            let mut cum_sum : usize = 0;
            for c in &alphabet {
                if let Some(n) = char_counts.get(c) {
                    cum_sum += n;
                    if cum_sum >= i {
                        rec_out.push(*c);
                        break;
                    }
                }
            }
        }      
        assert_eq!(rec_out.len(), args.order - 1);

        // walk through Markov chain
        for _ in args.order-1..*l {
            let mut prev_states = Vec::from_iter(rec_out[(rec_out.len() + 1 - args.order)..rec_out.len()].iter().cloned());
            prev_states.push(alphabet[0]);
            assert_eq!(prev_states.len(), args.order);

            // for some state e.g. AC gather the occurrence counts of k-mers ACA, ACC, ACG, ACT
            // then normalize to find the transition probabilities
            let mut state_sum : usize = 0;
            let mut next_count : Vec<usize> = Vec::new();
            for next in &alphabet {
                prev_states[args.order - 1] = *next;
                if let Some(count) = kmer_counts.get(&prev_states) {
                    next_count.push(*count);
                    state_sum += *count;
                }
            }
            
            if next_count.len() == 0 {
                // avoid terminating early by defaulting to uniform transition probabilities
                for _ in &alphabet {
                    next_count.push(1);
                    state_sum += 1;
                }
            }

            // apply decision border from random probability
            let p = rng.random_range(0.0..1.0);
            let mut cum_sum : usize = 0;
            let decision_border = ((state_sum as f64)* p) as usize;
            for i in 0..next_count.len() {
                cum_sum += next_count[i];
                if cum_sum >= decision_border {
                    rec_out.push(alphabet[i]);
                    break;
                }
            }
        }
        
        if args.verbose {
            print_record(&rec_out, id);
        }

        assert_eq!(rec_out.len(), *l);
        let _ = writer.as_mut().expect("Error writing record").write(&id.to_string(), None, rec_out.as_slice());
        id += 1;
    }
}

