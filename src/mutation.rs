use bio::io::fasta::{Writer};
use rand::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::args::MutateArgs;
use crate::io::{get_records, print_record, char_to_int, int_to_char};

pub fn run_mutation(args : &MutateArgs) {
    if args.verbose {
        println!("Input FASTA");
    }

    let mut alphabet_map = HashMap::new();
    let mut rank_map = HashMap::new();
    {
        let mut chars = HashSet::new();
        for result in get_records(args.input.clone()) {
            let record = result.as_ref().expect("Error during fasta record parsing");
            for c in record.seq() {
                chars.insert(char_to_int(&mut int_to_char(c))); // ignore case
            }
        }
        let mut rank = 0;
        for c in chars.iter() {
            alphabet_map.insert(c.clone(), rank);
            rank_map.insert(rank, c.clone());
            rank += 1;
        }
    }
    let sigma = alphabet_map.len();

    let mut writer = Writer::to_file(args.output.clone());
    for result in get_records(args.input.clone()) {
        let record = result.as_ref().expect("Error during fasta record parsing");

        if args.verbose {
            print_record(record.seq(), record.id());
        }
       
        let mut rng = StdRng::seed_from_u64(args.seed);
        let mut error_count: usize = 0;
       
        let mut rec_out: Vec<u8> = Vec::new();
        for i in 0..record.seq().len() {
            let mutation_state = rng.random_range(0.0..1.0);
            if mutation_state <= args.error {
                let var_state = rng.random_range(1..sigma); // for DNA4 there are 3 possible variants
                //!TODO: nothing is written if not DNA4 alphabet with only capital letters 
                if let Some(curr_rank) = alphabet_map.get(&record.seq()[i]) {
                    let mutated_rank = (curr_rank + var_state) % sigma;
                    if let Some(var_char) = rank_map.get(&mutated_rank) {
                        rec_out.push(*var_char);
                        error_count += 1;
                        /*
                        println!("{}", record.seq()[i].to_string() + "\t" + &var_char.to_string());
                        println!("{}", curr_rank.to_string() + "\t" + &mutated_rank.to_string());
                        */
                        assert!(record.seq()[i] != *var_char);
                    }
                }
            } else {
                rec_out.push(char_to_int(&mut int_to_char(&record.seq()[i])));
            }
        }
        let id = record.id().to_owned() + "_e" + &error_count.to_string();

        let _ = writer.as_mut().expect("Error writing record").write(&id.to_string(), None, rec_out.as_slice());
    }
}

