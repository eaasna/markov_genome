use bio::io::fasta::{Writer};
use rand::prelude::*;

use crate::sequence_model::{build_markov_model};
use crate::args::SimulateArgs;
use crate::io::{print_record};

pub fn run_markov_simulation(args : &SimulateArgs) {
    let sequence_model = build_markov_model(args);
    // seed for reproducible results
    let mut rng = StdRng::seed_from_u64(args.seed);

    // generate sequences
    let mut writer = Writer::to_file(args.output.clone());
    let mut id = 0;

    if args.verbose {
        println!("Output FASTA");
    }
    
    let mut alphabet: Vec<u8> = sequence_model.char_counts.clone().into_keys().collect();
    alphabet.sort_unstable(); // make deterministic
    
    for l in &args.lens {
        let mut rec_out: Vec<u8> = Vec::new();

        // initialize sequence by sampling from char probability distribution
        for _ in 0..args.order-1 {
            let i = rng.random_range(0..sequence_model.ref_len);
            //TODO: make cum probability distribution of chars
            // hold in memory instead of recalculating each time
            let mut cum_sum : usize = 0;
            for c in &alphabet {
                if let Some(n) = sequence_model.char_counts.get(c) {
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
                if let Some(count) = sequence_model.kmer_counts.get(&prev_states) {
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

