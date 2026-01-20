use bio::io::fasta::{Writer};
use rand::prelude::*;

use crate::args::SimulateArgs;
use crate::sequence_model::{SequenceModel};
use crate::io::{print_record};

fn simulate_markov(sequence_model : & impl SequenceModel, len : usize, seed : u64, rec_out : &mut Vec<u8>) {
    // seed for reproducible results
    let mut rng: StdRng = StdRng::seed_from_u64(seed);

    // initialize sequence by sampling from char probability distribution
    for _ in 0..sequence_model.order()-1 {
        let i = rng.random_range(0..sequence_model.ref_len());
        //TODO: make cum probability distribution of chars
        // hold in memory instead of recalculating each time
        let mut cum_sum : usize = 0;
        for c in sequence_model.alphabet() {
            if let Some(n) = sequence_model.char_counts().get(&c) {
                cum_sum += n;
                if cum_sum >= i {
                    rec_out.push(*c);
                    break;
                }
            }
        }
    }      
    assert_eq!(rec_out.len() + 1, sequence_model.order() as usize);

    // walk through Markov chain
    for _ in ((sequence_model.order() as usize) -1)..len {
        let mut prev_states = Vec::from_iter(rec_out[(rec_out.len() + 1 - (sequence_model.order() as usize))..rec_out.len()].iter().cloned());
        prev_states.push(sequence_model.alphabet()[0]);
        assert_eq!(prev_states.len(), sequence_model.order() as usize);

        // for some state e.g. AC gather the occurrence counts of k-mers ACA, ACC, ACG, ACT
        // then normalize to find the transition probabilities
        let mut state_sum : usize = 0;
        let mut next_count : Vec<usize> = Vec::new();
        for next in sequence_model.alphabet() {
            prev_states[(sequence_model.order() as usize) - 1] = *next;
            if let Some(count) = sequence_model.kmer_counts().get(&prev_states) {
                next_count.push(*count);
                state_sum += *count;
            }
        }
        
        if next_count.len() == 0 {
            // avoid terminating early by defaulting to uniform transition probabilities
            for _ in sequence_model.alphabet() {
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
                rec_out.push(sequence_model.alphabet()[i]);
                break;
            }
        }
    }    
}

pub fn run_simulation(sequence_model : & impl SequenceModel, args: &SimulateArgs) {
    // generate sequences
    let mut writer = Writer::to_file(args.output.clone());
    let mut id = 0;

    let mut rng: StdRng = StdRng::seed_from_u64(args.seed);
    for len in &args.lens {
        let mut rec_out: Vec<u8> = Vec::new();
        simulate_markov(sequence_model, *len, args.seed, &mut rec_out);
        if args.verbose {
            print_record(&rec_out, id);
        }
        assert_eq!(rec_out.len(), *len);

        sequence_model.apply_grammar(&mut rng, &mut rec_out);
        
        let _ = writer.as_mut().expect("Error writing record").write(&id.to_string(), None, rec_out.as_slice());
        id += 1;
    }
}

