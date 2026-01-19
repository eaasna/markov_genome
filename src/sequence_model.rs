use std::collections::HashMap;
use std::hash::Hash;

use crate::args::SimulateArgs;
use crate::io::{get_records, print_record, char_to_int, int_to_char};

/// Iterate the value for the given key by +1 if the key exists or add a new key with value 1 if the key does not exist.
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

pub struct SequenceModel {
    pub kmer_counts: HashMap<Vec<u8>, usize>,
    pub char_counts: HashMap<u8, usize>,
    pub ref_len: usize,
}   

impl SequenceModel {
    pub fn new(args : &SimulateArgs) -> Self {
    // hashmaps of k-mer and nucleotide frequencies
        let mut kmer_counts = HashMap::new();
        let mut char_counts = HashMap::new();

        if args.verbose {
            println!("Input FASTA");
        }

        // learn Markov probabilities
        let mut ref_len : usize = 0;
        for result in get_records(args.input.clone()) {
            let record = result.as_ref().expect("Error during fasta record parsing");
            
            if args.verbose {
                print_record(record.seq(), record.id());
            }
            
            for i in 0..record.seq().len()-args.order {
                let mut kmer = record.seq()[i..i+args.order].to_vec();
                for i in 0..kmer.len() {
                    let c = &mut kmer[i];
                    kmer[i] = char_to_int(&mut int_to_char(c)); // ignore case
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

        return SequenceModel { kmer_counts: kmer_counts, char_counts: char_counts, ref_len: ref_len };
    }
}
