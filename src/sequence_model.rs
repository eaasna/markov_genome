use std::collections::HashMap;
use std::hash::Hash;
use std::usize;
use bio::io::fasta::Record;

use crate::args::{SimulateArgs};
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

pub fn count_record(record : &Record, kmer_counts: &mut HashMap<Vec<u8>, usize>, char_counts: &mut HashMap<u8, usize> , ref_len: &mut usize, order: usize) {
    for i in 0..record.seq().len()-order {
        let mut kmer = record.seq()[i..i+order].to_vec();
        for i in 0..kmer.len() {
            let c = &mut kmer[i];
            kmer[i] = char_to_int(&mut int_to_char(c)); // ignore case
        }
        
        let c = kmer[0];
        // count chars and k-mers
        update_count_map( kmer_counts, kmer);
        update_count_map( char_counts, c);
        *ref_len += 1;
    }
}

pub trait SequenceModel {
    fn kmer_counts(&self) -> &HashMap<Vec<u8>, usize>;
    fn char_counts(&self) -> &HashMap<u8, usize>;
    fn ref_len(&self) -> usize;
}

pub struct MarkovModel {
    pub kmer_counts: HashMap<Vec<u8>, usize>,
    pub char_counts: HashMap<u8, usize>,
    pub ref_len: usize,
}   

impl SequenceModel for MarkovModel  {
    fn kmer_counts(&self) -> &HashMap<Vec<u8>, usize> {
        &self.kmer_counts
    }

    fn char_counts(&self) -> &HashMap<u8, usize> {
        &self.char_counts
    }

    fn ref_len(&self) -> usize {
        self.ref_len
    }
}

pub struct SequenceGrammarModel {
    markov_model: MarkovModel,
    pub repeat_counts: HashMap<Vec<u8>, usize>,
}

impl SequenceModel for SequenceGrammarModel {
    fn kmer_counts(&self) -> &HashMap<Vec<u8>, usize> {
        self.markov_model.kmer_counts()
    }

    fn char_counts(&self) -> &HashMap<u8, usize> {
        self.markov_model.char_counts()
    }

    fn ref_len(&self) -> usize {
        self.markov_model.ref_len
    }
}

impl MarkovModel {
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

            count_record(record, &mut kmer_counts, &mut char_counts, &mut ref_len, args.order);
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

        return MarkovModel {kmer_counts: kmer_counts, char_counts: char_counts, ref_len: ref_len };
    }
}

impl SequenceGrammarModel {
    pub fn new(args : &SimulateArgs) -> Self {
    // hashmaps of k-mer and nucleotide frequencies
        let mut kmer_counts = HashMap::new();
        let mut char_counts = HashMap::new();

        if args.verbose {
            println!("Input FASTA");
        }

        // learn Markov probabilities
        let mut ref_len : usize = 0;
        let mut seq_records = get_records(args.input.clone());
        let mut annotation_records = get_records(args.annotation.clone());
        
        while let Some(Ok(seq_record)) = seq_records.next() {
            let annotation_record = annotation_records.next().expect("No next").expect("Error during GenMap record parsing");
            // ref_len += seq_record.seq().len();
            // print!("{}", record.seq()[0])
            if args.verbose {
                print_record(seq_record.seq(), seq_record.id());
            }
            
            print!("{}", annotation_record.seq().len());
            count_record(&seq_record, &mut kmer_counts, &mut char_counts, &mut ref_len, args.order);
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

        let markov_model = MarkovModel {kmer_counts: kmer_counts, char_counts: char_counts, ref_len: ref_len };
        return SequenceGrammarModel { markov_model: markov_model, repeat_counts: HashMap::new() };
    }
}
