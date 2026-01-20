use bio::io::fasta::Record;
use rand::rngs::StdRng;
use std::hash::Hash;
use std::collections::HashMap;
use std::usize;

use crate::grammar::{Rule, find_best_pattern};
use crate::args::{SimulateArgs};
use crate::io::{char_to_int, get_annotations, get_records, int_to_char};

pub trait SequenceModel {
    fn kmer_counts(&self) -> &HashMap<Vec<u8>, usize>;
    fn char_counts(&self) -> &HashMap<u8, usize>;   // 8-bit alphabet
    fn ref_len(&self) -> usize;
    fn order(&self) -> u8;
    fn alphabet(&self) -> &Vec<u8>;
    fn apply_grammar(&self, rng : &mut StdRng, seq : &mut Vec<u8>);
}

pub struct MarkovModel {
    pub kmer_counts: HashMap<Vec<u8>, usize>,   // <kmer_seq, count>
    pub char_counts: HashMap<u8, usize>,
    pub ref_len: usize,
    pub order: u8,
    pub alphabet: Vec<u8>, 
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

    fn order(&self) -> u8 {
        self.order
    }

    fn alphabet(&self) -> &Vec<u8> {
        &self.alphabet
    }

    fn apply_grammar(&self, _rng : &mut StdRng, _seq : &mut Vec<u8>) {
        // do nothing
    }
}

pub struct SequenceGrammarModel {
    markov_model: MarkovModel,
    pub grammar: Vec<Rule>,
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

    fn order(&self) -> u8 {
        self.markov_model.order
    }

    fn alphabet(&self) -> &Vec<u8> {
        self.markov_model.alphabet()
    }

    fn apply_grammar(&self, rng : &mut StdRng, seq : &mut Vec<u8>) {
        for rule in &self.grammar {
            rule.apply(seq, rng);
        }
    }
}

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

pub fn count_record(record : &Record, kmer_counts: &mut HashMap<Vec<u8>, usize>, 
                    char_counts: &mut HashMap<u8, usize> , ref_len: &mut usize, order: u8) {
    for i in 0..record.seq().len()-(order as usize) {
        let mut kmer = record.seq()[i..i+(order as usize)].to_vec();
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

impl MarkovModel {
    pub fn from_counted(kmer_counts : &HashMap<Vec<u8>, usize> , char_counts : &HashMap<u8, usize>, ref_len : usize, order : u8) -> Self {
        let mut alphabet: Vec<u8> = char_counts.clone().into_keys().collect();
        alphabet.sort_unstable(); // make deterministic

        return MarkovModel {kmer_counts: kmer_counts.clone(), char_counts: char_counts.clone(), ref_len: ref_len, order : order, alphabet : alphabet };
    }

    pub fn new(args : &SimulateArgs) -> Self {
        // hashmaps of k-mer and nucleotide frequencies
        //TODO: initialize outside the constructor so no clone()
        let mut kmer_counts = HashMap::new();
        let mut char_counts = HashMap::new();

        // learn Markov probabilities
        let mut ref_len : usize = 0;
        for result in get_records(args.input.clone()) {
            let record = result.as_ref().expect("Error during fasta record parsing");
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

        return MarkovModel::from_counted(&kmer_counts, &char_counts, ref_len, args.order)
    }

}

impl SequenceGrammarModel {
    pub fn new(args : &SimulateArgs) -> Self {
    // hashmaps of k-mer and nucleotide frequencies
        let mut kmer_counts = HashMap::new();
        let mut char_counts = HashMap::new();

        // learn Markov probabilities
        let mut ref_len : usize = 0;
        let mut seq_records: bio::io::fasta::Records<std::io::BufReader<std::fs::File>> = get_records(args.input.clone());
        let annotation_records: Result<HashMap<String, Vec<u64>>, std::io::Error> = get_annotations(&args.annotation.clone());
        
        let mut grammar = Vec::new();
        
        //TODO: gather repeats across the whole reference not per sequence
        while let Some(Ok(seq_record)) = seq_records.next() {
            count_record(&seq_record, &mut kmer_counts, &mut char_counts, &mut ref_len, args.order);

            let annotation : &Vec<u64> = annotation_records.as_ref().expect("Error during GenMap record parsing")[seq_record.id()].as_ref();
            let mut window_pos = Vec::new();
            let pat = find_best_pattern(annotation, &mut window_pos, args.kmer);

            for pos in window_pos {
                let rule = Rule::new(&pat, seq_record.seq(), pos, ref_len);
                grammar.push(rule);
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

            if !grammar.is_empty() {
                println!("ref_len\trepeat_len\tmin_thresh\trepeat_seq");
            }
            for rule in &grammar {
                print!("{}\t", rule.ref_len);
                print!("{}\t", rule.pat.repeat_len);
                print!("{}\t", rule.pat.min_thresh);
                for i in &rule.repeat_seq {
                    let c = int_to_char(i);
                    print!("{c}");
                }
                println!();
            }
             
        }

        let markov_model = MarkovModel::from_counted(&kmer_counts, &char_counts, ref_len, args.order);
        return SequenceGrammarModel { markov_model: markov_model, grammar: grammar };
    }
}
