use bio::io::fasta::Record;
use std::hash::Hash;
use std::collections::HashMap;
use std::usize;

use crate::args::{SimulateArgs};
use crate::io::{char_to_int, get_annotations, get_records, int_to_char, print_record};

pub trait SequenceModel {
    fn kmer_counts(&self) -> &HashMap<Vec<u8>, usize>;
    fn char_counts(&self) -> &HashMap<u8, usize>;   // 8-bit alphabet
    fn ref_len(&self) -> usize;
}

pub struct MarkovModel {
    pub kmer_counts: HashMap<Vec<u8>, usize>,   // <kmer_seq, count>
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
                    char_counts: &mut HashMap<u8, usize> , ref_len: &mut usize, order: usize) {
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

pub fn find_repeat_positions(annotation: &Vec<u64>, window_pos: &mut Vec<usize>, kmer: u8, pat : &Pattern) {
    let mut max_value : u64 = 0;
    for c in annotation {
        if *c > max_value {
            max_value = *c;
        }
    }

    let mut pos : usize = 0;
    let annot_wind : Vec<&[u64]> = annotation.windows(pat.repeat_len - (kmer as usize) + 1).collect();
    for window in annot_wind {
        if window.iter().all(|&v| v > pat.min_thresh) {
            window_pos.push(pos);
        }
        pos+=1;
    }
}

#[derive(Copy, Clone)]
pub struct Pattern {
    repeat_len : usize,
    min_thresh : u64, 
}
pub struct Rule {
    pat : Pattern,
    seq : Vec<u8>,
    host_len : usize,
}

impl Rule {
    pub fn new(pat : &Pattern, seq : &[u8], pos : usize) -> Self {
        return Rule {pat : *pat, seq : seq[pos..pos + pat.repeat_len].to_vec(), host_len : seq.len()} 
    }
} 

fn get_repeat_score(pat : &Pattern, pos : &Vec<usize>) -> u64 {
    (pat.repeat_len as u64) * pat.min_thresh * (pos.len() as u64)
}

fn find_best_pattern(annotation: &Vec<u64>, window_pos : &mut Vec<usize>, kmer : u8) -> Pattern {
    let mut pat = Pattern{repeat_len : kmer as usize, min_thresh : average(annotation)};
    
    let mut prev_repeat_score = 0;
    let mut repeat_score : u64 = 1;
    while repeat_score > prev_repeat_score {
        // these variables determine how the parameter space is searched
        pat.min_thresh += 1;
        pat.repeat_len += 10;

        prev_repeat_score = repeat_score;
        window_pos.clear();
        find_repeat_positions(&annotation, window_pos, kmer, &pat);

        repeat_score = get_repeat_score(&pat, &window_pos);
    }
    return pat
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

fn average(numbers: &Vec<u64>) -> u64 {
    numbers.iter().sum::<u64>() / numbers.len() as u64
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
        let annotation_records: Result<HashMap<String, Vec<u64>>, std::io::Error> = get_annotations(&args.annotation.clone());
        
        let mut grammar = Vec::new();
        while let Some(Ok(seq_record)) = seq_records.next() {
            let annotation : &Vec<u64> = annotation_records.as_ref().expect("Error during GenMap record parsing")[seq_record.id()].as_ref();
            // ref_len += seq_record.seq().len();
            // print!("{}", record.seq()[0])
            if args.verbose {
                print_record(seq_record.seq(), seq_record.id());
            }

            let mut window_pos = Vec::new();
            let pat = find_best_pattern(annotation, &mut window_pos, args.kmer);

            for pos in window_pos {
                let rule = Rule::new(&pat, seq_record.seq(), pos);
                grammar.push(rule);    
            }

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

            if !grammar.is_empty() {
                println!("host_len\trepeat_len\tmin_thresh\tseq");
            }
            for rule in &grammar {
                print!("{}\t", rule.host_len);
                print!("{}\t", rule.pat.repeat_len);
                print!("{}\t", rule.pat.min_thresh);
                for i in &rule.seq {
                    let c = int_to_char(i);
                    print!("{c}");
                }
                println!();
            }
             
        }

        let markov_model = MarkovModel {kmer_counts: kmer_counts, char_counts: char_counts, ref_len: ref_len };
        return SequenceGrammarModel { markov_model: markov_model, grammar: grammar };
    }
}
