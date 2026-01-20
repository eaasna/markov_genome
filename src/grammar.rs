use rand::{rngs::StdRng, Rng};

use crate::io::{char_to_int, int_to_char};

#[derive(Copy, Clone)]
pub struct Pattern {
    pub repeat_len : usize,
    pub min_thresh : u64, 
}
pub struct Rule {
    pub pat : Pattern,
    pub repeat_seq : Vec<u8>,
    pub ref_len : usize,
}

impl Rule {
    pub fn new(pat : &Pattern, seq : &[u8], pos : usize, ref_len : usize) -> Self {
        let mut repeat_vec = Vec::new();
        for i in 0..pat.repeat_len {
            let mut c = seq[pos + i];
            repeat_vec.push(char_to_int(&mut int_to_char(&mut c)));  // ignore case
        }

        return Rule {pat : *pat, repeat_seq : repeat_vec, ref_len : ref_len} 
    }

    pub fn apply(&self, seq : &mut Vec<u8>, rng : &mut StdRng) {        
        assert!(self.pat.repeat_len == self.repeat_seq.len());
        assert!(self.repeat_seq.len() < seq.len());

        let mut repeat_positions = Vec::new();
        let mut count = 0;
        let repeats_to_simulate = self.pat.min_thresh / (self.ref_len * seq.len()) as u64;
        while count <= repeats_to_simulate {
            let i = rng.random_range(0..seq.len() - self.pat.repeat_len);
            repeat_positions.push(i);
            count += 1;
        }

        repeat_positions.sort();

        let mut offset = 0;
        for pos in repeat_positions {
            while offset < self.repeat_seq.len() - 1 {
                seq[pos + offset] = self.repeat_seq[offset];
                offset += 1;
            }
        }
    }
}

fn get_repeat_score(pat : &Pattern, pos : &Vec<usize>) -> u64 {
    (pat.repeat_len as u64) * pat.min_thresh * (pos.len() as u64)
}

fn average(numbers: &Vec<u64>) -> u64 {
    numbers.iter().sum::<u64>() / numbers.len() as u64
}

fn find_repeat_positions(annotation: &Vec<u64>, window_pos: &mut Vec<usize>, kmer: u8, pat : &Pattern) {
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

pub fn find_best_pattern(annotation: &Vec<u64>, window_pos : &mut Vec<usize>, kmer : u8) -> Pattern {
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
