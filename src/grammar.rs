#[derive(Copy, Clone)]
pub struct Pattern {
    pub repeat_len : usize,
    pub min_thresh : u64, 
}
pub struct Rule {
    pub pat : Pattern,
    pub seq : Vec<u8>,
    pub host_len : usize,
}

impl Rule {
    pub fn new(pat : &Pattern, seq : &[u8], pos : usize) -> Self {
        return Rule {pat : *pat, seq : seq[pos..pos + pat.repeat_len].to_vec(), host_len : seq.len()} 
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
