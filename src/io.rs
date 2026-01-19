use bio::io::fasta::{Reader};
use std::collections::HashMap;
use std::fmt::Display;
use std::fs::File;
use std::io::{self, Error, ErrorKind, BufRead, BufReader};

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

/// Read records from FASTA format.
pub fn get_records(input : String) -> bio::io::fasta::Records<BufReader<File>> { 
    let reader = Reader::from_file(input.clone());
    reader.expect("fasta reader: got an io::Error or could not read_line()").records()
}

// https://users.rust-lang.org/t/reading-integers-from-a-file-into-vector/17517/5
pub fn get_annotations(input: &str) -> Result<HashMap<String, Vec<u64>>, io::Error> {
    let fin = File::open(input)?;
    let br = BufReader::new(fin);
    let mut records  = HashMap::new();
    
    let mut record_id = String::new();
    for line in br.lines() {
        let line = line?;
        if line.starts_with(">") {
            record_id = line[1..].to_string();
        } else {
            if record_id.is_empty() {
                panic!{"Missing genmap record id"};
            }
            records.insert(record_id, parse_annotation_line(&line)?);
            record_id = String::new();
        }
    }
    Ok(records)
}

fn parse_annotation_line(line: &String) -> Result<Vec<u64>, io::Error> {
    let mut v = vec![];
    for annot in line.split(' ').collect::<Vec<&str>>() {
        v.push(annot.trim()
                    .parse::<u64>()
                    .map_err(|e| Error::new(ErrorKind::InvalidData, e))?); 
    }
    Ok(v)
}
