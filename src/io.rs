use bio::io::fasta::{Reader};
use std::fmt::Display;
use std::io::BufReader;
use std::fs::File;

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

pub fn get_records(input : String) -> bio::io::fasta::Records<BufReader<File>> { 
    let reader = Reader::from_file(input.clone());
    reader.expect("fasta reader: got an io::Error or could not read_line()").records()
}

