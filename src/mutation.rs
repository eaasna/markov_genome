use bio::io::fasta::{Writer};

use crate::args::MutateArgs;
use crate::io::{get_records, print_record};

pub fn run_mutation(args : &MutateArgs) {
    if args.verbose {
        println!("Input FASTA");
    }

    for result in get_records(args.input.clone()) {
        let record = result.as_ref().expect("Error during fasta record parsing");

        if args.verbose {
            print_record(record.seq(), record.id());
        }
    }
}

