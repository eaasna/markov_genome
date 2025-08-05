use clap::Parser;
use bio::io::fasta::{Reader, Writer};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[arg(short, long)]
    input: String,

    #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    output: String,
    
    #[arg(short, long)]
    lens: Vec<usize>,

    #[arg(long, default_value_t = 3)]
    order: usize,

    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

pub fn print_dna4_vec(v : &Vec<Vec<u8>>) {
    for i in v {
        for j in i {
            let c = match j {
                b'A' | b'a' => 'A',
                b'C' | b'c' => 'C',
                b'G' | b'g' => 'G',
                b'T' | b't' => 'T',
                _ => panic!("Unknown char\t{:b}", j),
            //TODO: read char from int code
            };
            print!("{}", c);
        }
        println!();
    }
}
fn main() {
    let mut args = Args::parse();

    if args.lens.len() == 0 {
        args.lens.push(1000);
    }
    let reader = Reader::from_file(args.input);
    let records = reader.expect("fasta reader: got an io::Error or could not read_line()").records();

    let mut rec_in_vec : Vec<Vec<u8>> = Vec::new(); 
    for result in records {
        let record = result.expect("Error during fasta record parsing");
        rec_in_vec.push(record.seq().to_vec());
    }
 
    if args.verbose {
        println!("Input FASTA");
        print_dna4_vec(&rec_in_vec);
    }

    use markov::Chain;
    let mut chain = Chain::of_order(args.order);
    // train the chain on some vectors
    for i in 0..rec_in_vec.len() {
        chain.feed(rec_in_vec[i].clone());
    }

    // generate sequences
    let mut rec_out_vec : Vec<Vec<u8>> = Vec::new();

    for l in args.lens {
        let mut gen = chain.iter();
        let mut rec_out : Vec<u8> = Vec::new();
        for _ in 0..l {
            let r = gen.next();
            match r {
                Some(c) => rec_out.push(c[0]),
                None    => panic!("Reached iterator end"),
            }
        }
        rec_out_vec.push(rec_out.clone());
    }
   
    let mut writer = Writer::to_file(args.output);
    let mut id = 0;
    for seq in &rec_out_vec {
        let _ = writer.as_mut().expect("Error writing record").write(&id.to_string(), None, seq.as_slice());
        id += 1;
    }
    
    if args.verbose {
        println!("Output FASTA");
        print_dna4_vec(&rec_out_vec);
    }
}

