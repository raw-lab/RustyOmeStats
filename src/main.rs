use bio::{alignment::pairwise::Aligner, io::fasta};
use bitvec::prelude as bv;
use clap::Parser;
use regex::Regex;
use std::{collections::HashMap, error::Error, fs, path::PathBuf};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Interval size in number of residues
    #[arg(short, long)]
    interval: usize,

    /// Fasta file or folder
    #[arg(short, long)]
    fasta: PathBuf,

    /// Reference genome
    #[arg(short, long)]
    reference: Option<PathBuf>,
}

struct MaskArray {
    data: bv::BitVec,
}

impl MaskArray {
    fn new(size: usize) -> Self {
        MaskArray {
            data: bv::bitvec![0; size],
        }
    }

    fn apply_mask(&mut self, start: usize, stop: usize) {
        for i in start..stop {
            self.data.set(i, true);
        }
    }

    fn to_string(&self) -> String {
        /* self.data
        .iter()
        .map(|b| if *b { '1' } else { '0' })
        .collect() */
        self.data.to_string()
    }
}

fn validate_input_dir(fasta: PathBuf) -> Vec<PathBuf> {
    let files = if fasta.is_file() {
        vec![fasta.clone()]
    } else if fasta.is_dir() {
        fs::read_dir(&fasta)
            .unwrap()
            .filter_map(|entry| {
                let path = entry.ok()?.path();
                if path.is_file()
                    && path
                        .extension()
                        .and_then(|ext| ext.to_str())
                        .map_or(false, |ext| {
                            [".fasta", ".fa", ".fna", ".ffn"].contains(&ext)
                        })
                {
                    Some(path)
                } else {
                    None
                }
            })
            .collect()
    } else {
        panic!("Not a valid fasta file or directory")
    };

    files
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let files = validate_input_dir(args.fasta.clone());
    let mut outputs: Vec<String> = vec![];

    let mut gc_count = 0;
    let mut at_count = 0;
    let mut length_inter: HashMap<usize, i32> = HashMap::new();
    let mut seq_lengths: Vec<usize> = vec![];

    let mut all_records: Vec<fasta::Record> = Vec::new();
    for file in files {
        let reader = fasta::Reader::from_file(file)?;
        let records: Vec<fasta::Record> = reader.records().map(|r| r.unwrap()).collect();
        all_records.extend(records);
    }
    let num_seq = all_records.len();

    if args.reference.is_some() {
        println!("arrived ref block");
        let reference_genome = fasta::Reader::from_file(args.reference.unwrap())?;
        let reference_sequences: Vec<fasta::Record> =
            reference_genome.records().map(|r| r.unwrap()).collect();
        let reference_string = reference_sequences
            .into_iter()
            .map(|seq| String::from_utf8(seq.seq().to_owned()).unwrap_or_default())
            .collect::<Vec<String>>()
            .join("");
        let reference = reference_string.as_bytes();
        let reference_length = reference.len();

        let score = |a: u8, b: u8| if a == b { 1 } else { -1 };
        let mut aligner = Aligner::new(-5, -1, score);
        // let mut unique_alignments: HashMap<String, Vec<fasta::Record>> = HashMap::new();
        println!("before mask array");
        let mut mask_array: MaskArray = MaskArray::new(reference_length);
        println!("after mask array");

        for record in &all_records {
            let alignment = aligner.local(&reference, record.seq());
            println!("got here");
            mask_array.apply_mask(alignment.ystart, alignment.yend);
        }

        let re = Regex::new(r"0+").unwrap();
        let mask_array_string = mask_array.to_string();
        let splits = re.split(&mask_array_string);
        for split in splits {
            println!("\"{}\"", split);
        }
    }

    for record in all_records {
        let record_length = record.seq().len();
        seq_lengths.push(record_length);
        *length_inter
            .entry(record_length / args.interval)
            .or_insert(0) += 1;

        gc_count += record
            .seq()
            .iter()
            .filter(|&&c| c == b'G' || c == b'C')
            .count();
        at_count += record
            .seq()
            .iter()
            .filter(|&&c| c == b'A' || c == b'T')
            .count();
    }

    // Calculate N25, N50, N75, N90 and counts
    seq_lengths.sort_unstable_by(|a, b| b.cmp(a));
    let max_seq = &seq_lengths.iter().max().unwrap_or(&0).clone();
    let min_seq = &seq_lengths.iter().min().unwrap_or(&0).clone();

    let total_bp: usize = (&seq_lengths).iter().sum();

    let [mut n25, mut n50, mut n75, mut n90] = [0; 4];
    let [mut l25, mut l50, mut l75, mut l90] = [0; 4];

    let mut cumulative_length = 0;
    let mut seq_count = 0;

    for len in seq_lengths {
        cumulative_length += len;
        seq_count += 1;

        if cumulative_length >= total_bp * 25 / 100 && n25 == 0 {
            n25 = cumulative_length.clone();
            l25 = seq_count.clone();
        }
        if cumulative_length >= total_bp * 50 / 100 && n50 == 0 {
            n50 = cumulative_length.clone();
            l50 = seq_count.clone();
        }
        if cumulative_length >= total_bp * 75 / 100 && n75 == 0 {
            n75 = cumulative_length.clone();
            l75 = seq_count.clone();
        }
        if cumulative_length >= total_bp * 90 / 100 && n90 == 0 {
            n90 = cumulative_length.clone();
            l90 = seq_count.clone();
        }
    }

    let gc_content = gc_count as f64 / total_bp as f64 * 100.0;
    let at_content = at_count as f64 / total_bp as f64 * 100.0;

    outputs.push(format!("Number of sequences: {}", num_seq));
    outputs.push(format!("Total base pairs: {}", total_bp));
    outputs.push(format!("GC content: {:.2}%", gc_content));
    outputs.push(format!("AT content: {:.2}%", at_content));
    outputs.push(format!("Max sequence {}", max_seq));
    outputs.push(format!("Min sequence {}", min_seq));
    outputs.push(format!("N25: {}", n25));
    outputs.push(format!("N50: {}", n50));
    outputs.push(format!("N75: {}", n75));
    outputs.push(format!("N90: {}", n90));
    outputs.push(format!("L25: {}", l25));
    outputs.push(format!("L50: {}", l50));
    outputs.push(format!("L75: {}", l75));
    outputs.push(format!("L90: {}", l90));
    outputs.push(format!("Length intervals ({} bp):", args.interval));

    let mut vec: Vec<_> = length_inter.iter().collect();
    vec.sort_by_key(|&(inter, _count)| inter);
    for &(inter, count) in &vec {
        outputs.push(format!(
            "{}-{}\t{}",
            inter * args.interval,
            (inter + 1) * args.interval - 1,
            count
        ));
    }

    let interval = args.interval;
    let path = args.fasta.to_string_lossy().to_string();
    let header = format!("Interval: {interval}\nPath: {path}\n");
    let out = outputs.join("\n");

    println!("{}", header);
    println!("{}", out);

    Ok(())
}
