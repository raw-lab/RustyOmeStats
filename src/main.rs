use clap::Parser;
use std::{
    collections::HashMap,
    error::Error,
    fs,
    io::{BufRead, BufReader},
    path::PathBuf,
};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Interval size in number of residues
    #[arg(short, long)]
    interval: usize,

    /// Fasta file or folder
    #[arg(short, long)]
    fasta: PathBuf,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let files = if args.fasta.is_file() {
        vec![args.fasta.clone()]
    } else if args.fasta.is_dir() {
        fs::read_dir(&args.fasta)?
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

    let mut num_seq = 0;
    let mut gc_count = 0;
    let mut at_count = 0;
    let mut length_inter = HashMap::new();
    let mut seq_lengths = vec![];

    for file in files {
        let file = fs::File::open(file)?;
        let reader = BufReader::new(file);

        let mut id: Option<String> = None;
        let mut seq_len = 0;

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('>') {
                // Next sequence in file
                if let Some(_id) = id.take() {
                    seq_lengths.push(seq_len);
                    let inter = seq_len / args.interval;
                    *length_inter.entry(inter).or_insert(0) += 1;
                }
                // Sequence Basic Info
                num_seq += 1;
                id = Some(line[1..].to_string());
                seq_len = 0;
            } else if id.is_some()
                && line
                    .chars()
                    .next()
                    .map_or(false, |c| c.is_ascii_alphanumeric())
            {
                // Sequence data
                seq_len += line.len();
                gc_count += line.chars().filter(|&c| c == 'G' || c == 'C').count();
                at_count += line.chars().filter(|&c| c == 'A' || c == 'T').count();
            }
        }

        // Incorporate totals from last sequence in file
        if let Some(_id) = id {
            seq_lengths.push(seq_len);
            let inter = seq_len / args.interval;
            *length_inter.entry(inter).or_insert(0) += 1;
        }
    }

    // Calculate N25, N50, N75, N90 and counts
    seq_lengths.sort_unstable_by(|a, b| b.cmp(a));
    let max_seq = &seq_lengths.iter().max().unwrap_or(&0).clone();
    let min_seq = &seq_lengths.iter().min().unwrap_or(&0).clone();

    let total_bp: usize = (&seq_lengths).iter().sum();

    let [mut n25, mut n50, mut n75, mut n90] = [0; 4];

    let mut cumulative_length = 0;

    for len in seq_lengths {
        cumulative_length += len;

        if cumulative_length >= total_bp * 25 / 100 && n25 == 0 {
            n25 = cumulative_length.clone();
        }
        if cumulative_length >= total_bp * 50 / 100 && n50 == 0 {
            n50 = cumulative_length.clone();
        }
        if cumulative_length >= total_bp * 75 / 100 && n75 == 0 {
            n75 = cumulative_length.clone();
        }
        if cumulative_length >= total_bp * 90 / 100 && n90 == 0 {
            n90 = cumulative_length.clone();
        }
    }

    let gc_content = gc_count as f64 / total_bp as f64 * 100.0;
    let at_content = at_count as f64 / total_bp as f64 * 100.0;

    let mut outputs: Vec<String> = vec![];
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
