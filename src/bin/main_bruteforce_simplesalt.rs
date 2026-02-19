use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use std::io;

#[derive(Parser, Debug)]
pub struct Args {
    #[arg(short, long)] 
    file: String,
    #[arg(short, long)] 
    patterns: String,
    /// Max Delta G threshold (kcal/mol)
    #[arg(short, long, default_value_t = -10.0)] 
    threshold: f64,
    /// Monovalent salt Na+ (mM)
    #[arg(long, default_value_t = 50.0)] 
    na: f64,
    /// Divalent salt Mg2+ (mM)
    #[arg(long, default_value_t = 1.5)] 
    mg: f64,
    /// dNTPs (mM)
    #[arg(long, default_value_t = 0.6)] 
    dntp: f64,
    /// Primer concentration (nM)
    #[arg(long, default_value_t = 200.0)] 
    dnac: f64,
    /// Temperature (C)
    #[arg(long, default_value_t = 37.0)] 
    temp: f64,
}

struct ThermoParams {
    dh: f64,
    ds: f64,
}

fn get_nn_params(a: u8, b: u8) -> ThermoParams {
    match (a, b) {
        (b'A', b'A') | (b'T', b'T') => ThermoParams { dh: -7.9, ds: -22.2 },
        (b'A', b'T') => ThermoParams { dh: -7.2, ds: -20.4 },
        (b'T', b'A') => ThermoParams { dh: -7.2, ds: -21.3 },
        (b'C', b'A') | (b'T', b'G') => ThermoParams { dh: -8.5, ds: -22.7 },
        (b'G', b'T') | (b'A', b'C') => ThermoParams { dh: -8.4, ds: -22.4 },
        (b'C', b'T') | (b'A', b'G') => ThermoParams { dh: -7.8, ds: -21.0 },
        (b'G', b'A') | (b'T', b'C') => ThermoParams { dh: -8.2, ds: -22.2 },
        (b'C', b'G') => ThermoParams { dh: -10.6, ds: -27.2 },
        (b'G', b'C') => ThermoParams { dh: -9.8, ds: -24.4 },
        (b'C', b'C') | (b'G', b'G') => ThermoParams { dh: -8.0, ds: -19.9 },
        _ => ThermoParams { dh: 0.0, ds: 0.0 },
    }
}

fn calculate_thermo(seq: &[u8], args: &Args) -> (f64, f64) {
    let mut total_dh = 0.2;
    let mut total_ds = -5.7;

    for i in 0..seq.len() - 1 {
        let p = get_nn_params(seq[i], seq[i+1]);
        total_dh += p.dh;
        total_ds += p.ds;
    }

    let na_eq = args.na + 120.0 * (args.mg - args.dntp).max(0.0).sqrt();
    let salt_corr = 0.368 * (seq.len() as f64 - 1.0) * (na_eq / 1000.0).ln();
    total_ds += salt_corr;

    let t_kelvin = args.temp + 273.15;
    let delta_g = total_dh - (t_kelvin * total_ds / 1000.0);
    
    let r = 1.9872;
    let c = args.dnac / 1e9;
    let tm = (1000.0 * total_dh) / (total_ds + r * (c / 4.0).ln()) - 273.15;

    (delta_g, tm)
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    // 1. Load patterns
    let mut pattern_reader = parse_fastx_file(&args.patterns).expect("Invalid patterns file");
    let mut patterns: Vec<Vec<u8>> = Vec::new();
    while let Some(record) = pattern_reader.next() {
        let rec = record.expect("Error reading pattern record");
        patterns.push(rec.seq().to_ascii_uppercase());
    }

    // 2. Process Genome
    let mut reader = parse_fastx_file(&args.file).expect("Invalid genome file");

    while let Some(record) = reader.next() {
        let rec = record.expect("FASTA error");
        let seq_id = String::from_utf8_lossy(rec.id()).to_string();
        let genome_fwd = rec.seq().to_ascii_uppercase();
        
        let mut genome_rev = rec.reverse_complement();
        genome_rev.make_ascii_uppercase();

        let genome_len = genome_fwd.len();

        // 3. Sliding window parallelized
        (0..genome_len).into_par_iter().for_each(|i| {
            for pattern in &patterns {
                let p_len = pattern.len();
                if i + p_len <= genome_len {
                    
                    // Check FWD strand window at this position
                    let window_fwd = &genome_fwd[i..i+p_len];
                    let (dg_fwd, tm_fwd) = calculate_thermo(window_fwd, &args);
                    
                    if dg_fwd <= args.threshold {
                        println!("{}\t{}\tFWD\t{:.2}\t{:.2}\t{}", 
                            seq_id, i, dg_fwd, tm_fwd, String::from_utf8_lossy(window_fwd));
                    }

                    // Check REV strand window at this position
                    let window_rev = &genome_rev[i..i+p_len];
                    let (dg_rev, tm_rev) = calculate_thermo(window_rev, &args);
                    
                    if dg_rev <= args.threshold {
                        // Position i on genome_rev is (genome_len - i - p_len) on forward strand
                        println!("{}\t{}\tREV\t{:.2}\t{:.2}\t{}", 
                            seq_id, i, dg_rev, tm_rev, String::from_utf8_lossy(window_rev));
                    }
                }
            }
        });
    }

    Ok(())
}
