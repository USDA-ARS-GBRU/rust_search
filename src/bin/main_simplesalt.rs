use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use aho_corasick::AhoCorasick;
use std::io;

#[derive(Parser, Debug)]
struct Args {
    #[arg(short, long)] file: String,
    #[arg(short, long)] patterns: String,
    /// Max Delta G threshold (kcal/mol)
    #[arg(short, long, default_value_t = -10.0)] threshold: f64,
    /// Monovalent salt Na+ (mM) - Primer3 default 50.0
    #[arg(long, default_value_t = 50.0)] na: f64,
    /// Divalent salt Mg2+ (mM) - Primer3 default 1.5
    #[arg(long, default_value_t = 1.5)] mg: f64,
    /// dNTPs (mM) - Primer3 default 0.6
    #[arg(long, default_value_t = 0.6)] dntp: f64,
    /// Primer concentration (nM) - Primer3 default 200.0
    #[arg(long, default_value_t = 200.0)] dnac: f64,
    /// Temperature (C) for Delta G - default 37.0
    #[arg(long, default_value_t = 37.0)] temp: f64,
}

struct ThermoParams {
    dh: f64,
    ds: f64,
}

// SantaLucia 1998 parameters
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
    let mut total_dh = 0.0;
    let mut total_ds = 0.0;

    // Initiation (SantaLucia 1998)
    total_dh += 0.2;
    total_ds += -5.7;

    // Nearest Neighbor sum
    for i in 0..seq.len() - 1 {
        let p = get_nn_params(seq[i], seq[i+1]);
        total_dh += p.dh;
        total_ds += p.ds;
    }

    // Salt correction (Santalucia 2004 / Primer3 default style)
    // Effect on Delta S: ds_corrected = ds + 0.368 * (N-1) * ln([Na_equivalent])
    let na_eq = args.na + 120.0 * (args.mg - args.dntp).sqrt();
    let salt_corr = 0.368 * (seq.len() as f64 - 1.0) * (na_eq / 1000.0).ln();
    total_ds += salt_corr;

    let t_kelvin = args.temp + 273.15;
    let delta_g = total_dh - (t_kelvin * total_ds / 1000.0);

    // Tm calculation for Heterodimer
    let r = 1.9872; // gas constant cal/(K*mol)
    let c = args.dnac / 1e9;
    let tm = (1000.0 * total_dh) / (total_ds + r * (c / 4.0).ln()) - 273.15;

    (delta_g, tm)
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    let mut pattern_reader = parse_fastx_file(&args.patterns).expect("Invalid pattern file");
    let mut all_motifs = Vec::new();
    let mut all_seeds = Vec::new();

    while let Some(record) = pattern_reader.next() {
        let rec = record.unwrap();
        let seq = rec.seq().to_ascii_uppercase();
        let rc = rec.reverse_complement().to_ascii_uppercase();
        for s in vec![seq, rc] {
            all_seeds.push(s[0..7].to_vec()); // 7-mer seed
            all_motifs.push(s);
        }
    }

    let ac = AhoCorasick::new(&all_seeds).unwrap();
    let mut reader = parse_fastx_file(&args.file).expect("Genome file error");
    let chunk_size = 1_000_000;
    let overlap = 100;

    while let Some(record) = reader.next() {
        let rec = record.unwrap();
        let seq_id = String::from_utf8_lossy(rec.id()).to_string();
        let full_seq = rec.seq();

        (0..full_seq.len()).into_par_iter().step_by(chunk_size - overlap).for_each(|start| {
            let end = (start + chunk_size).min(full_seq.len());
            let chunk = full_seq[start..end].to_ascii_uppercase();
            
            for mat in ac.find_iter(&chunk) {
                let motif = &all_motifs[mat.pattern().as_usize()];
                let hit_pos = mat.start();
                let v_end = (hit_pos + motif.len()).min(chunk.len());
                let vicinity = &chunk[hit_pos..v_end];

                if vicinity.len() == motif.len() {
                    let (dg, tm) = calculate_thermo(vicinity, &args);
                    if dg <= args.threshold {
                        println!("{}\t{}\t{:.2}\t{:.2}\t{}", 
                            seq_id, start + hit_pos, dg, tm, 
                            String::from_utf8_lossy(motif));
                    }
                }
            }
        });
    }
    Ok(())
}
