use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use aho_corasick::AhoCorasick;
use std::io;
use std::collections::{HashMap, HashSet};
use rust_search::{
    thal, ThalArgs, ThalAlignmentType, ThalMode, ABSOLUTE_ZERO,
};

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
    /// Primer concentration (nM) - Primer3 default 50.0
    #[arg(long, default_value_t = 50.0)] dnac: f64,
    /// Temperature (C) for Delta G - default 37.0
    #[arg(long, default_value_t = 37.0)] temp: f64,
    /// Maximum loop size (bp) - Primer3 default 30
    #[arg(long, default_value_t = 30)] max_loop: i32,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    // Initialize thermodynamic parameters from primer3_config
    rust_search::thal::ensure_parameters_loaded("primer3/src/primer3_config/")
        .expect("Failed to load thermodynamic parameters");

    let mut pattern_reader = parse_fastx_file(&args.patterns).expect("Invalid pattern file");
    let mut all_motifs = Vec::new();

    // Mapping from unique seed to list of (motif_idx, offset)
    let mut seed_map: HashMap<Vec<u8>, Vec<(usize, usize)>> = HashMap::new();

    while let Some(record) = pattern_reader.next() {
        let rec = record.unwrap();
        let seq = rec.seq().to_ascii_uppercase();
        let rc = rec.reverse_complement().to_ascii_uppercase();
        for s in vec![seq, rc] {
            let motif_idx = all_motifs.len();
            all_motifs.push(s.clone());

            if s.len() >= 7 {
                for offset in 0..=(s.len() - 7) {
                    let seed = s[offset..offset+7].to_vec();
                    seed_map.entry(seed).or_default().push((motif_idx, offset));
                }
            }
        }
    }

    let mut unique_seeds = Vec::new();
    let mut seed_to_motifs = Vec::new();
    for (seed, motifs) in seed_map {
        unique_seeds.push(seed);
        seed_to_motifs.push(motifs);
    }

    let ac = AhoCorasick::new(&unique_seeds).unwrap();
    let mut reader = parse_fastx_file(&args.file).expect("Genome file error");
    let chunk_size = 1_000_000;
    let overlap = 100;

    // Create thal_args for thermodynamic calculations
    let thal_args = ThalArgs {
        alignment_type: ThalAlignmentType::Any,
        max_loop: args.max_loop,
        mv: args.na,
        dv: args.mg,
        dntp: args.dntp,
        dna_conc: args.dnac,
        temp: args.temp + ABSOLUTE_ZERO,
        dimer: 1,
    };

    while let Some(record) = reader.next() {
        let rec = record.unwrap();
        let seq_id = String::from_utf8_lossy(rec.id()).to_string();
        let full_seq = rec.seq();

        (0..full_seq.len()).into_par_iter().step_by(chunk_size - overlap).for_each(|start| {
            let end = (start + chunk_size).min(full_seq.len());
            let chunk = full_seq[start..end].to_ascii_uppercase();
            
            let mut evaluated = HashSet::new();

            for mat in ac.find_overlapping_iter(&chunk) {
                let seed_idx = mat.pattern().as_usize();
                let motifs = &seed_to_motifs[seed_idx];
                let hit_pos = mat.start();

                for &(motif_idx, offset) in motifs {
                    let motif = &all_motifs[motif_idx];
                    let genome_start = hit_pos as isize - offset as isize;
                    let genome_end = genome_start + motif.len() as isize;

                    if genome_start >= 0 && genome_end <= chunk.len() as isize {
                        let is_last_chunk = end == full_seq.len();
                        if is_last_chunk || genome_start < (chunk_size - overlap) as isize {
                            if evaluated.insert((motif_idx, genome_start)) {
                                let vicinity = &chunk[genome_start as usize .. genome_end as usize];
                                // Use the thal function from the library
                                let result = thal::thal(motif, vicinity, &thal_args, ThalMode::Fast);

                                // ΔG is in cal/mol, convert to kcal/mol for threshold comparison
                                let dg_kcal = result.dg / 1000.0;

                                if dg_kcal <= args.threshold {
                                    println!("{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{}",
                                        seq_id, start + genome_start as usize, dg_kcal, result.temp,
                                        result.dh / 1000.0, result.ds,
                                        String::from_utf8_lossy(motif));
                                }
                            }
                        }
                    }
                }
            }
        });
    }
    Ok(())
}
