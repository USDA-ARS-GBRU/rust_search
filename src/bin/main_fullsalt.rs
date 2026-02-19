use clap::Parser;
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use aho_corasick::AhoCorasick;
use std::io;
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
            
            for mat in ac.find_iter(&chunk) {
                let motif = &all_motifs[mat.pattern().as_usize()];
                let hit_pos = mat.start();
                let v_end = (hit_pos + motif.len()).min(chunk.len());
                let vicinity = &chunk[hit_pos..v_end];

                if vicinity.len() == motif.len() {
                    // Use the thal function from the library
                    let result = thal::thal(motif, vicinity, &thal_args, ThalMode::Fast);
                    
                    // ΔG is in cal/mol, convert to kcal/mol for threshold comparison
                    let dg_kcal = result.dg / 1000.0;
                    
                    if dg_kcal <= args.threshold {
                        println!("{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{}",
                            seq_id, start + hit_pos, dg_kcal, result.temp, 
                            result.dh / 1000.0, result.ds,
                            String::from_utf8_lossy(motif));
                    }
                }
            }
        });
    }
    Ok(())
}
