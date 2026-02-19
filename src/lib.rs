/// Rust implementation of primer3 thermodynamic calculations
/// 
/// This library provides Rust bindings for primer3's thermodynamic alignment
/// calculations, based on the C source code from thal_main.c, thal.h, and thal.c.

pub mod thal {
    use std::f64;

    /// Alignment type enumeration
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum ThalAlignmentType {
        /// Any alignment (default for dimer)
        Any = 1,
        /// 3' end alignment (end1)
        End1 = 2,
        /// 5' end alignment (end2)
        End2 = 3,
        /// Hairpin structure
        Hairpin = 4,
    }

    /// Calculation mode enumeration
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum ThalMode {
        /// Fast mode - score only
        Fast = 0,
        /// General mode - without debug
        General = 1,
        /// Debug mode with fast
        DebugFast = 2,
        /// Debug mode
        Debug = 3,
        /// Calculate secondary structures as string
        Struct = 4,
    }

    /// Arguments for thermodynamic alignment calculation
    #[derive(Debug, Clone)]
    pub struct ThalArgs {
        /// Type of alignment (Any, End1, End2, Hairpin)
        pub alignment_type: ThalAlignmentType,
        /// Maximum size of loop to consider (0-30 bp)
        pub max_loop: i32,
        /// Concentration of monovalent cations (mM)
        pub mv: f64,
        /// Concentration of divalent cations (mM)
        pub dv: f64,
        /// Concentration of dNTPs (mM)
        pub dntp: f64,
        /// Concentration of DNA oligonucleotides (nM)
        pub dna_conc: f64,
        /// Temperature in Kelvin (for calculations)
        pub temp: f64,
        /// If non-zero, calculate dimer structure; otherwise hairpin
        pub dimer: i32,
    }

    /// Results from thermodynamic alignment calculation
    #[derive(Debug, Clone)]
    pub struct ThalResults {
        /// Error message (if any)
        pub msg: String,
        /// Melting temperature (°C)
        pub temp: f64,
        /// Gibbs free energy (cal/mol)
        pub dg: f64,
        /// Entropy (cal/mol/K)
        pub ds: f64,
        /// Enthalpy (cal/mol)
        pub dh: f64,
        /// Alignment end position in first sequence
        pub align_end_1: i32,
        /// Alignment end position in second sequence
        pub align_end_2: i32,
        /// Secondary structure representation (if requested)
        pub sec_struct: Option<String>,
    }

    /// Constants from primer3
    pub const THAL_MAX_ALIGN: usize = 60;
    pub const THAL_MAX_SEQ: usize = 10000;
    pub const MAX_LOOP: i32 = 30;
    pub const MIN_LOOP: i32 = 0;
    pub const ABSOLUTE_ZERO: f64 = 273.15;
    pub const THAL_ERROR_SCORE: f64 = f64::NEG_INFINITY;

    /// SantaLucia 1998 nearest neighbor parameters
    /// ΔH in kcal/mol, ΔS in cal/mol/K
    #[derive(Debug, Clone, Copy)]
    pub struct NNParams {
        pub dh: f64,
        pub ds: f64,
    }

    /// Get nearest neighbor parameters for a base pair
    pub fn get_nn_params(base1: u8, base2: u8) -> NNParams {
        match (base1, base2) {
            (b'A', b'A') | (b'T', b'T') => NNParams { dh: -7.9, ds: -22.2 },
            (b'A', b'T') => NNParams { dh: -7.2, ds: -20.4 },
            (b'T', b'A') => NNParams { dh: -7.2, ds: -21.3 },
            (b'C', b'A') | (b'T', b'G') => NNParams { dh: -8.5, ds: -22.7 },
            (b'G', b'T') | (b'A', b'C') => NNParams { dh: -8.4, ds: -22.4 },
            (b'C', b'T') | (b'A', b'G') => NNParams { dh: -7.8, ds: -21.0 },
            (b'G', b'A') | (b'T', b'C') => NNParams { dh: -8.2, ds: -22.2 },
            (b'C', b'G') => NNParams { dh: -10.6, ds: -27.2 },
            (b'G', b'C') => NNParams { dh: -9.8, ds: -24.4 },
            (b'C', b'C') | (b'G', b'G') => NNParams { dh: -8.0, ds: -19.9 },
            _ => NNParams { dh: 0.0, ds: 0.0 },
        }
    }

    /// Get initiation parameters based on terminal base pairs
    pub fn get_initiation_params(first_base: u8, last_base: u8) -> NNParams {
        match (first_base, last_base) {
            // A-T terminal pairs
            (b'A', b'T') | (b'T', b'A') => NNParams { dh: 2.3, ds: 4.1 },
            // G-C terminal pairs
            (b'G', b'C') | (b'C', b'G') => NNParams { dh: 0.1, ds: -2.8 },
            // Mixed terminal pairs
            (b'A', b'G') | (b'G', b'A') | (b'T', b'C') | (b'C', b'T') => NNParams { dh: 1.2, ds: 0.7 },
            (b'A', b'C') | (b'C', b'A') | (b'T', b'G') | (b'G', b'T') => NNParams { dh: 1.2, ds: 0.7 },
            _ => NNParams { dh: 0.2, ds: -5.7 },
        }
    }

    /// Calculate effective sodium concentration using SantaLucia 2004 model
    pub fn calculate_na_equivalent(na: f64, mg: f64, dntp: f64) -> f64 {
        let mg_eff = if mg > dntp { mg - dntp } else { 0.0 };
        na + 120.0 * mg_eff.sqrt()
    }

    /// Check if sequence is self-complementary (symmetric)
    pub fn is_self_complementary(seq: &[u8]) -> bool {
        let n = seq.len();
        for i in 0..n / 2 {
            let complement = match seq[n - 1 - i] {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                _ => return false,
            };
            if seq[i] != complement {
                return false;
            }
        }
        true
    }

    /// Create default ThalArgs
    pub fn create_default_args() -> ThalArgs {
        ThalArgs {
            alignment_type: ThalAlignmentType::Any,
            max_loop: MAX_LOOP,
            mv: 50.0,
            dv: 1.5,
            dntp: 0.6,
            dna_conc: 50.0,
            temp: 37.0 + ABSOLUTE_ZERO,
            dimer: 1,
        }
    }

    /// Calculate thermodynamic parameters for a sequence
    pub fn calculate_thermo(seq: &[u8], args: &ThalArgs) -> ThalResults {
        if seq.len() < 2 {
            return ThalResults {
                msg: "Sequence too short (minimum 2 bp)".to_string(),
                temp: THAL_ERROR_SCORE,
                dg: 0.0,
                ds: 0.0,
                dh: 0.0,
                align_end_1: 0,
                align_end_2: 0,
                sec_struct: None,
            };
        }

        let mut result = ThalResults {
            msg: String::new(),
            temp: 0.0,
            dg: 0.0,
            ds: 0.0,
            dh: 0.0,
            align_end_1: 0,
            align_end_2: 0,
            sec_struct: None,
        };

        let mut total_dh = 0.0;
        let mut total_ds = 0.0;

        // Initiation parameters based on terminal base pairs
        let init_params = get_initiation_params(seq[0], seq[seq.len() - 1]);
        total_dh += init_params.dh;
        total_ds += init_params.ds;

        // Nearest neighbor sum
        for i in 0..seq.len() - 1 {
            let p = get_nn_params(seq[i], seq[i + 1]);
            total_dh += p.dh;
            total_ds += p.ds;
        }

        // Salt correction (SantaLucia 2004)
        let na_eq = calculate_na_equivalent(args.mv, args.dv, args.dntp);
        let salt_corr = 0.368 * (seq.len() as f64 - 1.0) * (na_eq / 1000.0).ln();
        total_ds += salt_corr;

        // Calculate ΔG at specified temperature (in kcal/mol)
        let delta_g_kcal = total_dh - (args.temp * total_ds / 1000.0);
        
        // Store ΔH, ΔS, and ΔG (convert to cal/mol to match primer3-py output)
        result.dh = total_dh * 1000.0; // Convert kcal/mol to cal/mol
        result.ds = total_ds;
        result.dg = delta_g_kcal * 1000.0; // Convert kcal/mol to cal/mol

        // Calculate Tm
        let r = 1.9872; // gas constant cal/(K*mol)
        let c = args.dna_conc / 1e9; // convert nM to M

        let is_symmetric = is_self_complementary(seq);
        let c_factor = if is_symmetric { 2.0 } else { 4.0 };

        let c_term = if c > 0.0 { (c / c_factor).ln() } else { 0.0 };
        
        if (total_ds + r * c_term).abs() > 1e-10 {
            result.temp = (1000.0 * total_dh) / (total_ds + r * c_term) - ABSOLUTE_ZERO;
        } else {
            result.msg = "Invalid thermodynamic parameters".to_string();
            result.temp = THAL_ERROR_SCORE;
        }

        result.align_end_1 = seq.len() as i32;
        result.align_end_2 = seq.len() as i32;

        result
    }

    /// Perform thermodynamic alignment calculation
    pub fn thal(
        seq1: &[u8],
        seq2: &[u8],
        args: &ThalArgs,
        _mode: ThalMode,
    ) -> ThalResults {
        if seq1.len() > THAL_MAX_ALIGN || seq2.len() > THAL_MAX_ALIGN {
            return ThalResults {
                msg: format!(
                    "Sequence too long (max {} bp for one sequence)",
                    THAL_MAX_ALIGN
                ),
                temp: THAL_ERROR_SCORE,
                dg: 0.0,
                ds: 0.0,
                dh: 0.0,
                align_end_1: 0,
                align_end_2: 0,
                sec_struct: None,
            };
        }

        if seq1.len() > THAL_MAX_SEQ || seq2.len() > THAL_MAX_SEQ {
            return ThalResults {
                msg: format!("Sequence too long (max {} bp)", THAL_MAX_SEQ),
                temp: THAL_ERROR_SCORE,
                dg: 0.0,
                ds: 0.0,
                dh: 0.0,
                align_end_1: 0,
                align_end_2: 0,
                sec_struct: None,
            };
        }

        // For hairpin calculation, use the same sequence
        if args.dimer == 0 {
            calculate_thermo(seq1, args)
        } else {
            // For dimer, we need to find the best alignment
            // This is a simplified version that assumes perfect match
            let min_len = seq1.len().min(seq2.len());
            let mut best_result = calculate_thermo(&seq1[0..min_len], args);
            best_result.align_end_1 = min_len as i32;
            best_result.align_end_2 = min_len as i32;
            best_result
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_nn_params() {
            let params = get_nn_params(b'A', b'A');
            assert_eq!(params.dh, -7.9);
            assert_eq!(params.ds, -22.2);
        }

        #[test]
        fn test_initiation_params() {
            let params = get_initiation_params(b'A', b'T');
            assert_eq!(params.dh, 2.3);
            assert_eq!(params.ds, 4.1);
        }

        #[test]
        fn test_self_complementary() {
            assert!(is_self_complementary(b"ATCGAT"));
            assert!(!is_self_complementary(b"ATCGAA"));
        }

        #[test]
        fn test_calculate_thermo() {
            let args = create_default_args();
            let result = calculate_thermo(b"ATGCGATCGATCG", &args);
            assert!(result.temp > 0.0);
            assert!(result.dg < 0.0);
        }
    }
}

// Re-export commonly used types
pub use thal::{
    ThalAlignmentType, ThalMode, ThalArgs, ThalResults, NNParams,
    THAL_MAX_ALIGN, THAL_MAX_SEQ, MAX_LOOP, MIN_LOOP, ABSOLUTE_ZERO, THAL_ERROR_SCORE,
};
