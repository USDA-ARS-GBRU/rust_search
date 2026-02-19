use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_double};

pub mod thal {
    use super::*;

    #[repr(C)]
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum ThalAlignmentType {
        Any = 1,
        End1 = 2,
        End2 = 3,
        Hairpin = 4,
    }

    #[repr(C)]
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub enum ThalMode {
        Fast = 0,
        General = 1,
        DebugFast = 2,
        Debug = 3,
        Struct = 4,
    }

    #[repr(C)]
    #[derive(Debug, Clone)]
    pub struct ThalArgsC {
        pub alignment_type: ThalAlignmentType,
        pub max_loop: c_int,
        pub mv: c_double,
        pub dv: c_double,
        pub dntp: c_double,
        pub dna_conc: c_double,
        pub temp: c_double,
        pub dimer: c_int,
    }

    #[repr(C)]
    pub struct ThalResultsC {
        pub msg: [c_char; 255],
        pub temp: c_double,
        pub dg: c_double,
        pub ds: c_double,
        pub dh: c_double,
        pub align_end_1: c_int,
        pub align_end_2: c_int,
        pub sec_struct: *mut c_char,
    }

    #[repr(C)]
    pub struct ThalParametersC {
        pub dangle_dh: *mut c_char,
        pub dangle_ds: *mut c_char,
        pub loops_dh: *mut c_char,
        pub loops_ds: *mut c_char,
        pub stack_dh: *mut c_char,
        pub stack_ds: *mut c_char,
        pub stackmm_dh: *mut c_char,
        pub stackmm_ds: *mut c_char,
        pub tetraloop_dh: *mut c_char,
        pub tetraloop_ds: *mut c_char,
        pub triloop_dh: *mut c_char,
        pub triloop_ds: *mut c_char,
        pub tstack_tm_inf_ds: *mut c_char,
        pub tstack_dh: *mut c_char,
        pub tstack2_dh: *mut c_char,
        pub tstack2_ds: *mut c_char,
    }

    extern "C" {
        #[link_name = "thal"]
        pub fn thal_ffi(
            oligo1: *const u8,
            oligo2: *const u8,
            a: *const ThalArgsC,
            mode: ThalMode,
            o: *mut ThalResultsC,
        );
        pub fn thal_set_null_parameters(a: *mut ThalParametersC) -> c_int;
        pub fn thal_load_parameters(
            path: *const c_char,
            a: *mut ThalParametersC,
            o: *mut ThalResultsC,
        ) -> c_int;
        pub fn get_thermodynamic_values(tp: *const ThalParametersC, o: *mut ThalResultsC) -> c_int;
        pub fn thal_free_parameters(a: *mut ThalParametersC) -> c_int;
        pub fn destroy_thal_structures();
    }

    #[derive(Debug, Clone)]
    pub struct ThalArgs {
        pub alignment_type: ThalAlignmentType,
        pub max_loop: i32,
        pub mv: f64,
        pub dv: f64,
        pub dntp: f64,
        pub dna_conc: f64,
        pub temp: f64,
        pub dimer: i32,
    }

    #[derive(Debug, Clone)]
    pub struct ThalResults {
        pub msg: String,
        pub temp: f64,
        pub dg: f64,
        pub ds: f64,
        pub dh: f64,
        pub align_end_1: i32,
        pub align_end_2: i32,
        pub sec_struct: Option<String>,
    }

    pub const ABSOLUTE_ZERO: f64 = 273.15;
    pub const THAL_ERROR_SCORE: f64 = f64::NEG_INFINITY;

    use std::sync::Mutex;
    static INITIALIZED: Mutex<Option<Result<(), String>>> = Mutex::new(None);

    pub fn ensure_parameters_loaded(path: &str) -> Result<(), String> {
        let mut init = INITIALIZED.lock().unwrap();
        if let Some(ref res) = *init {
            return res.clone();
        }

        let res = (|| unsafe {
            let mut params: ThalParametersC = std::mem::zeroed();
            thal_set_null_parameters(&mut params);
            
            let mut results: ThalResultsC = std::mem::zeroed();
            results.sec_struct = std::ptr::null_mut();

            let c_path = CString::new(path).map_err(|e| e.to_string())?;
            if thal_load_parameters(c_path.as_ptr(), &mut params, &mut results) != 0 {
                return Err(format!("Failed to load parameters from {}: {}", path, 
                    CStr::from_ptr(results.msg.as_ptr()).to_string_lossy()));
            }

            if get_thermodynamic_values(&params, &mut results) != 0 {
                return Err(format!("Failed to initialize thermodynamic values: {}", 
                    CStr::from_ptr(results.msg.as_ptr()).to_string_lossy()));
            }
            
            Ok(())
        })();

        *init = Some(res.clone());
        res
    }

    pub fn thal_wrapper(
        seq1: &[u8],
        seq2: &[u8],
        args: &ThalArgs,
        mode: ThalMode,
    ) -> ThalResults {
        let c_seq1 = CString::new(seq1).unwrap();
        let c_seq2 = CString::new(seq2).unwrap();

        let c_args = ThalArgsC {
            alignment_type: args.alignment_type,
            max_loop: args.max_loop,
            mv: args.mv,
            dv: args.dv,
            dntp: args.dntp,
            dna_conc: args.dna_conc,
            temp: args.temp,
            dimer: args.dimer,
        };

        unsafe {
            let mut c_results: ThalResultsC = std::mem::zeroed();
            c_results.sec_struct = std::ptr::null_mut();
            
            thal_ffi(
                c_seq1.as_ptr() as *const u8,
                c_seq2.as_ptr() as *const u8,
                &c_args,
                mode,
                &mut c_results,
            );

            let msg = CStr::from_ptr(c_results.msg.as_ptr()).to_string_lossy().into_owned();
            let sec_struct = if !c_results.sec_struct.is_null() {
                let s = CStr::from_ptr(c_results.sec_struct).to_string_lossy().into_owned();
                libc::free(c_results.sec_struct as *mut libc::c_void);
                Some(s)
            } else {
                None
            };

            ThalResults {
                msg,
                temp: c_results.temp,
                dg: c_results.dg,
                ds: c_results.ds,
                dh: c_results.dh,
                align_end_1: c_results.align_end_1,
                align_end_2: c_results.align_end_2,
                sec_struct,
            }
        }
    }

    // Re-expose the thal function as the wrapper
    pub use thal_wrapper as thal;

    pub fn create_default_args() -> ThalArgs {
        ThalArgs {
            alignment_type: ThalAlignmentType::Any,
            max_loop: 30,
            mv: 50.0,
            dv: 0.0,
            dntp: 0.8,
            dna_conc: 50.0,
            temp: 37.0 + ABSOLUTE_ZERO,
            dimer: 1,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_thal_ffi() {
            ensure_parameters_loaded("primer3/src/primer3_config/").expect("Failed to load params");
            let args = create_default_args();
            // Test a simple perfect match
            let seq1 = b"ATGCGATCGATCG";
            let seq2 = b"ATGCGATCGATCG";
            let result = thal(seq1, seq2, &args, ThalMode::Fast);
            
            assert!(result.temp > 0.0);
            assert!(result.dg < 0.0);
            assert_eq!(result.msg, "");
        }

        #[test]
        fn test_thal_mismatch() {
            ensure_parameters_loaded("primer3/src/primer3_config/").expect("Failed to load params");
            let args = create_default_args();
            // Perfect match
            let seq1 = b"ATGCGATCGATCG";
            let seq2 = b"ATGCGATCGATCG";
            let result1 = thal(seq1, seq2, &args, ThalMode::Fast);

            // One mismatch
            let seq3 = b"ATGCGATCGATCT"; 
            let result2 = thal(seq1, seq3, &args, ThalMode::Fast);

            assert!(result2.temp < result1.temp);
            assert!(result2.dg > result1.dg);
        }
    }
}

pub use thal::*;
