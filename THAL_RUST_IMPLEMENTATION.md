# Rust Implementation of Primer3's thal() Function

## Overview

This document describes the Rust implementation of the primer3 thermodynamic alignment (`thal`) functionality, based on the C source code from `thal_main.c`, `thal.h`, and `thal.c`.

## Architecture

### File Structure

```
src/
├── lib.rs                    # Core thal module (thermodynamic calculations)
└── bin/
    └── main_fullsalt.rs     # Binary that uses the thal module
```

### Module Organization

The `src/lib.rs` file contains the `thal` module with:

1. **Enumerations**
   - `ThalAlignmentType`: Alignment type (Any, End1, End2, Hairpin)
   - `ThalMode`: Calculation mode (Fast, General, Debug, Struct)

2. **Structures**
   - `ThalArgs`: Input parameters for calculations
   - `ThalResults`: Output results from calculations
   - `NNParams`: Nearest neighbor parameters

3. **Functions**
   - `get_nn_params()`: SantaLucia 1998 nearest neighbor parameters
   - `get_initiation_params()`: Terminal base pair-dependent initiation
   - `calculate_na_equivalent()`: Salt correction calculation
   - `is_self_complementary()`: Symmetry detection
   - `calculate_thermo()`: Core thermodynamic calculation
   - `thal()`: Main alignment function

## Data Structures

### ThalArgs
```rust
pub struct ThalArgs {
    pub alignment_type: ThalAlignmentType,  // Type of alignment
    pub max_loop: i32,                      // Maximum loop size (0-30 bp)
    pub mv: f64,                            // Monovalent cation conc (mM)
    pub dv: f64,                            // Divalent cation conc (mM)
    pub dntp: f64,                          // dNTP concentration (mM)
    pub dna_conc: f64,                      // DNA concentration (nM)
    pub temp: f64,                          // Temperature (Kelvin)
    pub dimer: i32,                         // 1 for dimer, 0 for hairpin
}
```

### ThalResults
```rust
pub struct ThalResults {
    pub msg: String,                        // Error message
    pub temp: f64,                          // Melting temperature (°C)
    pub dg: f64,                            // ΔG (cal/mol)
    pub ds: f64,                            // ΔS (cal/mol/K)
    pub dh: f64,                            // ΔH (cal/mol)
    pub align_end_1: i32,                   // Alignment end in seq1
    pub align_end_2: i32,                   // Alignment end in seq2
    pub sec_struct: Option<String>,         // Secondary structure
}
```

## Thermodynamic Model

### 1. Nearest Neighbor Parameters (SantaLucia 1998)

The `get_nn_params()` function returns ΔH (kcal/mol) and ΔS (cal/mol/K) for each base pair:

```
AA/TT: ΔH = -7.9,  ΔS = -22.2
AT:    ΔH = -7.2,  ΔS = -20.4
TA:    ΔH = -7.2,  ΔS = -21.3
CA/TG: ΔH = -8.5,  ΔS = -22.7
GT/AC: ΔH = -8.4,  ΔS = -22.4
CT/AG: ΔH = -7.8,  ΔS = -21.0
GA/TC: ΔH = -8.2,  ΔS = -22.2
CG:    ΔH = -10.6, ΔS = -27.2
GC:    ΔH = -9.8,  ΔS = -24.4
CC/GG: ΔH = -8.0,  ΔS = -19.9
```

### 2. Initiation Parameters (Terminal Base Pair Dependent)

The `get_initiation_params()` function returns different penalties based on terminal bases:

```
A-T pairs: ΔH = 2.3,  ΔS = 4.1
G-C pairs: ΔH = 0.1,  ΔS = -2.8
Mixed:     ΔH = 1.2,  ΔS = 0.7
Default:   ΔH = 0.2,  ΔS = -5.7
```

### 3. Salt Correction (SantaLucia 2004)

```
[Na+]_eq = [Na+] + 120 * sqrt([Mg2+]_eff)
where [Mg2+]_eff = max(0, [Mg2+] - [dNTP])

ΔS_corrected = ΔS + 0.368 * (N-1) * ln([Na+]_eq / 1000)
```

### 4. Thermodynamic Calculations

```
ΔG = ΔH - T * ΔS / 1000

Tm = ΔH / (ΔS + R * ln(C/factor)) - 273.15
where:
  R = 1.9872 cal/(K*mol)
  C = DNA concentration (M)
  factor = 2 for self-complementary, 4 for heterodimer
```

## Comparison with C Implementation

### What IS Implemented

✓ SantaLucia 1998 nearest neighbor parameters
✓ Terminal base pair-dependent initiation
✓ SantaLucia 2004 salt correction
✓ Symmetry-aware Tm calculation
✓ Proper unit conversions
✓ Error handling

### What is NOT Implemented (Requires Full C Library)

✗ Dynamic programming alignment search
✗ Dangling end effects (database-driven)
✗ Mismatch penalties
✗ Loop entropy corrections
✗ Bulge penalties
✗ Internal loop penalties
✗ Optimal alignment discovery

## Usage Example

### Using the Library

```rust
use rust_search::thal::{self, ThalArgs, ThalAlignmentType, ThalMode, ABSOLUTE_ZERO};

let args = ThalArgs {
    alignment_type: ThalAlignmentType::Any,
    max_loop: 30,
    mv: 50.0,
    dv: 1.5,
    dntp: 0.6,
    dna_conc: 50.0,
    temp: 37.0 + ABSOLUTE_ZERO,
    dimer: 1,
};

let result = thal::thal(b"ATGCGATCGATCG", b"ATGCGATCGATCG", &args, ThalMode::Fast);

println!("Tm: {:.2}°C", result.temp);
println!("ΔG: {:.2} cal/mol", result.dg);
println!("ΔH: {:.2} cal/mol", result.dh);
println!("ΔS: {:.2} cal/mol/K", result.ds);
```

### Using the Binary

```bash
# Build
cargo build --release --bin main_fullsalt

# Run with default parameters
cargo run --release --bin main_fullsalt -- \
  --file genome.fasta --patterns primers.fasta

# Run with custom parameters
cargo run --release --bin main_fullsalt -- \
  --file genome.fasta --patterns primers.fasta \
  --na 50.0 --mg 1.5 --dntp 0.6 --dnac 50.0 --temp 37.0 \
  --max_loop 30 --threshold -10.0
```

## Constants

```rust
pub const THAL_MAX_ALIGN: usize = 60;      // Max length for one sequence
pub const THAL_MAX_SEQ: usize = 10000;     // Max length for other sequence
pub const MAX_LOOP: i32 = 30;              // Maximum loop size
pub const MIN_LOOP: i32 = 0;               // Minimum loop size
pub const ABSOLUTE_ZERO: f64 = 273.15;    // Kelvin conversion
pub const THAL_ERROR_SCORE: f64 = f64::NEG_INFINITY;  // Error indicator
```

## Default Parameters

Matching primer3-py defaults:

```
mv (Na+):       50.0 mM
dv (Mg2+):      1.5 mM
dntp:           0.6 mM
dna_conc:       50.0 nM
temp:           37.0 °C
max_loop:       30 bp
alignment_type: Any
dimer:          1 (calculate dimer)
```

## Unit Conversions

- **Input Temperature**: Celsius → Converted to Kelvin internally
- **ΔH**: kcal/mol (from parameters) → cal/mol (in results)
- **ΔS**: cal/mol/K (consistent throughout)
- **ΔG**: kcal/mol (calculated) → cal/mol (in results)
- **DNA Concentration**: nM → M (for calculations)

## Error Handling

Errors are indicated by:
1. `result.temp == THAL_ERROR_SCORE` (f64::NEG_INFINITY)
2. `result.msg` contains error description

Common errors:
- Sequence too short (< 2 bp)
- Sequence too long (> 60 bp for one sequence, > 10000 bp for other)
- Invalid thermodynamic parameters

## Testing

The module includes unit tests:

```bash
cargo test --lib thal
```

Tests cover:
- Nearest neighbor parameters
- Initiation parameters
- Self-complementary detection
- Thermodynamic calculations

## References

- SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95(4): 1460-1465.
- SantaLucia, J. & Hicks, D. (2004). "The thermodynamics of DNA structural motifs." Annual Review of Biophysics and Biomolecular Structure 33: 415-440.
- Primer3 Source: https://github.com/primer3-org/primer3
- Original C Implementation: thal_main.c, thal.h, thal.c

## Future Enhancements

To achieve full primer3 compatibility, the following would need to be implemented:

1. **Dynamic Programming Alignment**
   - Implement the DP algorithm for finding optimal alignments
   - Handle all alignment types (Any, End1, End2, Hairpin)

2. **Structural Features**
   - Dangling end effects (load from parameter files)
   - Mismatch penalties
   - Loop entropy corrections
   - Bulge penalties
   - Internal loop penalties

3. **Parameter Files**
   - Load thermodynamic parameters from primer3_config files
   - Support for custom parameter sets

4. **Performance Optimization**
   - Memoization of calculations
   - SIMD optimizations for alignment search
   - Parallel processing for multiple sequences

## Limitations

This Rust implementation is designed for **perfect match** calculations. For production use requiring full primer3 accuracy with all structural features, consider:

1. Using primer3-py directly
2. Creating FFI bindings to the C library
3. Implementing the full DP algorithm (complex undertaking)

## License

This code is based on primer3, which is licensed under the GNU General Public License v2.0.
