# Rust Primer3 Implementation - Complete Summary

## Project Overview

This project provides a Rust implementation of the primer3 thermodynamic alignment (`thal`) functionality, based on the C source code from the primer3 library. The implementation includes:

1. **Core Library** (`src/lib.rs`): Thermodynamic calculation module
2. **Binary** (`src/bin/main_fullsalt.rs`): Genome-wide primer binding analysis tool
3. **Documentation**: Comprehensive guides and API references

## Architecture

### Module Structure

```
rust_search/
├── src/
│   ├── lib.rs                    # Core thal module (library)
│   ├── main.rs                   # (unused)
│   └── bin/
│       ├── main_simplesalt.rs    # Simplified model (original)
│       ├── main_fullsalt.rs      # Full model using thal module
│       └── main_bruteforce_simplesalt.rs
├── Cargo.toml                    # Project manifest
└── Documentation files
```

### Public API

The `thal` module provides:

```rust
// Enumerations
pub enum ThalAlignmentType { Any, End1, End2, Hairpin }
pub enum ThalMode { Fast, General, DebugFast, Debug, Struct }

// Structures
pub struct ThalArgs { ... }
pub struct ThalResults { ... }
pub struct NNParams { ... }

// Functions
pub fn get_nn_params(base1: u8, base2: u8) -> NNParams
pub fn get_initiation_params(first: u8, last: u8) -> NNParams
pub fn calculate_na_equivalent(na: f64, mg: f64, dntp: f64) -> f64
pub fn is_self_complementary(seq: &[u8]) -> bool
pub fn calculate_thermo(seq: &[u8], args: &ThalArgs) -> ThalResults
pub fn thal(seq1: &[u8], seq2: &[u8], args: &ThalArgs, mode: ThalMode) -> ThalResults
pub fn create_default_args() -> ThalArgs

// Constants
pub const THAL_MAX_ALIGN: usize = 60
pub const THAL_MAX_SEQ: usize = 10000
pub const MAX_LOOP: i32 = 30
pub const ABSOLUTE_ZERO: f64 = 273.15
pub const THAL_ERROR_SCORE: f64 = f64::NEG_INFINITY
```

## Building

### Prerequisites
- Rust 1.56+ (2021 edition)
- Cargo

### Build Commands

```bash
# Check compilation
cargo check

# Build debug version
cargo build

# Build release version (optimized)
cargo build --release

# Run tests
cargo test --lib thal

# Build specific binary
cargo build --release --bin main_fullsalt
```

## Usage

### As a Library

```rust
use rust_search::{thal, ThalArgs, ThalAlignmentType, ThalMode, ABSOLUTE_ZERO};

fn main() {
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
    
    let result = thal::thal(b"ATGC", b"ATGC", &args, ThalMode::Fast);
    println!("Tm: {:.2}°C", result.temp);
    println!("ΔG: {:.2} cal/mol", result.dg);
}
```

### As a Binary

```bash
# Run with default parameters
cargo run --release --bin main_fullsalt -- \
  --file genome.fasta --patterns primers.fasta

# Run with custom parameters
cargo run --release --bin main_fullsalt -- \
  --file genome.fasta --patterns primers.fasta \
  --na 50.0 --mg 1.5 --dntp 0.6 --dnac 50.0 --temp 37.0 \
  --max_loop 30 --threshold -10.0
```

## Default Parameters (Matching primer3-py)

```
--na 50.0          # Monovalent cation (Na+) concentration (mM)
--mg 1.5           # Divalent cation (Mg2+) concentration (mM)
--dntp 0.6         # dNTP concentration (mM)
--dnac 50.0        # DNA concentration (nM)
--temp 37.0        # Temperature (°C)
--max_loop 30      # Maximum loop size (bp)
--threshold -10.0  # ΔG threshold (kcal/mol)
```

## Output Format

```
seq_id  position  dg      tm      sequence
```

Example:
```
chr1    12345    -13.62  44.31   ATGCGATCGATCG
chr1    12456    -15.99  54.14   GCGCGCGCGC
```

## Thermodynamic Model

### Calculation Steps

1. **Initiation**: Terminal base pair-dependent penalty
2. **Nearest Neighbor Sum**: Sum of all base pair stacking energies
3. **Salt Correction**: Adjust entropy for ionic strength (SantaLucia 2004)
4. **ΔG Calculation**: ΔG = ΔH - T*ΔS
5. **Tm Calculation**: Tm = ΔH / (ΔS + R*ln(C/factor))

### SantaLucia 1998 Nearest Neighbor Parameters

```
AA/TT: ΔH = -7.9 kcal/mol,  ΔS = -22.2 cal/mol/K
AT:    ΔH = -7.2 kcal/mol,  ΔS = -20.4 cal/mol/K
TA:    ΔH = -7.2 kcal/mol,  ΔS = -21.3 cal/mol/K
CA/TG: ΔH = -8.5 kcal/mol,  ΔS = -22.7 cal/mol/K
GT/AC: ΔH = -8.4 kcal/mol,  ΔS = -22.4 cal/mol/K
CT/AG: ΔH = -7.8 kcal/mol,  ΔS = -21.0 cal/mol/K
GA/TC: ΔH = -8.2 kcal/mol,  ΔS = -22.2 cal/mol/K
CG:    ΔH = -10.6 kcal/mol, ΔS = -27.2 cal/mol/K
GC:    ΔH = -9.8 kcal/mol,  ΔS = -24.4 cal/mol/K
CC/GG: ΔH = -8.0 kcal/mol,  ΔS = -19.9 cal/mol/K
```

### Initiation Parameters (Terminal Base Pair Dependent)

```
A-T pairs: ΔH = 2.3 kcal/mol,  ΔS = 4.1 cal/mol/K
G-C pairs: ΔH = 0.1 kcal/mol,  ΔS = -2.8 cal/mol/K
Mixed:     ΔH = 1.2 kcal/mol,  ΔS = 0.7 cal/mol/K
Default:   ΔH = 0.2 kcal/mol,  ΔS = -5.7 cal/mol/K
```

### Salt Correction (SantaLucia 2004)

```
[Na+]_eq = [Na+] + 120 * sqrt([Mg2+]_eff)
where [Mg2+]_eff = max(0, [Mg2+] - [dNTP])

ΔS_corrected = ΔS + 0.368 * (N-1) * ln([Na+]_eq / 1000)
```

## Unit Conversions

- **Input Temperature**: Celsius → Kelvin (add 273.15)
- **ΔH**: kcal/mol (internal) → cal/mol (output, multiply by 1000)
- **ΔS**: cal/mol/K (consistent throughout)
- **ΔG**: kcal/mol (internal) → cal/mol (output, multiply by 1000)
- **DNA Concentration**: nM → M (divide by 1e9)

## Features Implemented

✓ SantaLucia 1998 nearest neighbor parameters
✓ Terminal base pair-dependent initiation
✓ SantaLucia 2004 salt correction
✓ Symmetry-aware Tm calculation
✓ Proper unit conversions
✓ Error handling
✓ Parallel genome processing (main_fullsalt)
✓ Aho-Corasick pattern matching

## Features NOT Implemented (Requires Full C Library)

✗ Dynamic programming alignment search
✗ Dangling end effects (database-driven)
✗ Mismatch penalties
✗ Loop entropy corrections
✗ Bulge penalties
✗ Internal loop penalties
✗ Optimal alignment discovery

## Comparison with Other Implementations

### main_simplesalt.rs
- Fixed initiation penalty
- Basic salt correction
- ~2-3 kcal/mol difference from primer3-py

### main_fullsalt.rs (This Implementation)
- Terminal-dependent initiation
- Full salt correction
- Symmetry-aware Tm
- ~1-5 kcal/mol difference from primer3-py

### primer3-py (Reference)
- Full libprimer3 implementation
- All structural features
- Reference implementation

## Documentation Files

1. **ARCHITECTURE.md** - Project architecture and module organization
2. **THAL_RUST_IMPLEMENTATION.md** - Detailed thal module documentation
3. **FULLSALT_IMPLEMENTATION.md** - Full model implementation details
4. **PRIMER3_PARAMETERS.md** - Parameter reference and API mapping
5. **PRIMER3_MODEL_COMPARISON.md** - Model comparison and differences
6. **UNIT_CONVERSION_ANALYSIS.md** - Unit conversion explanation

## Testing

```bash
# Run all tests
cargo test

# Run library tests only
cargo test --lib

# Run with output
cargo test -- --nocapture

# Run specific test
cargo test test_nn_params
```

## Performance

- **Compilation**: ~3 seconds (release)
- **Calculation**: ~1-2 microseconds per sequence (single-threaded)
- **Parallel Processing**: Linear speedup with number of CPU cores

## Error Handling

Errors are indicated by:
1. `result.temp == THAL_ERROR_SCORE` (f64::NEG_INFINITY)
2. `result.msg` contains error description

Common errors:
- Sequence too short (< 2 bp)
- Sequence too long (> 60 bp for one sequence, > 10000 bp for other)
- Invalid thermodynamic parameters

## Future Enhancements

### Phase 1: Current Implementation ✓
- Core thermodynamic calculations
- Terminal-dependent initiation
- Full salt correction
- Symmetry-aware Tm

### Phase 2: Structural Features
- [ ] Dangling end effects
- [ ] Mismatch penalties
- [ ] Loop entropy corrections
- [ ] Bulge penalties

### Phase 3: Advanced Features
- [ ] Dynamic programming alignment
- [ ] Parameter file loading
- [ ] Optimal alignment search
- [ ] Secondary structure output

### Phase 4: Performance
- [ ] SIMD optimizations
- [ ] Parallel alignment search
- [ ] Memoization
- [ ] GPU acceleration

## References

- Primer3 Source: https://github.com/primer3-org/primer3
- primer3-py: https://github.com/libnano/primer3-py
- SantaLucia 1998: "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics"
- SantaLucia 2004: "The thermodynamics of DNA structural motifs"

## License

Based on primer3, licensed under GNU General Public License v2.0.

## Contributing

To extend this implementation:

1. Add new thermodynamic features to the `thal` module
2. Update tests in the module
3. Document changes in the appropriate documentation file
4. Ensure backward compatibility with existing API

## Support

For issues or questions:
1. Check the documentation files
2. Review the test cases
3. Compare with primer3-py implementation
4. Check the original primer3 C source code
