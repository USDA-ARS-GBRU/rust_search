# Rust Primer3 Implementation - Architecture Overview

## Project Structure

```
rust_search/
├── Cargo.toml                          # Project manifest
├── src/
│   ├── lib.rs                          # Core thal module (library)
│   ├── main.rs                         # (unused)
│   └── bin/
│       ├── main_simplesalt.rs          # Simplified model (original)
│       ├── main_fullsalt.rs            # Full model using thal module
│       └── main_bruteforce_simplesalt.rs
├── THAL_RUST_IMPLEMENTATION.md         # Detailed thal documentation
├── FULLSALT_IMPLEMENTATION.md          # Full model documentation
├── PRIMER3_PARAMETERS.md               # Parameter reference
├── PRIMER3_MODEL_COMPARISON.md         # Model comparison
└── compare_models_corrected.py         # Comparison script
```

## Module Architecture

### src/lib.rs - Core Thal Module

This is the heart of the implementation. It provides:

**Public API:**
```rust
pub mod thal {
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
    
    // Constants
    pub const THAL_MAX_ALIGN: usize = 60
    pub const THAL_MAX_SEQ: usize = 10000
    pub const MAX_LOOP: i32 = 30
    pub const ABSOLUTE_ZERO: f64 = 273.15
}
```

### src/bin/main_fullsalt.rs - Binary Implementation

Uses the thal module to:
1. Parse command-line arguments
2. Read primer and genome sequences
3. Find primer matches using Aho-Corasick
4. Calculate thermodynamics using `thal::thal()`
5. Output results

**Key Features:**
- Parallel processing with rayon
- Chunked genome processing
- Threshold-based filtering
- Matches primer3-py defaults

## Comparison of Implementations

### main_simplesalt.rs (Original)
- **Model**: Simplified SantaLucia 1998
- **Features**: Fixed initiation penalty, basic salt correction
- **Accuracy**: ~2-3 kcal/mol difference from primer3-py
- **Speed**: Fast
- **Code**: Inline calculations

### main_fullsalt.rs (New)
- **Model**: Full SantaLucia 1998 + 2004
- **Features**: Terminal-dependent initiation, full salt correction, symmetry-aware Tm
- **Accuracy**: ~1-5 kcal/mol difference from primer3-py
- **Speed**: Fast (same as simplesalt)
- **Code**: Uses thal module

### primer3-py (Reference)
- **Model**: Full libprimer3 implementation
- **Features**: All structural features (dangling ends, mismatches, loops, bulges)
- **Accuracy**: Reference implementation
- **Speed**: Slower (C library with DP algorithm)
- **Code**: Cython bindings to C

## Key Improvements in New Architecture

### 1. **Modularity**
- Core thermodynamic logic in `lib.rs`
- Can be imported by other binaries
- Easy to test and maintain

### 2. **Reusability**
- `thal` module can be used by multiple binaries
- Easy to create new analysis tools
- Clear API for thermodynamic calculations

### 3. **Maintainability**
- Centralized parameter definitions
- Consistent error handling
- Well-documented functions

### 4. **Extensibility**
- Easy to add new alignment types
- Can implement full DP algorithm later
- Can add parameter file loading

### 5. **Testing**
- Unit tests in lib.rs
- Easy to test individual functions
- Clear test cases

## Usage Examples

### As a Library

```rust
use rust_search::thal::{self, ThalArgs, ThalAlignmentType, ThalMode, ABSOLUTE_ZERO};

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
}
```

### As a Binary

```bash
# Build
cargo build --release --bin main_fullsalt

# Run
cargo run --release --bin main_fullsalt -- \
  --file genome.fasta --patterns primers.fasta \
  --na 50.0 --mg 1.5 --dntp 0.6 --dnac 50.0 --temp 37.0
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

## Thermodynamic Model Details

### Calculation Steps

1. **Initiation**: Terminal base pair-dependent penalty
2. **Nearest Neighbor Sum**: Sum of all base pair stacking energies
3. **Salt Correction**: Adjust entropy for ionic strength
4. **ΔG Calculation**: ΔG = ΔH - T*ΔS
5. **Tm Calculation**: Tm = ΔH / (ΔS + R*ln(C/factor))

### Unit Conversions

- **Input**: Celsius → Kelvin (add 273.15)
- **ΔH**: kcal/mol → cal/mol (multiply by 1000)
- **ΔS**: cal/mol/K (unchanged)
- **ΔG**: kcal/mol → cal/mol (multiply by 1000)
- **Output**: Kelvin → Celsius (subtract 273.15)

## Expected Differences from primer3-py

After unit conversion:
- **Magnitude**: 1-5 kcal/mol
- **Causes**: Missing dangling ends, mismatches, loops, bulges
- **When minimal**: Perfect matches, short sequences, no secondary structure

## Future Enhancements

### Phase 1: Current Implementation
✓ Core thermodynamic calculations
✓ Terminal-dependent initiation
✓ Full salt correction
✓ Symmetry-aware Tm

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

## Building and Testing

```bash
# Build library and binaries
cargo build --release

# Run tests
cargo test --lib thal

# Run main_fullsalt
cargo run --release --bin main_fullsalt -- \
  --file test.fasta --patterns test.fasta

# Compare with primer3-py
python3 compare_models_corrected.py
```

## Documentation Files

1. **THAL_RUST_IMPLEMENTATION.md** - Detailed thal module documentation
2. **FULLSALT_IMPLEMENTATION.md** - Full model implementation details
3. **PRIMER3_PARAMETERS.md** - Parameter reference and API mapping
4. **PRIMER3_MODEL_COMPARISON.md** - Model comparison and differences
5. **UNIT_CONVERSION_ANALYSIS.md** - Unit conversion explanation

## References

- Primer3 Source: https://github.com/primer3-org/primer3
- primer3-py: https://github.com/libnano/primer3-py
- SantaLucia 1998: "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics"
- SantaLucia 2004: "The thermodynamics of DNA structural motifs"

## License

Based on primer3, licensed under GNU General Public License v2.0.
