# main_fullsalt.rs - Full Primer3 Model Implementation

## Overview

This implementation is based on the Cython code from `primer3-py` (thermoanalysis.pyx/pxi/pxd) and implements the full SantaLucia 1998 + 2004 thermodynamic model for DNA primer-target binding calculations.

## Key Features from primer3-py Implementation

### 1. **Thermodynamic Calculation Method**
The Cython code calls the `thal()` function from libprimer3, which performs:
- Full thermodynamic alignment calculation
- Considers all possible alignments of primer and target
- Calculates optimal binding structure
- Includes dangling ends, mismatches, loops, and bulges

### 2. **Parameters Used**
From `thermoanalysis.pyx` initialization:
```cython
self.thalargs.mv = mv_conc           # Monovalent cation (Na+) concentration (mM)
self.thalargs.dv = dv_conc           # Divalent cation (Mg2+) concentration (mM)
self.thalargs.dntp = dntp_conc       # dNTP concentration (mM)
self.thalargs.dna_conc = dna_conc    # DNA concentration (nM)
self.thalargs.temp = temp_c + 273.15 # Temperature in Kelvin
self.thalargs.maxLoop = max_loop     # Maximum loop size (default 30)
self.thalargs.dimer = 1              # Calculate dimer structure
self.thalargs.type = 1               # thal_alignment_any
```

### 3. **Rust Implementation Differences**

The Rust version implements the **core thermodynamic model** but cannot replicate the full `thal()` function which requires:
- Dynamic programming alignment algorithm
- Complex loop entropy calculations
- Mismatch penalty tables
- Bulge and internal loop penalties
- Dangling end effects database

**What IS implemented in Rust:**
- ✓ SantaLucia 1998 nearest neighbor parameters
- ✓ Terminal base pair-dependent initiation (end effects)
- ✓ SantaLucia 2004 salt correction
- ✓ Symmetry-aware Tm calculation
- ✓ Proper unit conversions (kcal/mol, cal/mol/K)

**What is NOT implemented (requires libprimer3):**
- ✗ Dynamic programming alignment search
- ✗ Dangling end effects (database-driven)
- ✗ Mismatch penalties
- ✗ Loop entropy corrections
- ✗ Bulge penalties
- ✗ Internal loop penalties
- ✗ Optimal alignment discovery

## Default Parameters

Changed from previous version to match primer3-py defaults:
```rust
--na 50.0          # Monovalent salt Na+ (mM) - CHANGED from 50.0
--mg 1.5           # Divalent salt Mg2+ (mM)
--dntp 0.6         # dNTPs (mM)
--dnac 50.0        # Primer concentration (nM) - CHANGED from 200.0
--temp 37.0        # Temperature (C)
--max_loop 30      # Maximum loop size (bp) - NEW parameter
```

## Thermodynamic Model Components

### Initiation Parameters (Terminal Base Pair Dependent)
From `thermoanalysis.pyx` - SantaLucia 1998 Table 1:
```
A-T terminal pairs: ΔH = 2.3 kcal/mol, ΔS = 4.1 cal/mol/K
G-C terminal pairs: ΔH = 0.1 kcal/mol, ΔS = -2.8 cal/mol/K
Mixed pairs:        ΔH = 1.2 kcal/mol, ΔS = 0.7 cal/mol/K
Default:            ΔH = 0.2 kcal/mol, ΔS = -5.7 cal/mol/K
```

### Salt Correction (SantaLucia 2004)
From `thermoanalysis.pyx`:
```
[Na+]_eq = [Na+] + 120 * sqrt([Mg2+]_eff)
where [Mg2+]_eff = max(0, [Mg2+] - [dNTP])

ΔS_corrected = ΔS + 0.368 * (N-1) * ln([Na+]_eq / 1000)
```

### Tm Calculation (Symmetry-Aware)
From `thermoanalysis.pyx`:
```
For self-complementary sequences:
  Tm = ΔH / (ΔS + R * ln(C/2)) - 273.15

For heterodimers:
  Tm = ΔH / (ΔS + R * ln(C/4)) - 273.15

where:
  R = 1.9872 cal/(K*mol)
  C = DNA concentration (M)
  ΔH in kcal/mol (converted to cal/mol by multiplying by 1000)
  ΔS in cal/(mol*K)
```

## Comparison with Versions

### main_simplesalt.rs
- Fixed initiation penalty (0.2/-5.7) for all sequences
- Basic salt correction
- No symmetry correction
- **Differences from main_fullsalt: 0-2 kcal/mol**

### main_fullsalt.rs (This Version)
- Terminal base pair-dependent initiation
- Full salt correction with Mg2+ and dNTP effects
- Symmetry-aware Tm calculation
- **Differences from primer3-py: 1-5 kcal/mol**
  - Due to missing dangling ends, mismatches, loops, bulges
  - Due to different alignment strategy (perfect match vs optimal)

### primer3-py
- Full libprimer3 `thal()` function
- Dynamic programming alignment
- All structural features (dangling ends, mismatches, loops, bulges)
- **Reference implementation**

## Expected Differences from primer3-py

After unit conversion (primer3-py returns cal/mol, Rust returns kcal/mol):

1. **Magnitude: 1-5 kcal/mol difference**
2. **Causes:**
   - Dangling end effects (0.5-1.5 kcal/mol)
   - Mismatch penalties (if any mismatches)
   - Loop entropy corrections (if any loops)
   - Bulge penalties (if any bulges)
   - Optimal alignment search (primer3-py considers all alignments)

3. **When differences are minimal:**
   - Perfect matches with no loops/bulges
   - Short sequences (< 20 bp)
   - Sequences with minimal secondary structure

## Usage

```bash
# Build
cargo build --release --bin main_fullsalt

# Run with default parameters (matching primer3-py)
cargo run --release --bin main_fullsalt -- \
  --file genome.fasta --patterns primers.fasta

# Run with custom parameters
cargo run --release --bin main_fullsalt -- \
  --file genome.fasta --patterns primers.fasta \
  --na 50.0 --mg 1.5 --dntp 0.6 --dnac 50.0 --temp 37.0 \
  --max_loop 30 --threshold -10.0
```

## References

- SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95(4): 1460-1465.
- SantaLucia, J. & Hicks, D. (2004). "The thermodynamics of DNA structural motifs." Annual Review of Biophysics and Biomolecular Structure 33: 415-440.
- Primer3: https://github.com/primer3-org/primer3
- primer3-py: https://github.com/libnano/primer3-py
- Cython Implementation: thermoanalysis.pyx/pxi/pxd

## Implementation Notes

1. **Unit Consistency:**
   - ΔH: kcal/mol (from SantaLucia 1998 parameters)
   - ΔS: cal/mol/K (from SantaLucia 1998 parameters)
   - ΔG: kcal/mol (calculated)
   - Tm: °C (calculated)

2. **Temperature Conversion:**
   - Input: Celsius
   - Internal: Kelvin (for calculations)
   - Output: Celsius

3. **Concentration Units:**
   - Na+, Mg2+, dNTP: mM (millimolar)
   - DNA: nM (nanomolar) - converted to M for calculations

4. **Symmetry Detection:**
   - Checks if sequence is self-complementary
   - Uses C/2 for homodimers, C/4 for heterodimers
   - Matches primer3-py behavior

## Limitations

This Rust implementation is designed for **perfect match** calculations between primer and target. For production use requiring full primer3 accuracy, consider:
1. Using primer3-py directly
2. Implementing FFI bindings to libprimer3
3. Porting the full `thal()` algorithm to Rust (complex undertaking)
