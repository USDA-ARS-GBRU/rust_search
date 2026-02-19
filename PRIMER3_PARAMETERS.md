# Exact Parameters for primer3-py Model Comparison

## API Parameter Mapping

### primer3.bindings.calc_heterodimer() Parameters

| Parameter | Type | Default | Description | Rust Equivalent |
|-----------|------|---------|-------------|-----------------|
| `seq1` | str/bytes | Required | First DNA sequence | Pattern sequence |
| `seq2` | str/bytes | Required | Second DNA sequence | Target sequence |
| `mv_conc` | float | 50.0 | Monovalent cation (Na+) concentration (mM) | `--na` |
| `dv_conc` | float | 1.5 | Divalent cation (Mg2+) concentration (mM) | `--mg` |
| `dntp_conc` | float | 0.6 | dNTP concentration (mM) | `--dntp` |
| `dna_conc` | float | 50.0 | DNA concentration (nM) | `--dnac` |
| `temp_c` | float | 37.0 | Temperature (Celsius) | `--temp` |
| `max_loop` | int | 30 | Maximum loop size | Not implemented in Rust |
| `output_structure` | bool | False | Return ASCII structure | Not implemented in Rust |

## Critical Difference: DNA Concentration Default

**IMPORTANT**: The default DNA concentration differs:
- **primer3-py**: `dna_conc=50.0` nM
- **main_simplesalt.rs**: `--dnac 200.0` nM (default)
- **main_fullsalt.rs**: `--dnac 200.0` nM (default)

To match primer3-py exactly, use: `--dnac 50.0`

## Return Values

### primer3-py ThermoResult Object
```python
result.dg      # ΔG (kcal/mol) at specified temperature
result.tm      # Tm (°C) - melting temperature
result.dh      # ΔH (cal/mol) - enthalpy
result.ds      # ΔS (cal/mol/K) - entropy
result.structure  # ASCII structure (if output_structure=True)
```

**Note**: primer3-py reports ΔH in **cal/mol** (not kcal/mol)

### Rust Output Format
```
seq_id  position  dg      tm      sequence
```
Only outputs ΔG and Tm (not ΔH and ΔS)

**Note**: Rust scripts internally use:
- ΔH in **kcal/mol** (from SantaLucia 1998 parameters)
- ΔS in **cal/mol/K** (from SantaLucia 1998 parameters)
- ΔG in **kcal/mol** (calculated from ΔH and ΔS)

## Thermodynamic Model Details

### SantaLucia 1998 Nearest Neighbor Parameters (All Models)

Both Rust versions and primer3-py use identical NN parameters:

| Pair | ΔH (kcal/mol) | ΔS (cal/mol/K) |
|------|---------------|----------------|
| AA/TT | -7.9 | -22.2 |
| AT | -7.2 | -20.4 |
| TA | -7.2 | -21.3 |
| CA/TG | -8.5 | -22.7 |
| GT/AC | -8.4 | -22.4 |
| CT/AG | -7.8 | -21.0 |
| GA/TC | -8.2 | -22.2 |
| CG | -10.6 | -27.2 |
| GC | -9.8 | -24.4 |
| CC/GG | -8.0 | -19.9 |

### Initiation Parameters

#### main_simplesalt.rs (Fixed)
```
All sequences: ΔH = 0.2 kcal/mol, ΔS = -5.7 cal/mol/K
```

#### main_fullsalt.rs & primer3-py (Terminal Base Pair Dependent)
```
A-T terminal pairs: ΔH = 2.3 kcal/mol, ΔS = 4.1 cal/mol/K
G-C terminal pairs: ΔH = 0.1 kcal/mol, ΔS = -2.8 cal/mol/K
Mixed pairs:        ΔH = 1.2 kcal/mol, ΔS = 0.7 cal/mol/K
Default:            ΔH = 0.2 kcal/mol, ΔS = -5.7 cal/mol/K
```

### Salt Correction (SantaLucia 2004)

All three models use the same formula:

```
[Na+]_eq = [Na+] + 120 * sqrt([Mg2+]_eff)

where [Mg2+]_eff = max(0, [Mg2+] - [dNTP])

ΔS_corrected = ΔS + 0.368 * (N-1) * ln([Na+]_eq / 1000)
```

### Tm Calculation

#### main_simplesalt.rs (Fixed C/4)
```
Tm = ΔH / (ΔS + R * ln(C/4)) - 273.15

where:
  R = 1.9872 cal/(K*mol)
  C = DNA concentration (M)
  ΔH in kcal/mol (converted to cal/mol by multiplying by 1000)
  ΔS in cal/(mol*K)
```

#### main_fullsalt.rs & primer3-py (Symmetry-Aware)
```
For self-complementary sequences:
  Tm = ΔH / (ΔS + R * ln(C/2)) - 273.15

For heterodimers:
  Tm = ΔH / (ΔS + R * ln(C/4)) - 273.15
```

## Features Implemented in Each Model

### main_simplesalt.rs
- ✓ SantaLucia 1998 NN parameters
- ✓ Fixed initiation penalty
- ✓ SantaLucia 2004 salt correction
- ✓ Basic Tm calculation (C/4)
- ✗ Terminal base pair-dependent initiation
- ✗ Symmetry correction
- ✗ Dangling ends
- ✗ Mismatches
- ✗ Loops/bulges
- ✗ Optimal alignment search

### main_fullsalt.rs
- ✓ SantaLucia 1998 NN parameters
- ✓ Terminal base pair-dependent initiation
- ✓ SantaLucia 2004 salt correction
- ✓ Symmetry-aware Tm calculation
- ✗ Dangling ends
- ✗ Mismatches
- ✗ Loops/bulges
- ✗ Optimal alignment search

### primer3-py calc_heterodimer()
- ✓ SantaLucia 1998 NN parameters
- ✓ Terminal base pair-dependent initiation
- ✓ SantaLucia 2004 salt correction
- ✓ Symmetry-aware Tm calculation
- ✓ Dangling ends
- ✓ Mismatches
- ✓ Loops/bulges
- ✓ Optimal alignment search
- ✓ Maximum loop size constraint

## Expected Differences

### main_simplesalt vs main_fullsalt
- **Magnitude**: 0-2 kcal/mol (usually < 1 kcal/mol)
- **Cause**: Terminal base pair-dependent initiation
- **Visibility**: May not show at 2 decimal place precision
- **Tm difference**: Mainly for self-complementary sequences

### main_fullsalt vs primer3-py
- **Magnitude**: 1-5 kcal/mol (depends on sequence)
- **Causes**:
  - Dangling end effects (0.5-1.5 kcal/mol)
  - Mismatch penalties (if any mismatches)
  - Loop entropy corrections (if any loops)
  - Bulge penalties (if any bulges)
  - Optimal alignment search (primer3-py considers all alignments)
- **Visibility**: Usually visible at 2 decimal place precision

## How to Use for Comparison

### Command Line Examples

```bash
# Test sequence
SEQ="ATGCGATCGATCG"

# Create test files
echo ">test" > test.fasta
echo "$SEQ" >> test.fasta

# Test with primer3-py (requires installation)
python3 << 'EOF'
import primer3
result = primer3.bindings.calc_heterodimer(
    "ATGCGATCGATCG", "ATGCGATCGATCG",
    mv_conc=50.0,
    dv_conc=1.5,
    dntp_conc=0.6,
    dna_conc=50.0,
    temp_c=37.0
)
print(f"ΔG: {result.dg:.2f} kcal/mol")
print(f"Tm: {result.tm:.2f} °C")
print(f"ΔH: {result.dh:.2f} kcal/mol")
print(f"ΔS: {result.ds:.2f} cal/mol/K")
EOF

# Test with main_simplesalt
cargo run --release --bin main_simplesalt -- \
  --file test.fasta --patterns test.fasta \
  --na 50.0 --mg 1.5 --dntp 0.6 --dnac 50.0 --temp 37.0 \
  --threshold 0

# Test with main_fullsalt
cargo run --release --bin main_fullsalt -- \
  --file test.fasta --patterns test.fasta \
  --na 50.0 --mg 1.5 --dntp 0.6 --dnac 50.0 --temp 37.0 \
  --threshold 0
```

### Using the Comparison Script

```bash
# Make script executable
chmod +x compare_models.py

# Run comparison
python3 compare_models.py
```

This will test multiple sequences with all three models and show differences.

## Installation of primer3-py

```bash
# Using pip
pip3 install primer3-py

# Or using conda
conda install -c bioconda primer3-py
```

## References

- SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95(4): 1460-1465.
- SantaLucia, J. & Hicks, D. (2004). "The thermodynamics of DNA structural motifs." Annual Review of Biophysics and Biomolecular Structure 33: 415-440.
- Primer3: https://github.com/primer3-org/primer3
- primer3-py: https://github.com/libnano/primer3-py
