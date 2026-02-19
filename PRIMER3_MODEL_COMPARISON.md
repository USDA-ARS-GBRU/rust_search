# Primer3 Model Comparison: main_simplesalt.rs vs main_fullsalt.rs vs primer3-py

## API Comparison

### primer3-py calc_heterodimer()
```python
primer3.bindings.calc_heterodimer(
    seq1,                    # First DNA sequence
    seq2,                    # Second DNA sequence (target)
    mv_conc=50.0,           # Monovalent cation (Na+) concentration (mM)
    dv_conc=1.5,            # Divalent cation (Mg2+) concentration (mM)
    dntp_conc=0.6,          # dNTP concentration (mM)
    dna_conc=50.0,          # DNA concentration (nM)
    temp_c=37.0,            # Temperature (Celsius)
    max_loop=30,            # Maximum loop size
    output_structure=False  # Return ASCII structure
)
```

### Rust Implementation Parameters

Both `main_simplesalt.rs` and `main_fullsalt.rs` use:
```rust
--na 50.0          # Monovalent salt Na+ (mM)
--mg 1.5           # Divalent salt Mg2+ (mM)
--dntp 0.6         # dNTPs (mM)
--dnac 200.0       # Primer concentration (nM) - NOTE: Different default!
--temp 37.0        # Temperature (C)
```

## Key Differences

### 1. DNA Concentration Parameter
- **primer3-py**: `dna_conc=50.0` nM (default)
- **main_simplesalt.rs**: `--dnac 200.0` nM (default)
- **main_fullsalt.rs**: `--dnac 200.0` nM (default)

**To match primer3-py exactly, use: `--dnac 50.0`**

### 2. Thermodynamic Model Components

#### main_simplesalt.rs (Simplified)
- Fixed initiation penalty: ΔH = 0.2 kcal/mol, ΔS = -5.7 cal/mol/K
- SantaLucia 1998 nearest neighbor parameters
- Basic salt correction: [Na+]_eq = [Na+] + 120 * sqrt([Mg2+] - [dNTP])
- Tm calculation: Always uses C/4 (heterodimer formula)

#### main_fullsalt.rs (Full Model)
- **Terminal base pair-dependent initiation**:
  - A-T pairs: ΔH = 2.3, ΔS = 4.1
  - G-C pairs: ΔH = 0.1, ΔS = -2.8
  - Mixed pairs: ΔH = 1.2, ΔS = 0.7
  - Default: ΔH = 0.2, ΔS = -5.7
- SantaLucia 1998 nearest neighbor parameters
- Full salt correction: [Na+]_eq = [Na+] + 120 * sqrt([Mg2+] - [dNTP])
- **Symmetry-aware Tm calculation**:
  - Self-complementary sequences: C/2
  - Heterodimers: C/4

#### primer3-py (Full Model + Additional Features)
- **Terminal base pair-dependent initiation** (same as main_fullsalt)
- SantaLucia 1998 nearest neighbor parameters
- **Full salt correction with additional terms**:
  - Includes Mg2+ and dNTP binding effects
  - More sophisticated handling of ionic strength
- **Symmetry-aware Tm calculation** (same as main_fullsalt)
- **Additional features NOT in Rust versions**:
  - Dangling end effects
  - Mismatch penalties
  - Loop entropy corrections
  - Bulge penalties
  - Internal loop penalties
  - Maximum loop size constraint (max_loop parameter)

## To Create Perfect Comparison

### Step 1: Match DNA Concentration
Use `--dnac 50.0` instead of default 200.0 to match primer3-py default

### Step 2: Test with Same Sequences
```bash
# Test sequence
SEQ="ATGCGATCGATCG"

# primer3-py
python3 -c "
import primer3
result = primer3.bindings.calc_heterodimer(
    '$SEQ', '$SEQ',
    mv_conc=50.0,
    dv_conc=1.5,
    dntp_conc=0.6,
    dna_conc=50.0,
    temp_c=37.0
)
print(f'ΔG: {result.dg:.2f}')
print(f'Tm: {result.tm:.2f}')
"

# Rust main_simplesalt
cargo run --bin main_simplesalt --release -- \
  --file test.fasta --patterns test.fasta \
  --na 50.0 --mg 1.5 --dntp 0.6 --dnac 50.0 --temp 37.0

# Rust main_fullsalt
cargo run --bin main_fullsalt --release -- \
  --file test.fasta --patterns test.fasta \
  --na 50.0 --mg 1.5 --dntp 0.6 --dnac 50.0 --temp 37.0
```

### Step 3: Expected Differences

**main_simplesalt vs main_fullsalt:**
- Small differences (~0-2 kcal/mol) due to terminal effects
- Differences in Tm for self-complementary sequences
- Differences may not be visible at 2 decimal place precision

**main_fullsalt vs primer3-py:**
- Differences due to:
  - Dangling end effects (not implemented in Rust)
  - Mismatch penalties (not implemented in Rust)
  - Loop entropy corrections (not implemented in Rust)
  - Bulge and internal loop penalties (not implemented in Rust)
- These can cause 1-5 kcal/mol differences depending on sequence

## Implementation Notes

### What primer3-py calc_heterodimer() Does
1. Performs full thermodynamic calculation using SantaLucia 1998 + 2004
2. Considers all possible alignments of the two sequences
3. Calculates optimal binding structure
4. Includes dangling ends, mismatches, loops, bulges
5. Returns ΔG at specified temperature and Tm

### What main_fullsalt.rs Does
1. Calculates thermodynamics for exact sequence match
2. Uses terminal base pair-dependent initiation
3. Includes symmetry correction for Tm
4. Does NOT consider:
   - Dangling ends
   - Mismatches
   - Loops/bulges
   - Optimal alignment (assumes perfect match)

### What main_simplesalt.rs Does
1. Same as main_fullsalt but with fixed initiation penalty
2. Simpler model, faster calculation
3. Less accurate for sequences with different terminal bases

## Unit Clarification

### Internal Units Used

**primer3-py:**
- ΔH: **cal/mol** (reported)
- ΔS: **cal/mol/K** (reported)
- ΔG: **kcal/mol** (reported)

**Rust scripts (main_simplesalt.rs and main_fullsalt.rs):**
- ΔH: **kcal/mol** (internally, from SantaLucia 1998 parameters)
- ΔS: **cal/mol/K** (internally, from SantaLucia 1998 parameters)
- ��G: **kcal/mol** (calculated and reported)

**Conversion factor:** 1 kcal/mol = 1000 cal/mol

To convert Rust ΔH to match primer3-py format:
```
ΔH_cal_per_mol = ΔH_kcal_per_mol × 1000
```

## To Achieve Perfect Match with primer3-py

You would need to implement:
1. ✓ Terminal base pair-dependent initiation (done in main_fullsalt)
2. ✓ Symmetry correction (done in main_fullsalt)
3. ✗ Dangling end effects (NOT implemented)
4. ✗ Mismatch penalties (NOT implemented)
5. ✗ Loop entropy corrections (NOT implemented)
6. ✗ Bulge penalties (NOT implemented)
7. ✗ Internal loop penalties (NOT implemented)
8. ✗ Optimal alignment search (NOT implemented)

The Rust versions calculate thermodynamics for a **perfect match** between primer and target, while primer3-py considers **all possible alignments** and structural features.
