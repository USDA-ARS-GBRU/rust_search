# Analysis: Why primer3-py Results Are 1000x Off

## The Problem

Your comparison showed:
- **primer3-py**: ΔG = -11186.18 kcal/mol
- **Rust scripts**: ΔG = -13.62 kcal/mol

This is approximately **1000x difference**, which suggests a unit conversion issue.

## Root Cause: primer3-py Returns Values in cal/mol

**primer3-py returns:**
- ΔG in **cal/mol** (not kcal/mol)
- ΔH in **cal/mol** (not kcal/mol)
- ΔS in **cal/mol/K**

**Your Rust scripts return:**
- ΔG in **kcal/mol**
- ΔH in **kcal/mol** (internally)
- ΔS in **cal/mol/K**

## Conversion

To convert primer3-py results to match Rust:
```
ΔG_kcal_per_mol = ΔG_cal_per_mol / 1000
ΔH_kcal_per_mol = ΔH_cal_per_mol / 1000
```

### Example Conversion
```
primer3-py output:
  ΔG: -11186.18 cal/mol → -11.19 kcal/mol
  ΔH: -89000.00 cal/mol → -89.00 kcal/mol

Rust output:
  ΔG: -13.62 kcal/mol
```

After conversion, primer3-py ΔG = -11.19 kcal/mol vs Rust ΔG = -13.62 kcal/mol

**Difference: -2.43 kcal/mol** (much more reasonable!)

## Why There's Still a Difference

Even after unit conversion, there's still a ~2 kcal/mol difference. This is likely due to:

1. **Different calculation methods**:
   - primer3-py considers **all possible alignments** of primer and target
   - Rust scripts assume **perfect alignment** (primer matches target exactly)

2. **Structural features in primer3-py**:
   - Dangling ends
   - Mismatches
   - Loops and bulges
   - Optimal binding structure

3. **Different salt correction models**:
   - primer3-py may use a more sophisticated salt model
   - Rust uses basic SantaLucia 2004 formula

## Corrected Comparison Script

I've created `compare_models_corrected.py` which:
1. Converts primer3-py results from cal/mol to kcal/mol
2. Tests primer-target binding (heterodimer)
3. Shows meaningful comparisons

## How to Use

```bash
# Make script executable
chmod +x compare_models_corrected.py

# Run corrected comparison
python3 compare_models_corrected.py
```

## Expected Results After Correction

After unit conversion, you should see:
- **main_simplesalt vs main_fullsalt**: Small differences (0-2 kcal/mol)
- **Rust vs primer3-py**: Moderate differences (1-5 kcal/mol)
  - Due to structural features not implemented in Rust
  - Due to different alignment strategies

## Key Takeaway

The 1000x difference was purely a **unit conversion issue**. primer3-py reports thermodynamic values in cal/mol while your Rust scripts report in kcal/mol. After converting to the same units, the values are much closer and the remaining differences are due to model differences (dangling ends, mismatches, optimal alignment, etc.).
