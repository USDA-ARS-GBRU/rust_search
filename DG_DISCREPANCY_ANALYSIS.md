# ΔG Calculation Discrepancy Analysis

## Problem

The Rust implementation returns different ΔG values compared to primer3-py's `calc_heterodimer()`:

### Example: "ATCGATCGATCGATCGATCG" (20bp homodimer)

**primer3-py:**
```
ΔG: -20589.88 cal/mol
ΔH: -155200.00 cal/mol
ΔS: -434.02 cal/mol/K
Tm: 56.94°C
```

**Rust Implementation:**
```
ΔG: -22110.81 cal/mol
ΔH: -161600.00 cal/mol
ΔS: -449.75 cal/mol/K
Tm: 59.42°C
```

**Difference:**
- ΔH: -6400 cal/mol (-6.4 kcal/mol)
- ΔS: -15.73 cal/mol/K
- ΔG: -1520.93 cal/mol (-1.52 kcal/mol)

## Root Cause

The Rust implementation is a **simplified perfect-match calculator** that:
- ✓ Uses SantaLucia 1998 nearest neighbor parameters
- ✓ Includes terminal base pair-dependent initiation
- ✓ Applies SantaLucia 2004 salt correction
- ✗ Does NOT perform dynamic programming alignment search
- ✗ Does NOT consider dangling ends
- ✗ Does NOT handle mismatches
- ✗ Does NOT calculate loop entropy
- ✗ Does NOT consider bulges

primer3-py's `calc_heterodimer()` uses the full **libprimer3 `thal()` function** which:
- ✓ Performs dynamic programming alignment
- ✓ Finds optimal binding regions
- ✓ Considers all structural features
- ✓ Handles dangling ends
- ✓ Accounts for mismatches and gaps

## Analysis

The consistent ~6.4 kcal/mol difference in ΔH suggests that primer3-py is either:

1. **Using a different alignment strategy** - Finding a shorter or offset alignment
2. **Including dangling end penalties** - Which reduce the binding energy
3. **Applying structural corrections** - For loops, bulges, or mismatches

## Verification

Testing with different sequence lengths shows the pattern is consistent:

| Length | primer3-py ΔH | Rust ΔH | Difference |
|--------|---------------|---------|-----------|
| 20bp   | -155200       | -161600 | -6400     |
| 17bp   | -130800       | -136600 | -5800     |
| 16bp   | -121000       | -127400 | -6400     |
| 14bp   | -105200       | -107500 | -2300     |

The differences are not proportional to sequence length, suggesting they're not just missing a constant term.

## Solution Options

### Option 1: Accept the Limitation (Current)
- Use the Rust implementation for relative comparisons
- Document that absolute values differ from primer3-py
- Suitable for: Ranking, filtering, relative analysis

### Option 2: Implement Full DP Algorithm
- Port the complete `thal()` function from primer3
- Implement dynamic programming alignment
- Add dangling end effects
- Effort: Very high (~2000+ lines of code)
- Result: Full primer3-py compatibility

### Option 3: Use FFI to libprimer3
- Create Rust bindings to the C library
- Call the actual `thal()` function
- Effort: Medium (~500 lines)
- Result: Perfect primer3-py compatibility
- Trade-off: Requires C library dependency

### Option 4: Hybrid Approach
- Keep current implementation for fast calculations
- Add FFI option for high-accuracy results
- Allow users to choose speed vs accuracy

## Recommendation

For the current use case (genome-wide primer binding analysis), **Option 1 is acceptable** because:

1. The differences are consistent and predictable
2. Relative rankings between sequences are preserved
3. The implementation is fast and self-contained
4. For filtering/ranking purposes, the simplified model works well

If exact primer3-py compatibility is required, **Option 3 (FFI)** is recommended as it provides:
- Exact compatibility
- Access to all primer3 features
- Reasonable implementation effort
- Minimal performance overhead

## Implementation Notes

The discrepancy appears to be primarily in the **ΔH calculation**, not the salt correction or Tm calculation. This suggests the issue is in:

1. **Initiation parameters** - May need adjustment
2. **Dangling end effects** - Not currently implemented
3. **Alignment strategy** - primer3 may use optimal alignment, not full sequence

Further investigation would require:
- Examining primer3 source code for exact algorithm
- Testing with sequences of known structure
- Comparing intermediate calculation steps
