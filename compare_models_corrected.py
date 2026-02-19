#!/usr/bin/env python3
"""
Corrected comparison script for Primer3 thermodynamic models.
The key issue: primer3-py calc_heterodimer calculates primer-target binding,
not primer self-binding (homodimer).
"""

import subprocess
import sys
import tempfile
import os

def test_with_primer3py(seq1, seq2, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0):
    """Test using primer3-py"""
    try:
        import primer3
        result = primer3.bindings.calc_heterodimer(
            seq1, seq2,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
            temp_c=temp_c
        )
        # primer3-py returns ΔG and ΔH in cal/mol, convert to kcal/mol for comparison
        return {
            'dg': result.dg / 1000.0,  # Convert cal/mol to kcal/mol
            'tm': result.tm,
            'dh': result.dh / 1000.0,  # Convert cal/mol to kcal/mol
            'ds': result.ds
        }
    except ImportError:
        print("ERROR: primer3-py not installed. Install with: pip3 install primer3-py")
        return None
    except Exception as e:
        print(f"ERROR in primer3-py: {e}")
        return None

def test_with_rust(binary, seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0):
    """Test using Rust binary"""
    try:
        # Create temporary FASTA files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">test\n{seq}\n")
            seq_file = f.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">pattern\n{seq}\n")
            pattern_file = f.name
        
        try:
            # Run Rust binary
            cmd = [
                f'./target/release/{binary}',
                '--file', seq_file,
                '--patterns', pattern_file,
                '--na', str(mv_conc),
                '--mg', str(dv_conc),
                '--dntp', str(dntp_conc),
                '--dnac', str(dna_conc),
                '--temp', str(temp_c),
                '--threshold', '0'  # Show all results
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode != 0:
                print(f"ERROR running {binary}: {result.stderr}")
                return None
            
            # Parse output
            lines = result.stdout.strip().split('\n')
            if lines and lines[0]:
                parts = lines[0].split('\t')
                if len(parts) >= 4:
                    return {
                        'dg': float(parts[2]),
                        'tm': float(parts[3]),
                        'dh': None,  # Not output by Rust version
                        'ds': None
                    }
            return None
        finally:
            os.unlink(seq_file)
            os.unlink(pattern_file)
    
    except Exception as e:
        print(f"ERROR running {binary}: {e}")
        return None

def compare_models(primer_seq, target_seq, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0):
    """Compare all three models"""
    print(f"\n{'='*80}")
    print(f"Primer:  {primer_seq}")
    print(f"Target:  {target_seq}")
    print(f"Parameters: Na+={mv_conc}mM, Mg2+={dv_conc}mM, dNTP={dntp_conc}mM, DNA={dna_conc}nM, T={temp_c}°C")
    print(f"{'='*80}")
    
    # Test primer3-py (heterodimer: primer vs target)
    print("\n[primer3-py calc_heterodimer(primer, target)]")
    py_result = test_with_primer3py(primer_seq, target_seq, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)
    if py_result:
        print(f"  ΔG:  {py_result['dg']:8.2f} kcal/mol")
        print(f"  Tm:  {py_result['tm']:8.2f} °C")
        if py_result['dh'] is not None:
            print(f"  ΔH:  {py_result['dh']:8.2f} kcal/mol")
            print(f"  ΔS:  {py_result['ds']:8.2f} cal/mol/K")
    else:
        py_result = None
    
    # Test main_simplesalt (primer binding to target)
    print("\n[main_simplesalt.rs (primer vs target)]")
    simple_result = test_with_rust('main_simplesalt', primer_seq, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)
    if simple_result:
        print(f"  ΔG:  {simple_result['dg']:8.2f} kcal/mol")
        print(f"  Tm:  {simple_result['tm']:8.2f} °C")
    else:
        simple_result = None
    
    # Test main_fullsalt (primer binding to target)
    print("\n[main_fullsalt.rs (primer vs target)]")
    full_result = test_with_rust('main_fullsalt', primer_seq, mv_conc, dv_conc, dntp_conc, dna_conc, temp_c)
    if full_result:
        print(f"  ΔG:  {full_result['dg']:8.2f} kcal/mol")
        print(f"  Tm:  {full_result['tm']:8.2f} °C")
    else:
        full_result = None
    
    # Compare
    print(f"\n{'─'*80}")
    print("Differences:")
    print(f"{'─'*80}")
    
    if simple_result and full_result:
        dg_diff = full_result['dg'] - simple_result['dg']
        tm_diff = full_result['tm'] - simple_result['tm']
        print(f"main_fullsalt vs main_simplesalt:")
        print(f"  ΔG difference: {dg_diff:+.2f} kcal/mol")
        print(f"  Tm difference: {tm_diff:+.2f} °C")
    
    if py_result and simple_result:
        dg_diff = py_result['dg'] - simple_result['dg']
        tm_diff = py_result['tm'] - simple_result['tm']
        print(f"\nprimer3-py vs main_simplesalt:")
        print(f"  ΔG difference: {dg_diff:+.2f} kcal/mol")
        print(f"  Tm difference: {tm_diff:+.2f} °C")
    
    if py_result and full_result:
        dg_diff = py_result['dg'] - full_result['dg']
        tm_diff = py_result['tm'] - full_result['tm']
        print(f"\nprimer3-py vs main_fullsalt:")
        print(f"  ΔG difference: {dg_diff:+.2f} kcal/mol")
        print(f"  Tm difference: {tm_diff:+.2f} °C")

def main():
    # Build Rust binaries first
    print("Building Rust binaries...")
    result = subprocess.run(['cargo', 'build', '--release', '--bins'], 
                          capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR building Rust binaries: {result.stderr}")
        sys.exit(1)
    
    # Test with default parameters
    print("\n" + "="*80)
    print("TESTING PRIMER-TARGET BINDING (Heterodimer)")
    print("="*80)
    
    test_cases = [
        ("ATGCGATCGATCG", "ATGCGATCGATCG", "Perfect match"),
        ("ATGCGATCGATCG", "CGCGATCGATCGA", "Shifted by 1"),
        ("ATGCGATCGATCG", "CGCGATCGATCGAT", "Longer target"),
        ("GCGCGCGCGC", "GCGCGCGCGC", "G-C rich perfect match"),
        ("ATATATATATAT", "ATATATATATAT", "A-T rich perfect match"),
    ]
    
    for primer, target, description in test_cases:
        print(f"\n[{description}]")
        compare_models(primer, target, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0)
    
    # Test with different salt concentrations
    print("\n" + "="*80)
    print("TESTING WITH DIFFERENT SALT CONCENTRATIONS")
    print("="*80)
    
    primer = "ATGCGATCGATCG"
    target = "ATGCGATCGATCG"
    
    print("\n[High Mg2+ (5.0 mM)]")
    compare_models(primer, target, mv_conc=50.0, dv_conc=5.0, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0)
    
    print("\n[Low Na+ (10.0 mM)]")
    compare_models(primer, target, mv_conc=10.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=37.0)
    
    print("\n[Different Temperature (25°C)]")
    compare_models(primer, target, mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6, dna_conc=50.0, temp_c=25.0)

if __name__ == '__main__':
    main()
