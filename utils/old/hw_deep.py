#!/usr/bin/env python
"""
SIESTA ToolBox - HAM_WANTIBEXOS Interface
Convert SIESTA Hamiltonian to Wantibexos code
"""

import os
import sys
import re
import argparse
import textwrap
from time import sleep
import numpy as np
import sisl
from typing import Tuple

VERSION = "2.1.0"
COLORS = {
    'reset': '\033[0m',
    'cyan': '\033[96m',
    'blue': '\033[94m',
    'green': '\033[92m',
    'yellow': '\033[93m',
    'red': '\033[91m',
    'bold': '\033[1m',
    'underline': '\033[4m'
}

def color_text(text: str, color: str) -> str:
    """Format text with ANSI color codes"""
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
    """Display STB-SUITE introduction banner"""
    os.system('cls' if os.name == 'nt' else 'clear')
    logo = color_text(r"""
.----------------.  .----------------.  .----------------.
| .--------------. || .--------------. || .--------------. |
| |    _______   | || |  _________   | || |   ______     | |
| |   /  ___  |  | || | |  _   _  |  | || |  |_   _ \    | |
| |  |  (__ \_|  | || | |_/ | | \_|  | || |    | |_) |   | |
| |   '.___`-.   | || |     | |      | || |    |  __'.   | |
| |  |`\____) |  | || |    _| |_     | || |   _| |__) |  | |
| |  |_______.'  | || |   |_____|    | || |  |_______/   | |
| |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'
""", 'cyan')

    description = [
        "Siesta ToolBox Suite",
        "SIESTA-to-Wantibexos Hamiltonian Interface",
        f"Version {VERSION} | University of Brasilia - 2025",
        ""
    ]
    
    print(logo)
    print("\n" + "="*60)
    for line in description:
        print(line.center(60))
        sleep(0.1)
    print("="*60 + "\n")

def parse_arguments() -> argparse.Namespace:
    """Handle command-line arguments with argparse"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''
        Convert SIESTA Hamiltonian to TB format for wantibexos
        Requires SIESTA calculation with SaveHS = .true. and SaveRho = .true.
        '''))
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input FDF file from SIESTA calculation')
    parser.add_argument('-o', '--output', default=None,
                        help='SIESTA output file for Fermi energy extraction')
    parser.add_argument('-f', '--fermi-level', type=float, default=None,
                        help='Manual Fermi energy setting (optional)')
    parser.add_argument('-q', '--quiet', action='store_false',
                        help='Suppress terminal animations and colors')
    
    return parser.parse_args()

def get_system_label(fdf_path: str) -> str:
    """Extract SystemLabel from SIESTA .fdf file with fallback to 'siesta'"""
    default_label = "siesta"
    try:
        with open(fdf_path, 'r') as f:
            for line in f:
                # Strip comments and whitespace
                line = re.sub(r'#.*', '', line).strip()
                if not line:
                    continue
                
                # Split into key-value pairs
                key_value = re.split(r'\s*=\s*|\s+', line, maxsplit=1)
                if key_value[0].lower() == 'systemlabel':
                    return key_value[1].split()[0] if len(key_value) > 1 else default_label
    except Exception as e:
        pass
    return default_label

def get_fermi_energy(output_file: str) -> float:
    """Extract Fermi energy from SIESTA output"""
    try:
        with open(output_file, 'r') as f:
            for line in reversed(f.readlines()):
                if 'Fermi =' in line:
                    return float(line.split()[-1])
    except (FileNotFoundError, IndexError) as e:
        raise RuntimeError(f"Fermi energy extraction failed: {str(e)}")

from typing import Tuple

def spin_mapping(spin_obj) -> Tuple[str, int]:
    """Map SIESTA spin types to Wantibexos conventions (code, multiplicity)."""
    spin_str = str(spin_obj).lower()
    mapping = {
        'unpolarized': ('NP', 1),
        'polarized': ('SP', 2),
        'non-colinear': ('SOC', 2),
        'spin-orbit': ('SOC', 2)
    }
    for key in mapping:
        if key in spin_str:
            return mapping[key]
    raise ValueError(f"Unsupported spin type in: '{spin_str}'")



def format_float(value: float) -> str:
    """Format floating point numbers for TB files"""
    return f"{value:12.6f}"


def write_basis(hamiltonian, spin_suffix: str, spin_factor: int) -> None:
    """Write basis set information to file"""
    filename = f"basis_set-{spin_suffix}.basis"
    with open(filename, 'w') as f:
        f.write("# bindex aspecie ax ay az l m spin\n")
        
        for spin in range(spin_factor):
            io = 0
            for ia, atom in enumerate(hamiltonian.geometry.atoms):
                xyz = hamiltonian.geometry.xyz[ia, :]
                for orbital in atom:
                    spin_val = 1 if spin == 0 else -1
                    if spin_suffix == 'SOC':
                        spin_val = 1 - 2*spin  # 1 for up, -1 for down
                    f.write(f"{io+1 + spin*hamiltonian.no:6d} {atom.tag:6} "
                            f"{format_float(xyz[0])} {format_float(xyz[1])} "
                            f"{format_float(xyz[2])} {orbital.l:2} "
                            f"{orbital.m:2} {spin_val:3}\n")
                    io += 1

def write_hamiltonian(hamiltonian, spin_suffix: str, fermi: float) -> None:
    """Write Hamiltonian data in TB format"""
    filename = f"tb-{spin_suffix}.ham"
    nbasis = hamiltonian.no
    ncell = np.prod(hamiltonian.nsc)

    with open(filename, 'w') as f:
        # Header section
        f.write(f"{spin_suffix}\n")
        f.write(f"{format_float(0.0)}\n")  # Scissors operator
        f.write(f"{format_float(fermi)}\n")
        
        # Cell vectors
        for vec in hamiltonian.cell:
            f.write(" ".join([format_float(v) for v in vec]) + "\n")
        
        # Dimensions
        f.write(f"{hamiltonian.no * (2 if spin_suffix != 'NP' else 1)}\n")
        f.write(f"{ncell}\n")
        f.write("# rcell_x rcell_y rcell_z i j ReH ImH S\n")
        
        # Hamiltonian data handling based on spin type
        if spin_suffix == 'NP':
            # Unpolarized case
            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        sc_index = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        H = hamiltonian[i, j + icell*nbasis][0]
                        S = hamiltonian.S[i, j + icell*nbasis]
                        
                        line = (f"{' '.join(format_float(v) for v in sc_index)} "
                                f"{i+1} {j+1} {format_float(H)} "
                                f"{format_float(0.0)} {format_float(S)}\n")
                        f.write(line)
        
        elif spin_suffix == 'SP':
            # Spin polarized case
            S = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            reH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            a = np.zeros((ncell+1, 2*nbasis+2))
            b = np.zeros((ncell+1, 2*nbasis+2))
            c = np.zeros((ncell+1, 2*nbasis+2))

            # Populate matrices
            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        # Up-up block
                        reH[icell+1, i+1, j+1] = hamiltonian[i, j + icell*nbasis][0]
                        S[icell+1, i+1, j+1] = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Down-down block
                        reH[icell+1, i+1+nbasis, j+1+nbasis] = hamiltonian[i, j + icell*nbasis][1]
                        S[icell+1, i+1+nbasis, j+1+nbasis] = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Cell coordinates
                        sc = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        a[icell+1, j+1] = sc[0]
                        b[icell+1, j+1] = sc[1]
                        c[icell+1, j+1] = sc[2]
                        a[icell+1, j+1+nbasis] = sc[0]
                        b[icell+1, j+1+nbasis] = sc[1]
                        c[icell+1, j+1+nbasis] = sc[2]

            # Write to file
            for icell in range(1, ncell+1):
                for i in range(1, 2*nbasis+1):
                    for j in range(1, 2*nbasis+1):
                        line = (f"{format_float(a[icell,j])} {format_float(b[icell,j])} "
                                f"{format_float(c[icell,j])} {i} {j} "
                                f"{format_float(reH[icell,i,j])} {format_float(0.0)} "
                                f"{format_float(S[icell,i,j])}\n")
                        f.write(line)
        
        elif spin_suffix == 'SOC':
            # Non-colinear and spin-orbit cases
            S = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            reH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            imH = np.zeros((ncell+1, 2*nbasis+2, 2*nbasis+2))
            a = np.zeros((ncell+1, 2*nbasis+2))
            b = np.zeros((ncell+1, 2*nbasis+2))
            c = np.zeros((ncell+1, 2*nbasis+2))

            for icell in range(ncell):
                for i in range(nbasis):
                    for j in range(nbasis):
                        sc = hamiltonian.geometry.o2sc(j + icell*nbasis)
                        idx = icell+1
                        
                        # Common cell coordinates
                        a[idx, j+1] = sc[0]
                        b[idx, j+1] = sc[1]
                        c[idx, j+1] = sc[2]
                        a[idx, j+1+nbasis] = sc[0]
                        b[idx, j+1+nbasis] = sc[1]
                        c[idx, j+1+nbasis] = sc[2]

                        # Matrix elements
                        H = hamiltonian[i, j + icell*nbasis]
                        S_val = hamiltonian.S[i, j + icell*nbasis]
                        
                        # Up-up block
                        reH[idx, i+1, j+1] = H[0]
                        S[idx, i+1, j+1] = S_val
                        
                        # Down-down block
                        reH[idx, i+1+nbasis, j+1+nbasis] = H[1]
                        S[idx, i+1+nbasis, j+1+nbasis] = S_val
                        
                        if len(H) > 2:  # Spin-orbit case
                            # Up-down block
                            reH[idx, i+1, j+1+nbasis] = H[2]
                            imH[idx, i+1, j+1+nbasis] = H[3]
                            
                            # Down-up block
                            reH[idx, i+1+nbasis, j+1] = H[6]
                            imH[idx, i+1+nbasis, j+1] = H[7]
                            
                            # Down-down imaginary
                            imH[idx, i+1+nbasis, j+1+nbasis] = H[5]
                        else:  # Non-colinear case
                            # Up-down block
                            reH[idx, i+1, j+1+nbasis] = H[2]
                            
                            # Down-up block
                            reH[idx, i+1+nbasis, j+1] = H[3]

            # Write to file
            for icell in range(1, ncell+1):
                for i in range(1, 2*nbasis+1):
                    for j in range(1, 2*nbasis+1):
                        line = (f"{format_float(a[icell,j])} {format_float(b[icell,j])} "
                                f"{format_float(c[icell,j])} {i} {j} "
                                f"{format_float(reH[icell,i,j])} "
                                f"{format_float(imH[icell,i,j])} "
                                f"{format_float(S[icell,i,j])}\n")
                        f.write(line)



def main():
    """Main workflow controller"""
    args = parse_arguments()
    
    print("\n" + color_text("SIESTA-WANTIBEXOS Interface:", 'bold'))
    print("-"*60)

    if not args.quiet:
        show_intro()
        print(color_text("[INFO] Initializing Hamiltonian processing...", 'yellow'))
    
    try:
        # Get system label and construct HSX filename
        system_label = get_system_label(args.input)

        # Load SIESTA Hamiltonian
        geom = sisl.get_sile("calc.fdf").read_geometry()
        print("[OK] Read geometry")
        ham = sisl.get_sile(args.input).read_hamiltonian(geometry=geom)
        print("[OK] Read Hamiltonian from SIESTA")
        # Determine spin configuration
        spin_type = str(ham.spin).split('.')[0].strip()
        spin_suffix, spin_factor = spin_mapping(spin_type)
        print("[OK] Determine the spin type calculation")
        
        # Get Fermi energy
        
        if args.fermi_level is not None:
            fermi = args.fermi_level
        elif args.output is not None:
            fermi = get_fermi_energy(args.output)
        else:
            fermi = 0.0

        print(f"[OK] Read the Fermi Energy ( {fermi} eV )")

    # Generate output files
        print("[INFO] Writing the basis... wait!")
        write_basis(ham, spin_suffix, spin_factor)
        print("[OK] Wrote the basis")
        print("[INFO] Writing the Hamiltonian... wait!")
        write_hamiltonian(ham, spin_suffix, fermi)
        print("[OK] Wrote the Hamiltonian!")
        
    except Exception as e:
        sys.exit(color_text(f"\nCritical error: {str(e)}", 'red'))
    if not args.quiet:
        print(color_text("[INFO] Successfully generated TB files:", 'green'))
        print(f"[OK] Basis set: basis_set-{spin_suffix}.basis", 'green')
        print(f"[OK] Hamiltonian: tb-{spin_suffix}.ham", 'green')
    
    print("[INFO] Complete job!")
    print("\n"+"-"*60)
    print(color_text("Saving Hamiltonians before they collapse their own wavefunctions.\n\n", 'bold'))

if __name__ == "__main__":
    main()
