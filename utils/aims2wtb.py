

#Version
__version__ = "1.0.0"

# Import Library

import numpy as np
import os
import re
from glob import glob
from itertools import product
from sklearn.linear_model import LinearRegression
from tabulate import tabulate
import argparse

def print_header():
    print(fr"""
====================================================================
      WANTIBEXOS TOOLS  |  MODULE: fhi2wtb
--------------------------------------------------------------------
   ↪ Version     : {__version__}
   ↪ Purpose     : Convert FHI-aims outputs into real-space H(R)
                   and prepare data for tight-binding models.
--------------------------------------------------------------------
   ->  Ready to Fourier your way through the Brillouin zone?
   ->  Let’s discretize, transform, and export like pros.
====================================================================
""")

def print_footer():
    print(r"""
====================================================================
   ✅  fhi2wtb completed all tasks successfully
--------------------------------------------------------------------
   ->   Output file: fhi_tb.ham
   ->   Use this in Wantibexos code.
--------------------------------------------------------------------
       Thanks for using Wantibexos Tools!
       Stay symmetric. Stay localized.
====================================================================
""")

def print_footer_kp():
    print(r"""
====================================================================
--------------------------------------------------------------------
   ->   Output file: k-points_for_FHI.dat
   ->   Use k-points in control.in in FHI code.
   ->   Active the flags in FHI:
            output hamiltonian_matrix
            output overlap_matrix
--------------------------------------------------------------------
       Thanks for using Wantibexos Tools!
       Stay symmetric. Stay localized.
====================================================================
""")




def write_fhi_static_kpoints(kpoints_reduced,mesh_k,filename="k-points_for_FHI.dat"):
    with open(filename, 'w') as f:
        f.write("# Monkhorst-Pack grid as individual k-points for FHI-aims\n")
        f.write("# Each path is a single k-point (start == end), no interpolation\n\n")
        f.write(f"Mesh_k: {mesh_k[0]}  {mesh_k[1]}  {mesh_k[2]}\n\n")
        f.write("Kpoint_fractional:\n")
        for k in kpoints_reduced:
            f.write(
                f"output band  {k[0]: .8f} {k[1]: .8f} {k[2]: .8f}   {k[0]: .8f} {k[1]: .8f} {k[2]: .8f}   2\n")
        f.write("End_Kpoint_fractional\n\n")
        f.write("!!!!!!!! Remember!!!!!!!!!\n")
        f.write("1- Copy and paste the points in control.in in FHI-aims\n")        
        f.write("2- Active the flags: \n    output hamiltonian_matrix\n    output overlap_matrix\n")
        f.write("!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
        f.write("# End of file")



def read_lattice_vectors_from_geometry(filename="geometry.in"):
    """
    Read lattice vectors a1, a2, a3 from FHI-aims geometry.in file.

    Returns:
        a1, a2, a3: numpy arrays (shape = (3,))
    """
    lattice_vectors = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith('lattice_vector'):
                parts = line.strip().split()
                if len(parts) != 4:
                    raise ValueError(f"Invalid lattice_vector line: {line}")
                vec = [float(x) for x in parts[1:]]
                lattice_vectors.append(vec)
    if len(lattice_vectors) != 3:
        raise ValueError("Expected exactly 3 lattice_vector entries in geometry.in")
    a1, a2, a3 = map(np.array, lattice_vectors)
    return a1, a2, a3

def generate_monkhorst_pack_grid(nk,shift=False):
    """
    Create the Monkhorst-Pack a grid for k.
    """
    #verify if nk is odd
    k_mesh=[]
    icont=1
    for n in nk:
        if n % 2 == 0:
            k_mesh.append(n+1)
            print(f"\n[INFO] Adjust N{icont} to correct grid in Fourier Transform.")
            icont=icont+1
        else:
            k_mesh.append(n)
    if icont!= 1:
        print(f"\n[INFO] New grid: {k_mesh}\n")
    mesh = []
    for n in k_mesh:
        if shift==False:
            ks = (np.arange(n) + 0.5) / n - 0.5
        else:
            ks = np.arange(n) / n - 0.5
        mesh.append(ks)
    grid = np.meshgrid(*mesh, indexing='ij')
    h = np.stack([g.ravel() for g in grid], axis=-1)
    return h,k_mesh

def generate_R_vectors(nk, a1, a2, a3):
    """
    Generate the symmetric vectors R reduced Li and cartesianas R = L1*a1 + L2*a2 + L3*a3.
    """
    freq_lists = [np.fft.fftfreq(n, d=1) * n for n in nk]
    freq_lists = [freq.astype(int) for freq in freq_lists]
    R_red = list(product(*freq_lists))
    R_cart = [L[0]*a1 + L[1]*a2 + L[2]*a3 for L in R_red]
    return R_red, R_cart

def check_R_symmetry(R_list,k_list):
    R_set = {tuple(np.round(R, 8)) for R in R_list}
    for R in R_list:
        neg_R = tuple(np.round(-R, 8))
        if neg_R not in R_set:
            print(f"\n[FAIL] R = {R} does not have a corresponding -R!")
            return False
    print("\n[OK] All R vectors are symmetric with respect to the origin.")
    
    k_set = {tuple(np.round(k, 8)) for k in k_list}
    for k in k_list:
        neg_k = tuple(np.round(-k, 8))
        if neg_k not in k_set:
            print(f"\n[FAIL] k = {k} does not have a corresponding -k!")
            return False
    print("\n[OK] All R vectors are symmetric with respect to the origin.")
    return

def print_correspondence(kpoints_red, R_red, R_cart):
    """
    Sort and Print the R and correspondente k point
    """
    R_norm = [np.linalg.norm(R) for R in R_cart]
    data = list(zip(R_norm, R_red, kpoints_red))
    data.sort()
    print("\nIndex |      R (int)         | |R|     ||     k (reduced)")
    print("-------------------------------------------------------------")
    for i, (rnorm, R, k) in enumerate(data):
        R_str = f"({R[0]:2d}, {R[1]:2d}, {R[2]:2d})"
        k_str = f"({k[0]: .6f}, {k[1]: .6f}, {k[2]: .6f})"
        print(f"{i:5d} | {R_str} | {rnorm:6.3f} || {k_str}")
    print(f"\nTotal points: {len(kpoints_red)}")
    return

def positons_define(directory,k_grid,shift=False):
    """
    calculate the Li and R from discretization of the k grid N1xN2xN3: 
    Li=Ni/2 --> even  Li=(Ni-1)/2 --> impar and R=L1xa1+L2xa2+L3xa3
    Returns:
       L1,L2,L3: numpy arrays (shape=3) and R1,R2, R3: numpy array (shape=3)
    """
    print(f"\n[INFO] Generating Monkhorst-Pack and R-points for nk = {k_grid}")
    
    a1,a2,a3=read_lattice_vectors_from_geometry(filename=directory+"/geometry.in")
    # Gerar malha k em coordenadas reduzidas
    kpoints_red,k_grid = generate_monkhorst_pack_grid(k_grid,shift)
    # Gerar malha R correspondente
    R_red, R_cart = generate_R_vectors(k_grid, a1, a2, a3)

    # Verificar correspondência 1-1
    if len(kpoints_red) != len(R_red):
        raise ValueError(f"[FAIL] k and R numbers don't have correspondence dim(k)={len(kpoints_red)} and dim(R)={len(R_red)}!")
    # Imprimir
    print_correspondence(kpoints_red, R_red, R_cart)
    check_R_symmetry(np.array(R_red),np.array(kpoints_red))
    return R_cart,R_red,k_grid




def reciprocal_lattice_from_direct():
    """
    Return the reciprocal lattice from direct lattice
    """
    a1,a2,a3=read_lattice_vectors_from_geometry()
    volume = np.dot(a1, np.cross(a2, a3))
    b1 = 2 * np.pi * np.cross(a2, a3) / volume
    b2 = 2 * np.pi * np.cross(a3, a1) / volume
    b3 = 2 * np.pi * np.cross(a1, a2) / volume
    return b1, b2, b3

def kpoints_reduced_to_cartesian_from_lattice(kpoints_reduced):
    """
    Convert k-points reduced in lattice cartesian
    """
    b1, b2, b3 = reciprocal_lattice_from_direct()
    B = np.vstack([b1, b2, b3])
    kpoints_reduced = np.array(kpoints_reduced)
    ckp= kpoints_reduced @ B
    return tuple(round(x, 8) for x in ckp)

def parse_k_point_from_header(line):
    """
    Extract k-point from the header comment line.
    Example line:
    # k-point 1 : at relative reciprocal-space coordinates : 0.00000000 0.00000000 0.00000000
    """
    match = re.search(r"coordinates\s*:\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)", line)
    if match:
        kx, ky, kz = map(float, match.groups())
        return (kx, ky, kz)
    else:
        raise ValueError("Failed to extract k-point from line: " + line)

def read_complex_hk_matrix(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Find k-point from first comment line
    for line in lines:
        if line.strip().startswith("#") and "coordinates" in line:
            k_point = parse_k_point_from_header(line)
            break
    else:
        raise ValueError(f"No k-point line found in {filename}")
    # Locate the start of the matrix
    try:
        spin_start = lines.index("# spin channel        1\n")
    except ValueError:
        raise ValueError(f"'spin channel 1' section not found in {filename}")
    # Read matrix data after spin_start
    data = []
    for line in lines[spin_start+1:]:
        if line.strip() == "":
            continue
        values = list(map(float, line.strip().split()))
        if len(values) % 2 != 0:
            raise ValueError(f"Unexpected number of columns in line: {line}")
        # Combine real and imaginary parts
        row = [complex(values[i], values[i+1]) for i in range(0, len(values), 2)]
        data.append(row)
    matrix = np.array(data, dtype=complex)
    return k_point, matrix

def read_complex_ok_matrix(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Find k-point from first comment line
    for line in lines:
        if line.strip().startswith("#") and "coordinates" in line:
            k_point = parse_k_point_from_header(line)
            break
    else:
        raise ValueError(f"No k-point line found in {filename}")
    # Locate the start of the matrix
    try:
        spin_start = lines.index("#\n")
    except ValueError:
        raise ValueError(f"'spin channel 1' section not found in {filename}")
    # Read matrix data after spin_start
    data = []
    for line in lines[spin_start+1:]:
        if line.strip() == "":
            continue
        values = list(map(float, line.strip().split()))
        if len(values) % 2 != 0:
            raise ValueError(f"Unexpected number of columns in line: {line}")
        # Combine real and imaginary parts
        row = [complex(values[i], values[i+1]) for i in range(0, len(values), 2)]
        data.append(row)
    matrix = np.array(data, dtype=complex)
    return k_point, matrix

def read_all_hk_matrices(name,directory="."):
    """
    Reads all *_matrix.* files in the directory and returns a dict:
    { (kx, ky, kz): matrix }
    """
    k_matrices = {}
    for filepath in glob(os.path.join(directory, name)):
        try:
            k, mat = read_complex_hk_matrix(filepath)
            k_matrices[k] = mat
            print(f"[OK] Read matrix for k = {k} from {os.path.basename(filepath)}")
        except Exception as e:
            print(f"[ERROR] Skipping {filepath}: {e}")
    return k_matrices

def read_all_ok_matrices(name,directory="."):
    """
    Reads all *_matrix.* files in the directory and returns a dict:
    { (kx, ky, kz): matrix }
    """
    k_matrices = {}
    for filepath in glob(os.path.join(directory, name)):
        try:
            k, mat = read_complex_ok_matrix(filepath)
            k_matrices[k] = mat
            print(f"[OK] Read matrix for k = {k} from {os.path.basename(filepath)}")
        except Exception as e:
            print(f"[ERROR] Skipping {filepath}: {e}")
    return k_matrices



def fourier_transform_hamiltonian_fast(tp,Hk_dict, n_points,directory,force_h=False, imag_threshold=None,simmetrization=False):
    """
    Fast version of Fourier transform from H(k) to H(R) using numpy vectorization.
    Args:
        Hk_dict: dict where keys are (kx, ky, kz) and values are complex matrices H(k)
        n_points: list of R integer vectors (n1, n2, n3)
        imag_threshold: threshold to remove small imaginary parts
    Returns:
        HR_dict: dict where keys are R integer tuples and values are complex matrices H(R)
    """
    Hartree= 27.211386245988 #Convert hartree to eV.
    k_list = np.array(list(Hk_dict.keys()))   # (Nk, 3)
    Hk_list = np.array(list(Hk_dict.values())) # (Nk, dim, dim)
    n_list = np.array(n_points)                # (Nr, 3)
    Nk = len(k_list)
    Nr = len(n_list)
    dim = Hk_list.shape[1]
    # Phase Matrix: (Nr, Nk)
    phases = np.exp(-2j * np.pi * np.dot(n_list, k_list.T))  # (Nr, Nk)
    # Apply Fourier Transform: (Nr, dim, dim)
    HR_array = np.einsum('rk,kij->rij', phases, Hk_list) / Nk
    # Force Hermicity
    if force_h==True:
        HR_array = 0.5 * (HR_array + HR_array.conj().transpose(0,2,1))  # (rij) -> (rij)*
        print(f"\n[OK] Forced the Hermitian {'Hamiltonian' if tp=='ham' else 'Overlap'} matrix")
    # Clear numerical residus
    if imag_threshold!=None:
        HR_array.imag[np.abs(HR_array.imag) < imag_threshold] = 0.0
        print(f"\n[OK] Removed imaginary values < {imag_threshold} in  {'Hamiltonian' if tp=='ham' else 'Overlap'} matrix")
    # Calculate the lattice vectors    
    a1,a2,a3=read_lattice_vectors_from_geometry(filename=directory+"/geometry.in")
    list_r=[]
    for n in n_list:
        list_r.append(np.dot(n,[a1,a2,a3]))
    # Convert to dictionary
    if tp=='ham':
        HR_dict = {tuple(n): HR*Hartree for n, HR in zip(list_r, HR_array)}
    if tp=='ove':
        HR_dict = {tuple(n): HR for n, HR in zip(list_r, HR_array)}

    # Simetrization
    if simmetrization==True:
        for R in HR_dict.keys():
            HR = HR_dict[R]
            HR_minus = HR_dict.get(tuple(-np.array(R)), np.zeros_like(HR))
            HR_sym = 0.5 * (HR + HR_minus.conj().T)
            HR_dict[R] = HR_sym
        print(f"\n[OK] Symmetrize the {'Hamiltonian' if tp=='ham' else 'Overlap'} matrix")

    return HR_dict


def fourier_transform_hamiltonian(tp,Hk_dict,n_points,directory,force_h=False, imag_threshold=None,simmetrization=False):
    """
    Performs a discrete Fourier transform from H(k) to H(R),
    and sets small imaginary parts (abs < imag_threshold) to zero.

    Args:
        Hk_dict: dict where keys are (kx, ky, kz) and values are complex matrices H(k)
        rposition: list of R vectors (tuples or arrays)
        imag_threshold: threshold below which imaginary parts are discarded

    Returns:
        HR_dict: dict where keys are R (as tuple) and values are complex matrices H(R)
    """
    Hartree= 27.211386245988 #Convert hartree to eV.
    HR_dict = {}
    k_points = np.array(list(Hk_dict.keys()))
    Hk_matrices = list(Hk_dict.values())
    for R in n_points:
        R = np.array(R)
        HR = np.zeros_like(Hk_matrices[0], dtype=complex)
        for k, Hk in zip(k_points, Hk_matrices):
            phase = np.exp(-2j * np.pi * np.dot(k, R))
            HR += phase * Hk
        HR /= len(k_points)
        #force Hermicity
        if force_h==True:
            HR = 0.5 * (HR + HR.conj().T)  # (rij) -> (rij)*
            print(f"\n[OK] Forced the Hermitian {'Hamiltonian' if tp=='ham' else 'Overlap'} matrix")
        # Clear numerical residus
        if imag_threshold!=None:
            HR.imag[np.abs(HR.imag) < imag_threshold] = 0.0
            print(f"\n[OK] Removed imaginary values < {imag_threshold} in  {'Hamiltonian' if tp=='ham' else 'Overlap'} matrix")
        # Calculate the lattice vectors and cartesian R points  
        a1,a2,a3=read_lattice_vectors_from_geometry(filename=directory+"/geometry.in")
        cart=np.dot(R,[a1,a2,a3])
        if tp=='ham':
            HR_dict[tuple(cart)] = HR*Hartree
        if tp=='ove':
            HR_dict[tuple(cart)] = HR
    # simetrization
    if simmetrization==True:
        for R in HR_dict.keys():
            HR = HR_dict[R]
            HR_minus = HR_dict.get(tuple(-np.array(R)), np.zeros_like(HR))
            HR_sym = 0.5 * (HR + HR_minus.conj().T)
            HR_dict[R] = HR_sym
        print(f"\n[OK] Symmetrize the {'Hamiltonian' if tp=='ham' else 'Overlap'} matrix")
    return HR_dict

def check_if_hamiltonians_are_real(HR_dict, tol=1e-10):
    """
    Verrify if all matrix are reals (do not have imaginary part)
    """
    all_real = True
    for R, H in HR_dict.items():
        if np.any(np.abs(H.imag) > tol):
            print(f"\n[WARNING] H(R) at R = {R} has non-negligible imaginary components.")
            all_real = False
    if all_real:
        print("\n[OK] All matrices are real (within tolerance).")
    return all_real

def write_wantibexos(h_real,s_real,scisor,directory):
    lattice=read_lattice_vectors_from_geometry(directory+"/geometry.in")
    outfile = open("fhi_tb.ham", "w")
    outfile.write(f"NP\n")
    outfile.write(f"{scisor}\n")
    outfile.write(f"0.000000\n")
    for i in range(3):
        outfile.write(f"{lattice[i][0]:.8f} {lattice[i][1]:.8f} {lattice[i][2]:.8f}\n")
    fk= next(iter(h_real))
    outfile.write(f"{len(h_real[fk][0])}\n")
    outfile.write(f"{len(h_real)}\n")
    outfile.write(f"#rcell x   rcell y   rcell z   i   j   ReH   ImH   S\n")
    for k in h_real.keys():
        h_matrix=h_real[k]
        s_matrix=s_real[k]
        for (i, j), el in np.ndenumerate(h_matrix):
            s_el=s_matrix[i,j]
            outfile.write(f"{k[0]:.8f}  {k[1]:.8f}  {k[2]:.8f} {i+1}  {j+1}  {el.real:.8f}  {el.imag:.8f}  {s_el.real:.8f}\n")
    outfile.close()
    return

def is_hermitian(H, tol=1e-10):
    """
    Verifica se uma matriz é hermitiana: H == H†

    Args:
        H: matriz complexa
        tol: tolerância para comparar (default: 1e-10)

    Returns:
        True se H for hermitiana, False caso contrário
    """
    return np.allclose(H, H.conj().T, atol=tol)

def check_number_matrix_and_grid(matrix,matrix2,nr):
    if len(matrix)==len(nr) and len(matrix2)==len(nr):
        print("\n[OK] Check: The number of matrix is the same of k-points and R-points")
    else:
        print("\n[FAIL] The number of matrix is different of k-points and R-points!")
        raise ValueError("Verify the grid of k-points")


def analyze_matrices(matrix_dict):
    """
    Analyze a dictionary of matrices: check hermiticity, nullity, and linear trend (R^2)
    between real and imaginary parts. Also compute distance from origin and sort accordingly.

    Args:
        matrix_dict (dict): Keys are 3D positions (tuples), values are 2D numpy arrays (matrices).

    Prints:
        A formatted table with position, distance, Hermitian check, null check, 
        and R^2 score of linear regression (real vs. imaginary).
    """
    results = []

    for pos, mat in matrix_dict.items():
        mat = np.array(mat)
        is_hermitian = np.allclose(mat, mat.conj().T)
        is_null = np.allclose(mat, 0)

        real = mat.real.flatten().reshape(-1, 1)
        imag = mat.imag.flatten()

        # R² from linear regression: imag ~ real
        if np.all(imag == 0) or np.all(real == 0):
            r_squared = np.nan
        else:
            model = LinearRegression().fit(real, imag)
            r_squared = model.score(real, imag)

        # Distance from origin
        distance = np.linalg.norm(pos)

        results.append([
            str(pos),
            distance,
            "yes" if is_hermitian else "no",
            "yes" if is_null else "no",
            f"{r_squared:.3f}" if not np.isnan(r_squared) else "n/a"
        ])

    # Sort by distance
    results.sort(key=lambda x: x[1])

    # Print table
    print(tabulate(
        results,
        headers=["Position", "Distance", "Hermitian", "Null", "Linear Trend (R²)"],
        tablefmt="github",
        floatfmt=".3f"
    ))
    return


def main_generator_kpoints(nk,shift=False):
    print_header()
    print(f"\n[INFO] Monkhorst-Pack mesh: {nk[0]} x {nk[1]} x {nk[2]}")
    print(f"\n[INFO] Shifted grid: {'Yes' if shift else 'No (centered at 0)'}\n")
    print("[INFO] Generating Monkhorst-Pack k-points...")
    kpoints_reduced,mesh_k = generate_monkhorst_pack_grid(nk,shift)
    print(f"[OK]   Total k-points generated: {len(kpoints_reduced)}\n")
    print(f"[INFO] Writing k-points to FHI-aims file: k-points_for_FHI.dat")
    write_fhi_static_kpoints(kpoints_reduced,mesh_k)
    print("\n[OK]   Output written successfully.")
    print("\n[DONE] All tasks completed.")
    print_footer_kp()
    return

def main_generator_hamiltonian(nk,directory,shift=False,method='vec',force_h=False, imag_threshold=None,simmetrization=False):
    print_header()
    # Determine R points
    print(f"[INFO] Monkhorst-Pack mesh: {nk[0]} x {nk[1]} x {nk[2]}")
    print(f"[INFO] Shifted grid: {'Yes' if shift else 'No (centered at 0)'}\n")
    rposition,n_points,k_grid=positons_define(directory,nk,shift)
    #extract the hamiltonian H(k) and k points
    Hk_dict = read_all_hk_matrices("KS_hamiltonian_matrix.*.kpt_1.out",directory)
    print(f"[INFO] Total Hamiltian matrix and k-points read: {len(Hk_dict)}")
    #extract the hamiltonian S(k) and k points
    Sk_dict = read_all_ok_matrices("KS_overlap_matrix.*.kpt_1.out",directory)
    print(f"[INFO] Total Overlap matrix and k-points read: {len(Sk_dict)}")
    #Check the numer of matrix is the same of kpoints and R.
    check_number_matrix_and_grid(Hk_dict,Sk_dict,n_points)
    #Make the Fourier Transform
    if method=="vec":
        HR_dict = fourier_transform_hamiltonian_fast('ham',Hk_dict, n_points,directory,force_h, imag_threshold,simmetrization)
        SR_dict = fourier_transform_hamiltonian_fast('ove',Sk_dict, n_points,directory,force_h, imag_threshold,simmetrization)
    if method=="std":
        HR_dict = fourier_transform_hamiltonian('ham',Hk_dict, n_points,directory,force_h, imag_threshold,simmetrization)
        SR_dict = fourier_transform_hamiltonian('ove',Sk_dict, n_points,directory,force_h, imag_threshold,simmetrization)
    print("\n[INFO] Check if the Hamiltonian matrix are real...")
    check=check_if_hamiltonians_are_real(HR_dict)
    print("\n[INFO] Check if the Overlap matrix are real...")
    check=check_if_hamiltonians_are_real(SR_dict)
    print("\n[INFO] Hamiltonian Matrix in Real Space:")
    check=analyze_matrices(HR_dict)
    print("\n[INFO] Overlap Matrix in Real Space:")
    check=analyze_matrices(SR_dict)
    print("\n[INFO] Write in Wantibexos format...")
    write_wantibexos(HR_dict,SR_dict,scisor,directory)
    print("\n[INFO] Complete. Created a file fhi_tb.ham !")
    print("\n[DONE] All tasks completed.\n")
    print_footer()
    return

def main():
    parser = argparse.ArgumentParser(description="Wantibexos Toolbox - FHI-aims k-points & Wantibexos export Hamiltonian. Choose between --kpoints or --wtb")

    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--kpoints', action='store_true', help='[req] Run k-points generator (Monkhorst-Pack)')
    mode_group.add_argument('--wtb', action='store_true', help='[req] Run Fourier transform to generate H(R) and S(R)')

    parser.add_argument('--nk', type=int, nargs=3, metavar=('N1', 'N2', 'N3'),
                        help='[req] Monkhorst-Pack mesh size (e.g. --nk 9 9 3)', required=True)
    parser.add_argument('--directory', type=str, default='.', help='[req-wtb]Directory with FHI-aims outputs')
    parser.add_argument('--shift', action='store_true', help='[opt] Use shifted Monkhorst-Pack mesh. Default no shift')
    parser.add_argument('--method', choices=['vec', 'std'], default='vec',
                        help='[opt-wtb] Fourier transform method: "vec" (vectorized) or "std" (standard). Default is vec.')
    parser.add_argument('--scissor', type=str, default="0.00000", help='[opt-wtb] Scissor operator value (string format).Default 0.00.')
    parser.add_argument('--force-hermitian',  action='store_true', help='[opt-wtb] Force Hermitian H(R) matrices. Default is not activate. Warning: Can be generate a wrong results due to force a symmetry!!!)')
    parser.add_argument('--imag-threshold',type=float, help='[opt-wtb] Threshold to remove small imaginary components (e.g. 1e-10). Default is none')
    parser.add_argument('--symmetrize',  action='store_true', help='[opt-wtb] Symmetrize matrices between R and -R. Default is not activate.')

    args = parser.parse_args()

    if args.kpoints:
        main_generator_kpoints(args.nk,args.shift)

    elif args.wtb:
        # Set global variables if needed by subroutines
        global scisor
        scisor = args.scissor
        global shift
        shift = args.shift

        main_generator_hamiltonian(
            nk=args.nk,
            directory=args.directory,
            method=args.method,
            force_h=args.force_hermitian,
            imag_threshold=args.imag_threshold,
            simmetrization=args.symmetrize)

if __name__ == "__main__":
    main()

