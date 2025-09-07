import numpy as np
import os
import psi4

# Set output file (ensure folder exists)
psi4.core.set_output_file(os.path.join("Output", "CO.dat"), False)

# Define molecule
co = psi4.geometry("""
                   0 1
                   C
                   O 1 1.128
                   units angstrom
                   symmetry c1
                   """)

# Set computational options
psi4.set_options({
    "scf_type": "pk",
    "e_convergence": 1e-8,
    'reference': 'uhf'
})

# Define multiple methods and basis sets for demonstration
methods = ["scf", "mp2", "b3lyp"]
basis_sets = ["sto-3g", "6-31G", "cc-pVDZ", "aug-cc-pVDZ"]

# Run calculations and print results
with open(os.path.join("Output", "CO.txt"), "w") as f:
    f.write("Quantum Chemistry with Psi4 - CO Properties\n")
    for method in methods:
        for basis in basis_sets:
            psi4.set_options({"basis": basis})
            print(f"Running {method.upper()} with {basis} ...")

            # Run calculation and get wavefunction
            energy, wfn = psi4.energy(method, return_wfn=True)
            f.write(f"\nMethod: {method.upper():<5}, Basis: {basis:<7} --> Energy = {energy:.8f} Eh\n")

            # Calculate properties
            psi4.oeprop(wfn, 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'DIPOLE', 'QUADRUPOLE',
                        'MAYER_INDICES', 'WIBERG_LOWDIN_INDICES', 'ESP_AT_NUCLEI')
            
            # Access and print properties
            f.write(f"  1. Mulliken Charges: {wfn.atomic_point_charges().np}\n")
            f.write(f"  2. Lowdin Charges: {wfn.variable('LOWDIN CHARGES')}\n")
            f.write(f"  3. Dipole Moment (a.u.): {wfn.variable('DIPOLE')}\n")
            f.write(f"  4. Quadrupole Moment (a.u.): {wfn.variable('QUADRUPOLE')}\n")
            f.write(f"  5. Mayer Indices: {np.array(wfn.variable('MAYER INDICES'))}\n")
            f.write(f"  6. Wiberg Lowdin Indices: {np.array(wfn.variable('WIBERG LOWDIN INDICES'))}\n")
            f.write(f"  7. Electrostatic Potential at Nuclei 1: {np.array(wfn.variable('ESP AT CENTER 1'))}\n")
            f.write(f"  8. Electrostatic Potential at Nuclei 2: {np.array(wfn.variable('ESP AT CENTER 2'))}\n")

print("End of CO property calculations")