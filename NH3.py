import numpy as np
import os
import psi4

# Set output file (Psi4 raw output)
psi4.core.set_output_file(os.path.join("Output", "NH3.dat"), False)

# Define NH3 molecule (approx geometry)
nh3 = psi4.geometry("""
                    0 1
                    N
                    H 1 1.01
                    H 1 1.01 2 107.0
                    H 1 1.01 2 107.0 3 120.0
                    units angstrom
                    symmetry c1
                    """)

# Set computational options
psi4.set_options({
    "e_convergence": 1e-8,
    "g_convergence": "gau_tight"
})

# Define multiple methods and basis sets for exploration
methods = ["b3lyp"]
basis_sets = ["sto-3g", "6-31G", "cc-pVDZ", "cc-pCVDZ", "aug-cc-pVDZ", "heavy-aug-cc-pVDZ", "cc-pVTZ"]

# Run calculations and print results
with open(os.path.join("Output", "NH3.txt"), "w") as f:
    f.write("Quantum Chemistry with Psi4 - NH3 Geometry & Optimization\n")
    for method in methods:
        for basis in basis_sets:
            psi4.set_options({"basis": basis})
            print(f"Running {method.upper()} with {basis} ...")

            # Geometry optimization
            opt_energy, opt_wfn = psi4.optimize(method, return_wfn=True)
            geom = opt_wfn.molecule().geometry().to_array()
            natoms = opt_wfn.molecule().natom()

            # Frequency
            freq_energy, freq_wfn = psi4.frequency(method, return_wfn=True)

            # Write
            f.write(f"\nMethod: {method.upper():<6}, Basis: {basis:<7}\n")
            f.write(f"  Optimized Energy   = {opt_energy:.8f} Eh\n")
            f.write("  Optimized Geometry (Å):\n")
            for i in range(natoms):
                atom = opt_wfn.molecule().symbol(i)
                x, y, z = geom[i]
                f.write(f"    {atom:<2} {x: .6f} {y: .6f} {z: .6f}\n")
            f.write(f"  Frequencies (cm^-1): {np.array(freq_wfn.frequencies())}\n")

print("End of NH3 geometry optimization and property calculations")