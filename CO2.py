import numpy as np
import os
import psi4

# Set output file (Psi4 raw output)
psi4.core.set_output_file(os.path.join("Output", "CO2.dat"), False)

# Define CO2 molecule (approx geometry)
co2 = psi4.geometry("""
                    0 1
                    O
                    C 1 1.16
                    O 2 1.16 1 180.0
                    units angstrom
                    symmetry c1
                    """)

# Set computational options
psi4.set_options({
    "e_convergence": 1e-8,
})

# Define multiple methods and basis sets for exploration
methods = ["scf"]
basis_sets = ["sto-3g", "6-31G", "cc-pVDZ", "cc-pCVDZ", "aug-cc-pVDZ", "heavy-aug-cc-pVDZ", "cc-pVTZ"]

# Run calculations and print results
with open(os.path.join("Output", "CO2.txt"), "w") as f:
    f.write("Quantum Chemistry with Psi4 - CO2 Frequency Analysis\n")
    for method in methods:
        for basis in basis_sets:
            psi4.set_options({"basis": basis})
            print(f"Running {method.upper()} with {basis} ...")

            # Frequency
            freq_energy, freq_wfn = psi4.frequency(method, return_wfn=True)
            
            # Write
            f.write(f"\nMethod: {method.upper():<6}, Basis: {basis:<7}\n")
            f.write(f"  Energy = {freq_energy:.8f} Eh\n")
            f.write(f"  Frequencies (cm^-1): {np.array(freq_wfn.frequencies())}\n")
            f.write(f"  Infrared intensity (km/mol): {np.array(freq_wfn.frequency_analysis['IR_intensity'].data)}\n")
            f.write(f"  Degeneracy: {np.array(freq_wfn.frequency_analysis['degeneracy'].data)}\n")
            f.write(f"  Char temp (K): {np.array(freq_wfn.frequency_analysis['theta_vib'].data)}\n")

print("End of CO2 calculations")