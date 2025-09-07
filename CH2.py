import os
import psi4

# Set output file (ensure folder exists)
psi4.core.set_output_file(os.path.join("Output", "CH2.dat"), False)

# Define molecule (CH2 with simple geometry)
R, A = 1.075, 133.93
ch2 = psi4.geometry("""
                    0 3
                    C
                    H 1 {0}
                    H 1 {0} 2 {1}
                    """.format(R, A)
                    )

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
with open(os.path.join("Output", "CH2.txt"), "w") as f:
    f.write("Quantum Chemistry with Psi4 - CH2 Properties\n")
    for method in methods:
        for basis in basis_sets:
            psi4.set_options({"basis": basis})
            print(f"Running {method.upper()} with {basis} ...")

            # Run calculation and get wavefunction
            energy, wfn = psi4.energy(method, return_wfn=True)
            f.write(f"\nMethod: {method.upper():<5}, Basis: {basis:<7} --> Energy = {energy:.8f} Eh\n")

            # Calculate properties: Mulliken, Lowdin charges, Dipole, Quadrupole
            psi4.oeprop(wfn, 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'DIPOLE')
            
            # Access and print properties
            f.write(f"  1. Mulliken Charges: {wfn.atomic_point_charges().np}\n")
            f.write(f"  2. Lowdin Charges: {wfn.variable('LOWDIN CHARGES')}\n")
            f.write(f"  3. Dipole Moment (a.u.): {wfn.variable('DIPOLE')}\n")
print("End of CH2 property calculations")