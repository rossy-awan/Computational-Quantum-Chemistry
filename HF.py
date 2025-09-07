import os
import psi4

# Set output file (ensure folder exists)
psi4.core.set_output_file(os.path.join("Output", "HF.dat"), False)

# Molecule definition
hf = psi4.geometry("""
                   0 1
                   H
                   F 1 0.916
                   units angstrom
                   symmetry c1
                   """)

# Set computational options
psi4.set_options({"scf_type": "pk", "e_convergence": 1e-8})

# Run calculations and print results
with open(os.path.join("Output", "HF.txt"), "w") as f:
    f.write("Quantum Chemistry with Psi4 - HF Non-Relativistic vs Relativistic\n")

    # Non-relativistic basis sets
    basis_nr = ["sto-3g", "6-31G", "cc-pVDZ", "aug-cc-pVDZ"]
    for basis in basis_nr:
        psi4.set_options({"basis": basis})
        print(f"Running SCF/{basis} (non-relativistic)...")
        energy, wfn = psi4.energy("scf", return_wfn=True)
        f.write(f"\nNon-Relativistic SCF/{basis} Energy = {energy:.8f} Eh\n")

        # Calculate properties: Mulliken, Lowdin charges, Dipole, Quadrupole
        psi4.oeprop(wfn, 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'DIPOLE')
        
        # Access and print properties
        f.write(f"  1. Mulliken Charges: {wfn.atomic_point_charges().np}\n")
        f.write(f"  2. Lowdin Charges: {wfn.variable('LOWDIN CHARGES')}\n")
        f.write(f"  3. Dipole Moment (a.u.): {wfn.variable('DIPOLE')}\n")

    # Relativistic (ECP)
    basis_rel = ["def2-SVP", "def2-TZVP", "def2-QZVP"]
    for basis in basis_rel:
        psi4.set_options({"basis": basis})
        print(f"Running SCF/{basis} (with ECP, relativistic)...")
        energy, wfn = psi4.energy("scf", return_wfn=True)
        f.write(f"\nRelativistic SCF/{basis} (ECP) Energy = {energy:.8f} Eh\n")

        # Calculate properties: Mulliken, Lowdin charges, Dipole, Quadrupole
        psi4.oeprop(wfn, 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'DIPOLE')
        
        # Access and print properties
        f.write(f"  1. Mulliken Charges: {wfn.atomic_point_charges().np}\n")
        f.write(f"  2. Lowdin Charges: {wfn.variable('LOWDIN CHARGES')}\n")
        f.write(f"  3. Dipole Moment (a.u.): {wfn.variable('DIPOLE')}\n")

    # Relativistic all-electron (DKH)
    basis_dkh = ["cc-pVDZ-DK", "cc-pVTZ-DK", "cc-pVQZ-DK"]
    for basis in basis_dkh:
        psi4.set_options({"basis": basis, "relativistic": "dkh"})
        print(f"Running SCF/{basis} (with DKH, relativistic)...")
        energy, wfn = psi4.energy("scf", return_wfn=True)
        f.write(f"\nRelativistic SCF/{basis} (DKH) Energy = {energy:.8f} Eh\n")

        psi4.oeprop(wfn, 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'DIPOLE')
        f.write(f"  1. Mulliken Charges: {wfn.atomic_point_charges().np}\n")
        f.write(f"  2. Lowdin Charges: {wfn.variable('LOWDIN CHARGES')}\n")
        f.write(f"  3. Dipole Moment (a.u.): {wfn.variable('DIPOLE')}\n")

    # Relativistic all-electron (X2C)
    basis_x2c = ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ"]
    for basis in basis_x2c:
        psi4.set_options({"basis": basis, "relativistic": "x2c"})
        print(f"Running SCF/{basis} (with X2C, relativistic)...")
        energy, wfn = psi4.energy("scf", return_wfn=True)
        f.write(f"\nRelativistic SCF/{basis} (X2C) Energy = {energy:.8f} Eh\n")

        psi4.oeprop(wfn, 'MULLIKEN_CHARGES', 'LOWDIN_CHARGES', 'DIPOLE')
        f.write(f"  1. Mulliken Charges: {wfn.atomic_point_charges().np}\n")
        f.write(f"  2. Lowdin Charges: {wfn.variable('LOWDIN CHARGES')}\n")
        f.write(f"  3. Dipole Moment (a.u.): {wfn.variable('DIPOLE')}\n")

print("End of HF non-relativistic vs relativistic calculations")