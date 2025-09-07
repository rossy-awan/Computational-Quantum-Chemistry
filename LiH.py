import matplotlib.pyplot as plt
import numpy as np
import os
import psi4
import time
from scipy.interpolate import CubicSpline

# Set Psi4 output file
psi4.core.set_output_file(os.path.join("Output", "LiH.dat"), False)

# Grid of bond lengths (20 points from 0.5 to 3.5 Å)
R_region1 = np.linspace(0.50, 1.45, 5, endpoint=False) # coarse left
R_region2 = np.linspace(1.45, 1.65, 8, endpoint=False) # dense around equilibrium
R_region3 = np.linspace(1.65, 3.50, 7)                 # coarse right
R_vals = np.unique(np.concatenate([R_region1, R_region2, R_region3]))

# Storage
results, summary = {}, []

# Run calculations
for basis in ["cc-pVDZ", "cc-pCVDZ", "aug-cc-pVDZ", "heavy-aug-cc-pVDZ", "cc-pVTZ"]:
    print(f"Running MP4 with {basis} ...")
    start_time = time.time()
    energies = []
    for R in R_vals:
        lih = psi4.geometry(f"""
                            0 1
                            Li
                            H 1 {R}
                            units angstrom
                            symmetry c1
                            """)
        psi4.set_options({"basis": basis, "e_convergence": 1e-8})
        energies.append(psi4.energy("mp4"))
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()
    energies = np.array(energies)
    elapsed = time.time() - start_time

    # Interpolation
    cs = CubicSpline(R_vals, energies)
    R_fine = np.linspace(0.5, 3.5, 1001)
    E_fine = cs(R_fine)

    # Find equilibrium bond length and energy
    idx_min = np.nanargmin(E_fine)
    R_min = R_fine[idx_min]
    E_min = E_fine[idx_min]

    # Store results
    results[basis] = (R_fine, E_fine)
    summary.append(("MP4", basis, E_min, R_min, elapsed))

# Write summary to file
with open(os.path.join("Output", "LiH.txt"), "w") as f:
    f.write("Quantum Chemistry with Psi4 - LiH Potential Energy Curve\n\n")
    f.write(f"{'Method':<10} {'Basis':<10} {'E_min (Eh)':<15} {'R_min (Å)':<15} {'Time (s)':<10}\n")
    for method, basis, E_min, R_min, elapsed in summary:
        f.write(f"{method:<10} {basis:<10} {E_min:<15.8f} {R_min:<15.4f} {elapsed:<10.2f}\n")
print("All calculations completed")

# Visualization
plt.rcParams.update({'font.family': 'Times New Roman', 'mathtext.fontset': 'cm', 'font.size': 18})
colors = ['#FA74A6', '#05BAB6', '#6FA4FF', '#6713F4', '#000000']
fig, ax = plt.subplots(figsize=(5.25, 5))
for i, (basis, (R_fine, E_fine)) in enumerate(results.items()):
    ax.plot(R_fine, E_fine, label=basis, color=colors[i])
ax.set_xlim(0.5, 3.5)
ax.set_ylim(-8.1, -7.5)
ax.set_xticks(np.linspace(0.5, 3.5, 7))
ax.set_yticks(np.linspace(-8.1, -7.5, 7))
ax.set_xlabel("Bond Length (Å)")
ax.set_ylabel("Energy (Hartree)")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join("Output", "LiH.svg"), dpi=300)
plt.show()