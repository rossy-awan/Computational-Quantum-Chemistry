# Quantum Chemistry with Psi4

This repository is designed for **self-learners** who want to explore **computational quantum chemistry** using [Psi4](https://www.psicode.org). It provides Python scripts to calculate energies, electronic properties, geometry optimizations, vibrational frequencies, and to analyze the effect of different methods and basis sets.

## Getting Started
Explore basic computations for small molecules:
- H₂O, O₂, CH₂
- Compute total energies and basic electronic properties
- Use a variety of methods and basis sets

## Electronic Properties
Detailed electronic property calculations using CO:
- Mulliken & Lowdin charges
- Dipole and quadrupole moments
- Mayer & Wiberg Lowdin bond indices
- Electrostatic potential at nuclei (ESP)

## Methods
Study total energy vs bond distance for H₂:
- Vary methods: scf, b3lyp, mp2, mp4, fci
- Compute energies at 20 initial bond distances
- Interpolate using cubic spline (SciPy) to generate 1000-point smooth curves

## Basis Set Effects
Study total energy vs bond distance for LiH:
- Vary basis sets: cc-pVDZ, cc-pCVDZ, aug-cc-pVDZ, heavy-aug-cc-pVDZ, cc-pVTZ
- Compute energies at 20 bond distances
- Interpolate using cubic spline to 1000 points

## Non-Relativistic vs Relativistic Effects
Compare Hartree-Fock (HF) calculations using:
- ECP (Effective Core Potential)
- DKH (Douglas-Kroll-Hess)
- X2C (Exact Two-Component)

## Geometry Optimization
Example: NH₃
- Optimize geometry with different methods and basis sets
- Store optimized energies and Cartesian coordinates

## Frequency Analysis
Example: CO₂
- Compute vibrational frequencies, IR intensities, and degeneracy
- Identify imaginary frequencies to check for true minima

## Python Libraries & Environment

To run the scripts, install the following:

- [Psi4](https://psicode.org)
- [NumPy](https://numpy.org)
- [SciPy](https://scipy.org)

Recommended using **Conda** environment.