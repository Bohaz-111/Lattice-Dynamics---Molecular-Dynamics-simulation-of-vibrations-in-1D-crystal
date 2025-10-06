# Lattice-Dynamics — Molecular-Dynamics simulation of vibrations in a 1D crystal

Classical MD code for a 1D optical-phonon chain that computes the vibrational density of states (VDOS) from velocity autocorrelation functions (VACFs). Dynamics are run in the NVT ensemble using the Bussi–Donadio–Parrinello (CSVR) thermostat.

## Features
- 1D chain with periodic boundary conditions; velocity-Verlet integrator.
- CSVR thermostat (canonical sampling through velocity rescaling) with user-set time constant $$\tau$$.
- Real-space <-> normal-mode transforms; VACF $$_k(t)$$; Fourier transform to spectra $$S_k(\omega)$$ and total VDOS.
- Exports frequency–VDOS data for plotting.

## Model (optical phonon with on-site anharmonicity)
Default units set mass and lattice spacing to 1. The Hamiltonian (Cyclic boundary condition: $$q_{N+1} = q_1$$) is:

$$H = \sum_{j=1}^N \Bigl[ \frac{p_j^2}{2} + h q_j^2 + g q_j^4 \Bigr] + \frac{1}{2} * \sum_{j=1}^N (q_j - q_{j+1})^2$$

## Quick start
1. Build with a modern Fortran compiler (e.g. `gfortran`). If enabled, link FFTW3 for fast transforms.
2. Adjust `N`, `dt`, `steps`, `tau`, `h`, `g` in the source or input file.
3. Run the executable; plot `VDOS.dat` (frequency vs VDOS) with your preferred tool.

## Outputs
- `VCAF.dat` — two columns: time (t), VCAF.
- `VDOS.dat` — two columns: frequency (`omega`), VDOS.

## Reference (thermostat)
- G. Bussi, D. Donadio, M. Parrinello, **Canonical sampling through velocity rescaling**, *J. Chem. Phys.* **126**, 014101 (2007). DOI: https://doi.org/10.1063/1.2408420

> TL;DR: NVT MD of an anharmonic optical chain (quartic on-site + harmonic coupling) using CSVR; VDOS obtained from normal-mode VACFs.
