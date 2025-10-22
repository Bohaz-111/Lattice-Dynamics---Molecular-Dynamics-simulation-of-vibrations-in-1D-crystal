# Simulation of vibrations in a 1D crystal

Classical and semiclassical MD code for a 1D optical-phonon chain that computes the vibrational density of states (VDOS) from velocity autocorrelation functions (VACFs). Dynamics are run in the NVT ensemble using the Bussi–Donadio–Parrinello (CSVR) thermostat. The effective potential for centroid molecular dynamics (CMD) is calculated using the method in Feynman and Kleinert's paper on effective classical partition function.

## Features
- 1D chain with periodic boundary conditions; velocity-Verlet integrator.
- CSVR thermostat (canonical sampling through velocity rescaling) with user-set time constant $$\tau$$.
- CMD evolves the atoms on Feynman-Kleinert effective potential, the method is based on path integral formulation of quantum mechanics.
- CSPD samples each normal mode at an effective temperature, the method is based on the diagonal coherent state representation of Kubo-transformed momentum operator. (see attached CSPD document for detail)
- Real-space <-> normal-mode transforms; VACF $$C_k(t)$$; Fourier transform to spectra $$S_k(\omega)$$ and total VDOS.
- Exports frequency–VDOS data for plotting.

## Model (optical phonon with on-site anharmonicity)
Default units set mass and lattice spacing to 1. The Hamiltonian (Cyclic boundary condition: $$q_{N+1} = q_1$$) is:

$$H = \sum_{j=1}^N \Bigl[ \frac{p_j^2}{2} + h q_j^2 + g q_j^4 \Bigr] + \frac{1}{2} \sum_{j=1}^N (q_j - q_{j+1})^2$$

## Quick start
1. Build with a modern Fortran compiler (e.g. `gfortran`).
2. Adjust `N`, `dt`, `steps`, `tau`, `h`, `g` in the source or input file.
3. Run the executable; plot `VDOS.dat` (frequency vs VDOS) with your preferred tool.

## Outputs
- `VCAF.dat` — two columns: time (t), VCAF.
- `VDOS.dat` — two columns: frequency (`omega`), VDOS.

## Reference
- G. Bussi, D. Donadio, M. Parrinello, **Canonical sampling through velocity rescaling**, *J. Chem. Phys.* **126**, 014101 (2007). DOI: https://doi.org/10.1063/1.2408420
- R. P. Feynman, H. Kleinert, **Effective classical partition functions**, *Phys. Rev. A* **34**, 5080–5084 (1986). DOI: https://doi.org/10.1103/PhysRevA.34.5080
