# Exact Diagonalization of a Spin Hamiltonian
Diagonalization of spin-1/2 J-Q Hamiltonian on a 2D periodic rectangular lattice, implemented in Julia.

Uses Mz conservation, translational symmetry, and inversion symmetry to reduce the Hilbert space dimension, and either eigen() package or Lanczos Diagonalization based on input parameter "Sparse" to find eigenvalues and eigenstates.

## Requirements
Julia 1.10.5 or advanced

## Some Points
The "main" code uses path based storing for efficient operation. Please change the paths accordingly on your system.

The kx, ky in input fields are actually mx, my (kx = 2 pi mx/Nx)

## Author
* **Himadri Halder**