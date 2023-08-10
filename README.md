## AtomicBethe: a program for calculation of Bethe logarithm for many-electron atoms in the mean-field approximation

### Citations and references

This code is a companion to the paper:
- (will be filled later)

Some additional key literature references:
- C. Schwartz, Phys. Rev. 123, 1700 (1961);
- K. Pachucki and J. Komasa, Phys. Rev. Lett. 92, 213001 (2004);

The program uses the open-shell Hartree-Fock theory:
- C. C. J. Roothaan, Rev. Mod. Phys. 32, 179 (1960);

### Installation

To install the program, you require:
- a C++ compiler compliant with the C++17 standard (or newer);
- Boost library, https://www.boost.org/;
- Eigen library, https://eigen.tuxfamily.org/index.php?title=Main_Page;

The latter two libraries are available in repositories of most Linux-based systems. Moreover, they are header-only libraries, so alternatively you can just download them and include their directory in the Makefile (via `-I/path_to_library/ option`). Besides, the program uses the mINI library (https://github.com/pulzed/mINI/tree/master) which is included along the source code. The program is compiled using the `make` command; you need to create `object` directory manually before compiling.

### Usage

Examples of the input files are available in the `sample-inputs` directory. The inputs have the standard structure of an `INI` file:

```
[section1]
key1 = value1
key2 = value2

[section2]
key3 = value3
```

Below, each section and the corresponding keywords relevant for the users are listed with a brief explanation.

Section `[job]`:
- `mode`: 0 means basis set optimization, 1 means automatic calculation of the Bethe logarithm, 2 means calculation of the response function for a given value of `k`;
Section `[orbital_basis]`:
- `orbital_momentum`: number of functions of `s`, `p`, `d`, `...` symmetries (separated by spaces);
- `basis_params_s`, `basis_params_p`, `basis_params_d`, and so forth: basis set tempering parameters (geometric progression);
Section `[system]`:
- `nuclear_charge`: charge of the atomic nucleus;
- `n_electrons`: number of electrons in the atom;
- `n_closed`: number of closed (doubly-occupied) orbitals;
- `n_open`: number of open orbitals;
Section `[thresh]`:
- `linear_dependent`: cutoff for dropping linearly dependent functions from the basis;
- `eri_threshold`: cutoff for dropping negligible integrals;
Section `[scf]`:
- `maxit`: maximum number of SCF iteractions;
- `n_diis`: size of the DIIS subspace;
- `n_diis_turn_on`: turn on DIIS after this number of iterations;
- `rohf_coupling_f`, `rohf_coupling_a`, `rohf_coupling_b` - coupling parameters defining the atomic state;
