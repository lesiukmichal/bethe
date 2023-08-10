## AtomicBethe: a program for calculation of Bethe logarithm for many-electron atoms in the mean-field approximation

### Citations and references

This code is a companion to the paper:
- (will be filled later)

Some additional key literature references:
- C. Schwartz, Phys. Rev. 123, 1700 (1961);
- K. Pachucki and J. Komasa, Phys. Rev. Lett. 92, 213001 (2004);

### Installation

To install the program, you require:
- a C++ compiler compliant with the C++17 standard;
- Boost library, https://www.boost.org/;
- Eigen library, https://eigen.tuxfamily.org/index.php?title=Main_Page;

The latter two libraries are available in repositories of most Linux-based systems. Moreover, they are header-only libraries, so alternatively you can just download them and include their directory in the Makefile (via -I/path_to_library/ option). Besides,
the program uses the mINI library (https://github.com/pulzed/mINI/tree/master) which is included along the source code.
The program is compiled using the "make" command; you need to create "object" directory manually before compiling.
