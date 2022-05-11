# About

The Abundance Matching Box for the Epoch of Reionization (AMBER) is a semi-numerical code for modeling the cosmic dawn and reionization. The new algorithm is not based on the excursion set formalism for reionization, but takes the novel approach of calculating the reionization-redshift field assuming that hydrogen gas encountering higher radiation intensity are photoionized earlier. Redshift values are assigned while matching the abundance of ionized mass according to a given mass-weighted ionization fraction. The code has the unique advantage of allowing users to directly specify the reionization history through the redshift midpoint, duration, and asymmetry input parameters. The reionization process is further controlled through the minimum halo mass for galaxy formation and the radiation mean free path for radiative transfer. There are improved methods for constructing density, velocity, halo, and radiation fields, which are essential components for modeling reionization observables.


# Install

AMBER is written in modern Fortran and parallelized using OpenMP to run on a multicore, shared-memory node. In addition to modules for simulating reionization, there are preliminary modules to compute CMB and 21cm observables.

Contents:
- Source files in [src](https://github.com/hytrac/amber/tree/main/src)
- Documentation READMEs in [docs](https://github.com/hytrac/amber/tree/main/docs)
- Example inputs in [examples](https://github.com/hytrac/amber/tree/main/examples)

Requirements:
- [Intel ifort compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.zjewa9)
- [Intel Math Kernel Library (MKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) for FFT, random numbers, spline interpolation

# Contact

Please submit comments and questions to [issues](https://github.com/hytrac/amber/issues).

The current lead developers are Hy Trac (hytrac@andrew.cmu.edu) and Nianyi Chen (nianyic@andrew.cmu.edu).


# References

H. Trac, N. Chen, I. Holst, M.A. Alvarez, R. Cen, AMBER: A Semi-numerical Abundance Matching Box for the Epoch of Reionization, 2022, [ApJ, 927, 186](https://iopscience.iop.org/article/10.3847/1538-4357/ac5116), [arXiv:2109.10375](https://arxiv.org/abs/2109.10375)

N. Chen, H. Trac, S. Mukherjee, R. Cen, Patchy Kinetic Sunyaev-Zel'dovich Effect with Controlled Reionization History and Morphology, 2022, submitted to ApJ, [arXiv:2203.04337](https://arxiv.org/abs/2203.04337)
