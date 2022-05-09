# AMBER

## Source

amber.f90 is the main source file.

input.txt is the main parameter file.

<br> 

## Compiling

Makefile is for compiling the code.

Requires Intel ifort compiler and Math Kernel Library (MKL).

<br>

To compile from scratch:  make clean; make

To execute and save to log file:  ./amber.x < input/input.txt > output/log

<br>

## Simulation

Sim input parameters:
1) ncore = number of cores (OpenMP threads)
2) input = input directory
3) output = output directory

<br>

## Cosmology

Cosmo input parameters:
1) Lbox = comoving sim box length [Mpc/h]
2) om = matter density omega_m
3) ol = lambda density omega_l
4) ob = baryon density omega_b
5) or = radiation density omega_r
6) h = hubble constant h_0
7) s8 = sigma_8
8) ns = spectral index n_s
9) w = dark energy equation of state
20) Tcmb = CMB temperature [uK]
21) XH = hydrogen mass fraction
22) YHe = helium mass frction
23) linpowspec filename = linear power spectrum file
24) dir = output directory


<br>

input/linpowspec.txt is the input linear power spectrum

First line is a header. The two columns are:
1) k = comoving wavenumber [h/Mpc]
2) Delta^2 = dimensionless matter power

<br>

output/cosmo/power_z=xx.xx.txt are output power spectra

First line is a header. The various columns are:
1) k = comoving wavenumber [h/Mpc]
2) D^2_lin = theory matter power linearly evolved to redshift z
3) D^2_lpt = initial matter power linearly evolved to redshift z
4) D^2_mm = matter power
5) D^2_hh = halo power
6) D^2_zz = reionization-redshift power
7) b_hm = halo-matter bias
8) r_hm = halo-matter cross correlation
9) b_zm = redshift-matter bias
10) r_zm = redshift-matter cross correlation
11) b_zh = redshift-halo bias
12) r_zh = redshift-halo cross correlation

<br>

## Reionization

Reion input parameters:
1) make = make, read, or write
2) z_mid = redshift midpoint
3) z_del = redshift duration
4) z_asy = redshift asymmetry
5) x_early = ionization fraction early
6) x_mid = ionization fraction mid
7) x_late = ionization fraction late
8) M_min = minimum halo mass [Msolar/h]
9) l_mfp = RT comoving mean free path [Mpc/h]
10) dir = output directory

<br>

output/reion/zre.dat is a single-precision binary file containing the zre field.

<br>

## Gaussian Random Field (GRF)

GRF is used to generate initial conditions.

GRF input parameters:
1) make = make, read, or write
2) seed = integer seed
3) dir = output directory

<br>

output/grf/grf.dat is a double-precision file containing the Fourier transform of the GRF.

<br>

## Lagrangian Perturbation Theory (LPT)

LPT is used to generate particles.

LPT input parameters:
1) make = make, read, or write
2) order = 1 or 2 (2lpt is recommended)
3) assign = tsc, cic, or ngp (tsc is recommended)
4) dir = output directory

<br>

output/lpt/power_linear.txt is the output linear power spectrum

First line is a header. The three columns are:
1) k = comoving wavenumber [h/Mpc]
2) Delta^2_lin = dimensionless matter power (linear theory)
3) Delta^2_lpt = dimensionless matter power (lpt realization)

<br>

## Excursion Set Formalism (ESF)

ESF is used to model the halo mass density field.

ESF input parameters:
1) make = make, read, or write
2) filter = sharpk or tophat (sharpk is recommended)
3) assign = tsc, cic, or ngp (tsc is recommended)
4) dir = output directory

<br>

output/esf/rhoh1.dat is a single-precision binary file containing the halo mass density field.

output/esf/rhoh2.dat is a single-precision binary file containing the interaced halo mass density field.

<br>

## Mesh/Part

Mesh/part input parameters:
1) make = make, read, or write
2) nm1d = number of cells/particles per side length
3) dir = output directory

<br>

output/mesh/rho1.dat is a single-precision binary file containing the matter density field.

output/mesh/rho2.dat is a single-precision binary file containing the interaced matter density field.
