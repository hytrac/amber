# AMBER

## Source

amber.f90 is the main source file

input.txt is the main parameter file

<br> 

## Compiling, executing

Makefile is for compiling the code

To compile from scratch:  make clean; make

To execute and save to log file:  ./amber.x < input/input.txt > output/log

<br>

## Cosmology

Cosmological parameters are set in input.txt

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

## Lagrangian Perturbation Theory (LPT)

LPT is used to evolve particles in space and time. 2LPT is recommended.

<br>

output/lpt/power_linear.txt is the output linear power spectrum

First line is a header. The three columns are:
1) k = comoving wavenumber [h/Mpc]
2) Delta^2_lin = dimensionless matter power (linear theory)
3) Delta^2_lpt = dimensionless matter power (lpt realization)

<br>

## Excursion Set Formalism (ESF)

ESF is used to model the halo mass density field. ESF-L with sharp k-space filter is recommended.

<br>

## Reionization


