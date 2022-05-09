# CMB

The CMB module computes tau(z), P_ee(k), P_qq(k), C_tau(l), and C_ksz(l). Note electron ionization fraction x_e = n_e/n_e,tot.

We will add subroutines for making Healpix maps soon.

<br>

## Source

cmb.f90, cmbreion.f90

<br>

## Input

CMB parameters are set in input/input.txt.

1) make = make, write, or null
2) zmin = minimum redshift
3) zmax = maximum redshift
4) zdel = redshift interval
5) zspacing = linear or logarithm redshift spacing
6) lmin = minimum multipole
7) lmax = maximum multipole
8) dir = output directory

<br>

## Output

output/cmb/tau.txt is the Thomson optical depth tau(z).

First line is a header. The three columns are:
1) z = redshift
2) xe = electron ionization fraction x_e = n_e/n_e,tot
3) tau = integrated Thomson optical depth

<br>

output/cmb/power_z=xx.xx.txt are power spectra P(k). 

Note electron ionization fraction x_e = n_e/n_e,tot.

First line is a header. The three columns are:
1) k = comoving wavenumber [h/Mpc]
2) Pee = electron power for tau [(Mpc/h)^3]
3) Pqq = electron projected momentum power for KSZ [(Mpc/h)^3(km/s)^2]

<br>

output/cmb/Cl_tau_z=xx.xx.txt is the patchy tau angular power spectrum C_l.

First line is a header. The two columns are:
1) l = multipole
2) C = angular power

<br>

output/cmb/Cl_ksz_z=xx.xx.txt is the patchy KSZ angular power spectrum C_l.

First line is a header. The two columns are:
1) l = multipole
2) C = angular power
