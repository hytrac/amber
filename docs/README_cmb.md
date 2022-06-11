# About

The CMB module computes tau(z), P_ee(k), P_qq(k), C_tau(l), and C_ksz(l). There are also subroutines to make Healpix maps of patchy tau and KSZ. Note electron ionization fraction x_e = n_e/n_e,tot.


# Source

cmb.f90, cmbreion.f90


# Input

CMB parameters are set in input/input.txt.

1) make = make, write, or null
2) mapmake = make, write, or null
3) zmin = minimum redshift
4) zmax = maximum redshift
5) zdel = redshift interval
6) zspacing = linear or logarithm redshift spacing
7) lmin = minimum multipole
8) lmax = maximum multipole
9) nside = Healpix nside
10) dir = output directory


# Output

output/cmb/tau.txt is the Thomson optical depth tau(z).

First line is a header. The three columns are:
1) z = redshift
2) xe = electron ionization fraction x_e = n_e/n_e,tot
3) tau = integrated Thomson optical depth

<br>

output/cmb/power_z=xx.xx.txt are power spectra P(k). 

Note electron ionization fraction x_e = n_e/n_e,tot.

First line is a header. The various columns are:
1) k = comoving wavenumber [h/Mpc]
2) P_mm = matter power [(Mpc/h)^3]
3) P_ee = electron power for tau [(Mpc/h)^3]
4) P_qq = electron projected momentum power for KSZ [(Mpc/h)^3(km/s)^2]
5) b_em = electron-matter bias
6) r_em = electron-matter cross correlation

<br>

output/cmb/cl_tau_z=xx.xx.txt is the patchy tau angular power spectrum C_l.

First line is a header. The two columns are:
1) l = multipole
2) C = angular power

<br>

output/cmb/cl_ksz_z=xx.xx.txt is the patchy KSZ angular power spectrum C_l.

First line is a header. The two columns are:
1) l = multipole
2) C = angular power

<br>

output/cmb/map_tau_nside=xxxx.fits is the patchy tau Healpix map.

output/cmb/map_ksz_nside=xxxx.fits is the patchy KSZ Healpix map.


# References

Chen, Trac, Mukherjee, Cen, Patchy Kinetic Sunyaev-Zel'dovich Effect with Controlled Reionization History and Morphology, 2022, submitted to ApJ, [arXiv:2203.04337](https://arxiv.org/abs/2203.04337)

Dvorkin & Smith, Reconstructing Patchy Reionization from the Cosmic Microwave Background, 2009, [arXiv:0812.1566](https://arxiv.org/abs/0812.1566)

Park et al, The Kinetic Sunyaev-Zel'dovich Effect as a Probe of the Physics of Cosmic Reionization: The Effect of Self-regulated Reionization, 2013, [arXiv:1301.3607](https://arxiv.org/abs/1301.3607)
