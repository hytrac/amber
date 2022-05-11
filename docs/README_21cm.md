# About
The 21cm module contains example subroutines to approximately compute the global brightness temperature T_b(z) and power spectrum P_TT(k), assuming T_s >> T_cmb and ignoring peculiar vel, lightcone effects.


# Source

h21cm.f90, h21cmreion.f90


# Input

21cm parameters are set in input/input.txt.

1) make = make, write, or null
2) zmin = minimum redshift
3) zmax = maximum redshift
4) zdel = redshift interval
5) zspacing = linear or logarithmic redshift spacing
6) dir = output directory


# Output

output/21cm/21cm_global.txt is the global brightness temperature T_b(z).

First line is a header. The three columns are:
1) z = redshift
2) x_H = neutral hydrogen fraction (mass-weighted)
3) T_b = temperature brightness [mK]

<br>

output/21cm/power_z=xx.xx.txt is the power spectrum P(k).

First line is a header. The two columns are:
1) k = comoving wavenumber [h/Mpc]
2) P = temperature power [(mK)^2(Mpc/h)^3]


# References
