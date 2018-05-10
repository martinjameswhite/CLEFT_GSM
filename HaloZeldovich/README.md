# Halo-Zeldovich/ZEFT power spectrum code

Python(3) code to compute the Halo-Zeldovich and ZEFT models' power spectra
of biased tracers, as described in Appendix B of Modi, White & Vlah (2017):

Modi, White & Vlah,
JCAP, 08(2017)009, [https://arxiv.org/abs/1706.03173]

This code can be used to model the Fourier, real-space, auto- and
cross-correlations of biased tracers at large scales and high redshift,
for example in modeling cross-correlations of galaxies or quasars
and CMB lensing or quasar-galaxy cross-correlations.  It implements
the core calculations required for the Zeldovich model, the ZEFT model
and the Halo-Zeldovich model (as described in the paper).

This code requires NumPy and SciPy and makes use of the
"multiplicative convolutional fast integral transforms"
library:

https://github.com/eelregit/mcfit

which you will need to install.
