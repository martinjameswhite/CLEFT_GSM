# Halo-Zeldovich power spectrum code

Python(3) code to compute the Halo-Zeldovich power spectrum of biased
tracers, as described in Appendix B of Modi, White & Vlah (2017):

Modi, White & Vlah, JCAP, submitted.
https://arxiv.org/abs/1706.03173

This code can be used to model the Fourier, real-space, auto- and
cross-correlations of biased tracers at large scales and high redshift,
for example in modeling cross-correlations of galaxies or quasars
and CMB lensing or quasar-galaxy cross-correlations.

This code requires NumPy and SciPy and makes use of the
"multiplicative convolutional fast integral transforms"
library:

https://github.com/eelregit/mcfit

which you will need to install.
