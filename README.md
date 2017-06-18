# CLEFT_GSM

This code implements the Gaussian Streaming Model using components from
Convolution Lagrangian Effective Field Theory as described in:

Z.Vlah, E.Castorina, M.White

The Gaussian streaming model and Convolution Lagrangian effective field theory

JCAP 12(2016)007, [https://arxiv.org/abs/1609.02908]

The code is written (mostly) in C++.  It can be run from the command line, or
called from Python (wrappers provided).

The C++ version currently only implements the configuration-space statistics
(i.e. the correlation function).  The Fortran routines in "ps_fortran" provide
an implementation of the power spectrum routines.

We additionally provide fast and simple Python routines for computing the ZEFT
and Halo-Zeldovich models' real-space auto- and cross-correlations of biased tracers
as described in

C.Modi, M.White, Z.Vlah

Modeling CMB Lensing Cross Correlations with CLEFT

JCAP, submitted. [https://arxiv.org/abs/1706.03173]
