# Power spectrum code

Python code to compute the real-space auto- and cross-correlations of biased tracers
for the ZEFT, Halo-Zeldovich and GSM models as described in

C.Modi, M.White, Z.Vlah

Modeling CMB Lensing Cross Correlations with CLEFT

JCAP, 08(2017)009, [https://arxiv.org/abs/1706.03173]

The underlying theory is based upon Convolution Lagrangian Effective Field
Theory, as described in:

Z.Vlah, E.Castorina, M.White

The Gaussian streaming model and Convolution Lagrangian effective field theory

JCAP 12(2016)007, [https://arxiv.org/abs/1609.02908]



The code is parallelized using 'pool' object in python. 

The main code is in cleftpool.py, which has CLEFT class to create PT kernels. make_table function
in the same file uses these kernels to create a table of P(k) where different columns correspond
to contribution of different bias parameters. The CLEFT class takes in an 'order' argument which can
be 1 or 2 (default) to do ZA or 1-loop perturbation theory.


- We provide a script - main.py - which can be directly run as follows: <br>
python main.py --pfile path_to_linear_ps_file

- Other arguments that can be provided are can be seen by calling: <br>
python main.py --help

- Some of the important arguments are - <br>
'pfile' is the Linear Power Spectrum file, 'npool' is the number of cores, 'z' is the redshift, 
'M' is Omega-Matter, 'nk' is the number of log-spaced k-values between 'k_min' and 'k_max'
at which power spectrum is evaluated.

- For jupyter notebooks, the most basic call for the code is: <br>
import cleftpool <br>
cl = cleftpool.CLEFT(pfile = pfile,  npool=32, order=2) <br>
pk = cleftpool.make_table(cl, kmin = 0.002, kmax = 1, nk = 200, npool=32, z = 1, M = 0.3)

- Timing: With 32 cores on a single node of Cori-jupyter hub, it takes ~35 seconds to
compute power spectra at 200 k-values.


We have also provided a C++ version of the tree-level code (zeldovich.cpp), though this is far less thoroughly debugged and less robust in its current implementation (see test_cpp.sh for an example) largely because it does not careful treat extrapolation of "short" P(k) input files.  


A fast Python package (VelociLPTors) to compute real- and redshift-space
power spectra and correlation functions using LPT is also available at

https://github.com/sfschen/velocilptors
