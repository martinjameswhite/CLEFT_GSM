import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.integrate import trapz, simps
from itertools import product
import multiprocessing as mp
## Define functions

#constants

#Define Q:
def Q1approx( r, x):
    return 0.25*(1 + x)**2

def Q1( r, x):
    y =  1 + r**2 - 2*r*x
    return (r**2 *(1 - x**2)**2)/y**2

def Q2approx( r, x):
    return 0.25*x*(1 + x)

def Q2( r, x):
    y =  1 + r**2 - 2*r*x
    return ((1- x**2)*r*x*(1 - r*x))/y**2

def Q3approx( r, x):
    return 0.25*x**2

def Q3( r, x):
    y =  1 + r**2 - 2*r*x
    return (x**2 *(1 - r*x)**2)/y**2

def Q5approx( r, x):
    y =  1 + r**2 - 2*r*x
    return x*(1 + x)/2.

def Q5( r, x):
    y =  1 + r**2 - 2*r*x
    return r*x*(1 - x**2.)/y

def Q8approx( r, x):
    y =  1 + r**2 - 2*r*x
    return (1 + x)/2.

def Q8( r, x):
    y =  1 + r**2 - 2*r*x
    return r**2.*(1 - x**2.)/y

def Qs2approx( r, x):
    y =  1 + r**2 - 2*r*x
    return (1 + x)*(1- 3*x)/4.

def Qs2( r, x):
    y =  1 + r**2 - 2*r*x
    return r**2.*(x**2 - 1.)*(1 - 2*r**2 + 4*r*x - 3*x**2)/y**2

dictQapprox = {1:Q1approx, 2:Q2approx, 3:Q3approx, 5:Q5approx, 8:Q8approx, -1:Qs2approx}
dictQ = {1:Q1, 2:Q2, 3:Q3, 5:Q5, 8: Q8, -1:Qs2}


#####Define double integral here
def Q_internal(k, n, r, approx = False):
        
    if approx:
        func = lambda r, x: ilpk(k*numpy.sqrt(1 + r**2 - 2*r*x))*dictQapprox[n](r, x)
    else:
        func = lambda r, x: ilpk(k*numpy.sqrt(1 + r**2 - 2*r*x))*dictQ[n](r, x)
    return (glwval*func(r, glxval)).sum(axis = -1)

def Q_external(n, k):
    fac = k**3 /(2*numpy.pi)**2
    r = (kint/k)
    absr1 = abs(r-1)
    mask = absr1 < tol
    y = numpy.zeros_like(r)
    y[~mask]  = Q_internal(k, n, r[~mask].reshape(-1, 1))
    if mask.sum():
        y[mask] = Q_internal(k, n,  r[mask].reshape(-1, 1), approx = True)
    y *= pkint
    return fac*trapz(y, r)


def Q(kv, pk, npool = 4, ns = None, kintv = None):

    global ilpk, tol, kint, pkint, glxval, glwval

    ilpk = pk
    if kintv is None:
        kint = numpy.logspace(-6, 3, 1e3)
    else:
        kint = kintv
    pkint = ilpk(kint)

    tol = 10**-5. 
    #glxval, glwval = numpy.loadtxt("gl_128.txt", unpack = True)
    glxval, glwval = numpy.array(gl_128_X), numpy.array(gl_128_W)
    glxval = glxval.reshape(1, -1)
    glwval = glwval.reshape(1, -1)

    if ns is None:
        ns = [1, 2, 3, 5, 8, -1]

    pool = mp.Pool(npool)
    prod = product(ns, kv)
    Qv = pool.starmap(Q_external, list(prod))
    pool.close()
    pool.join()
    Qv = numpy.array(Qv).reshape(len(ns), -1)
    toret = numpy.zeros([Qv.shape[0] + 1, Qv.shape[1]])
    toret[0] = kv
    toret[1:, :] = Qv

    del ilpk, tol, kint, pkint, glxval, glwval

    return toret

###############################################################


def R1approx(r, x):
    return 0.5*(1- x)*(1 + x)**2

def R1(r, x):
    y =  1 + r**2 - 2*r*x
    return (r**2 *(1 - x**2)**2)/y

def R2approx(r, x):
    return 0.5*x*(1 + x)*(1-x)

def R2(r, x):
    y =  1 + r**2 - 2*r*x
    return ((1- x**2)*r*x*(1 - r*x))/y

listRapprox = [R1approx, R2approx]
listR = [R1, R2]


#Double integral here

def R_internal(n, r, approx = False):

    if approx:
        func = lambda r, x: listRapprox[n-1](r, x)
    else:
        func = lambda r, x: listR[n-1](r, x)
    return (glwval*func(r, glxval)).sum(axis = -1)

def R_external(n, k):

    fac = k**3 /(2*numpy.pi)**2 *ilpk(k)
    r = (kint/k)
    absr1 = abs(r-1)
    mask = absr1 < tol
    y = numpy.zeros_like(r)
    y[~mask]  = R_internal(n, r[~mask].reshape(-1, 1))
    if mask.sum():
        y[mask] = R_internal(n,  r[mask].reshape(-1, 1), approx = True)
    y *= pkint
    return fac*trapz(y, r)



def R(kv, pk, npool = 4, ns = None, kintv = None):

    global ilpk, tol, kint, pkint, glxval, glwval

    ilpk = pk
    if kintv is None:
        kint = numpy.logspace(-6, 3, 1e3)
    else:
        kint = kintv
    pkint = ilpk(kint)

    tol = 10**-5. 
    #glxval, glwval = numpy.loadtxt("gl_128.txt", unpack = True)
    glxval, glwval = numpy.array(gl_128_X), numpy.array(gl_128_W)
    glxval = glxval.reshape(1, -1)
    glwval = glwval.reshape(1, -1)

    if ns is None:
        ns = [1, 2]

    pool = mp.Pool(npool)
    prod = product(ns, kv)
    Rv = pool.starmap(R_external, list(prod))
    pool.close()
    pool.join()
    Rv = numpy.array(Rv).reshape(len(ns), -1)
    toret = numpy.zeros([Rv.shape[0] + 1, Rv.shape[1]])
    toret[0] = kv
    toret[1:, :] = Rv

    del ilpk, tol, kint, pkint, glxval, glwval

    return toret


#####################################################################################


gl_128_X = [9.99825e-01, 9.99077e-01, 9.97733e-01, 9.95793e-01, 9.93257e-01, 9.90128e-01, 9.86407e-01, 9.82096e-01, 9.77198e-01, 9.71717e-01, 9.65654e-01, 9.59015e-01, 9.51802e-01, 9.44020e-01, 9.35674e-01, 9.26769e-01, 9.17310e-01, 9.07303e-01, 8.96753e-01, 8.85668e-01, 8.74053e-01, 8.61915e-01, 8.49263e-01, 8.36103e-01, 8.22443e-01, 8.08292e-01, 7.93657e-01, 7.78548e-01, 7.62974e-01, 7.46944e-01, 7.30468e-01, 7.13554e-01, 6.96215e-01, 6.78459e-01, 6.60298e-01, 6.41742e-01, 6.22802e-01, 6.03490e-01, 5.83818e-01, 5.63797e-01, 5.43438e-01, 5.22755e-01, 5.01760e-01, 4.80464e-01, 4.58881e-01, 4.37025e-01, 4.14906e-01, 3.92540e-01, 3.69940e-01, 3.47118e-01, 3.24088e-01, 3.00865e-01, 2.77463e-01, 2.53894e-01, 2.30174e-01, 2.06316e-01, 1.82334e-01, 1.58244e-01, 1.34059e-01, 1.09794e-01, 8.54636e-02, 6.10820e-02, 3.66638e-02, 1.22237e-02, -9.99825e-01, -9.99077e-01, -9.97733e-01, -9.95793e-01, -9.93257e-01, -9.90128e-01, -9.86407e-01, -9.82096e-01, -9.77198e-01, -9.71717e-01, -9.65654e-01, -9.59015e-01, -9.51802e-01, -9.44020e-01, -9.35674e-01, -9.26769e-01, -9.17310e-01, -9.07303e-01, -8.96753e-01, -8.85668e-01, -8.74053e-01, -8.61915e-01, -8.49263e-01, -8.36103e-01, -8.22443e-01, -8.08292e-01, -7.93657e-01, -7.78548e-01, -7.62974e-01, -7.46944e-01, -7.30468e-01, -7.13554e-01, -6.96215e-01, -6.78459e-01, -6.60298e-01, -6.41742e-01, -6.22802e-01, -6.03490e-01, -5.83818e-01, -5.63797e-01, -5.43438e-01, -5.22755e-01, -5.01760e-01, -4.80464e-01, -4.58881e-01, -4.37025e-01, -4.14906e-01, -3.92540e-01, -3.69940e-01, -3.47118e-01, -3.24088e-01, -3.00865e-01, -2.77463e-01, -2.53894e-01, -2.30174e-01, -2.06316e-01, -1.82334e-01, -1.58244e-01, -1.34059e-01, -1.09794e-01, -8.54636e-02, -6.10820e-02, -3.66638e-02, -1.22237e-02]


gl_128_W = [4.49381e-04, 1.04581e-03, 1.64250e-03, 2.23829e-03, 2.83275e-03, 3.42553e-03, 4.01625e-03, 4.60458e-03, 5.19016e-03, 5.77264e-03, 6.35166e-03, 6.92689e-03, 7.49798e-03, 8.06459e-03, 8.62638e-03, 9.18301e-03, 9.73415e-03, 1.02795e-02, 1.08187e-02, 1.13514e-02, 1.18773e-02, 1.23961e-02, 1.29076e-02, 1.34113e-02, 1.39070e-02, 1.43943e-02, 1.48731e-02, 1.53430e-02, 1.58037e-02, 1.62550e-02, 1.66966e-02, 1.71281e-02, 1.75495e-02, 1.79603e-02, 1.83604e-02, 1.87496e-02, 1.91275e-02, 1.94940e-02, 1.98489e-02, 2.01919e-02, 2.05228e-02, 2.08414e-02, 2.11476e-02, 2.14412e-02, 2.17219e-02, 2.19897e-02, 2.22443e-02, 2.24857e-02, 2.27135e-02, 2.29278e-02, 2.31284e-02, 2.33152e-02, 2.34881e-02, 2.36469e-02, 2.37916e-02, 2.39220e-02, 2.40382e-02, 2.41400e-02, 2.42273e-02, 2.43002e-02, 2.43586e-02, 2.44024e-02, 2.44316e-02, 2.44462e-02, 4.49381e-04, 1.04581e-03, 1.64250e-03, 2.23829e-03, 2.83275e-03, 3.42553e-03, 4.01625e-03, 4.60458e-03, 5.19016e-03, 5.77264e-03, 6.35166e-03, 6.92689e-03, 7.49798e-03, 8.06459e-03, 8.62638e-03, 9.18301e-03, 9.73415e-03, 1.02795e-02, 1.08187e-02, 1.13514e-02, 1.18773e-02, 1.23961e-02, 1.29076e-02, 1.34113e-02, 1.39070e-02, 1.43943e-02, 1.48731e-02, 1.53430e-02, 1.58037e-02, 1.62550e-02, 1.66966e-02, 1.71281e-02, 1.75495e-02, 1.79603e-02, 1.83604e-02, 1.87496e-02, 1.91275e-02, 1.94940e-02, 1.98489e-02, 2.01919e-02, 2.05228e-02, 2.08414e-02, 2.11476e-02, 2.14412e-02, 2.17219e-02, 2.19897e-02, 2.22443e-02, 2.24857e-02, 2.27135e-02, 2.29278e-02, 2.31284e-02, 2.33152e-02, 2.34881e-02, 2.36469e-02, 2.37916e-02, 2.39220e-02, 2.40382e-02, 2.41400e-02, 2.42273e-02, 2.43002e-02, 2.43586e-02, 2.44024e-02, 2.44316e-02, 2.44462e-02]
