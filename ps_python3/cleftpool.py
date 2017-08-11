import numpy
import numpy as np
from mcfit import SphericalBessel as sph
#mcfit multiplies by sqrt(2/pi)*x**2 to the function. 
#Divide the funciton by this to get the correct form 

from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from qfunc import Qfunc
from itertools import repeat
from functools import partial
import multiprocessing as mp


class CLEFT():
    '''
    Class to evaluate CLEFT Power Spectrum Components given a linear power spectrum k, p.
    Call hzpt.make_table() to create the table of power spectra
    The order is 
    k, ZA, A, W, b1, b1^2, b2, b2^2, b1b2, bs2, b1bs2, b2bs2, bs2^2, bn, b1bn
    #bn and b1bn are not implemented yet
    '''
    
    def __init__(self, k = None, p = None, pfile = None, qfile = None, rfile = None):
        if pfile is None:
            if p is None:
                print("Specify the power sepctrum file or the array")
        else:
            k, p = np.loadtxt(pfile, unpack = True)
        self.kp = k
        self.p = p
        self.qf = Qfunc(k, p, Qfile=qfile, Rfile = rfile)
        print("Q & R kernels created")

        self.renorm = numpy.sqrt(numpy.pi/2.) #mcfit normaliztion
        self.tpi2 = 2*numpy.pi**2.
        self.jn = 10 #number of bessels to sum over

        self.setup_dm()
        print("Matter q-functions created")
        self.setup_blocal()
        print("Bias(local) q-functions created")
        self.setup_bshear()
        print("Shear q-functions created")
    
    def setup_dm(self):
        qf = self.qf

        #Linear
        self.xi00lin = qf.xi0lin0()
        self.qv, xi0lin = qf.xi0lin() #qv determined here
        xi2lin = qf.xi2lin()[1]
        self.Xlin = 2/3.*(self.xi00lin - xi0lin - xi2lin)
        ylinv = 2*xi2lin

        #Since we divide by ylin, check for zeros
        mask = (ylinv == 0)
        ylinv[mask] = interpolate(self.qv[~mask], ylinv[~mask])(self.qv[mask])
        self.Ylin = ylinv

        #Useful variables here
        self.XYlin = (self.Xlin + self.Ylin)
        self.sigma = self.XYlin[-1]
        self.yq = (1*self.Ylin/self.qv)

        #Loop
        self.xi00loop = qf.xi0loop0()
        xi0loop = qf.xi0loop()[1] 
        xi2loop = qf.xi2loop()[1]
        self.Xloop = 2/3.*(self.xi00loop - xi0loop - xi2loop)
        self.Yloop = 2*xi2loop
        self.XYloop = (self.Xloop + self.Yloop)
        self.sigmaloop = self.XYloop[-1]
        
        #Loop2
        self.xi1loop = qf.xi1loop(tilt = 0.5)[1]
        self.xi3loop = qf.xi3loop(tilt = 0.5)[1]

        
    def setup_blocal(self):
        qf = self.qf
        self.x10 = qf.x10()[1]
        self.y10 = qf.y10()[1]
        self.u10 = qf.u1()[1]
        self.u11 = qf.u11()[1]
        self.u20 = qf.u20()[1]
        self.corr = qf.corr()[1]

    def setup_bshear(self):
        qf = self.qf
        js = qf.jshear()
        js2 = qf.js2()

        self.chi = 4/3.*js2[1]**2
        self.v12 = 4*js[1]*js2[1]
        self.x20 = 4*js[2]**2
        self.y20 = 2*(3*js[1]**2 + 4*js[1]*js[2] + 2*js[1]*js[3] + 2*js[2]**2 + 4*js[2]*js[3] + js[3]**2)
        self.v10 = qf.v10()[1]
        self.zeta = qf.zeta()[1]

#####Do Bessel Integrals Here

#Declaring global variables to be able to use pool without
#passing around huge arrays which will have to be pickled
#causing lot of overhead        


def template0(k, func, extrap = False):
    '''Template for j0 integral
    '''
    expon0 = numpy.exp(-0.5*k**2 * (XYlin - sigma))
    suppress = numpy.exp(-0.5*k**2 *sigma)
    Fq = expon0*func(k = k, l= 0)
    if abs(func(k =k, l=0)[-1] ) > 1e-7:
        #print("Subtracting large scale constant = ", func(0)[-1], k)
        sigma2 = func(k = k, l= 0)[-1]
        Fq -= sigma2
    ktemp, ftemp = \
            sph(qv, nu= 0, q=1.5)(Fq*renorm,extrap = extrap)
    ftemp *= suppress
    return 1* numpy.interp(k, ktemp, ftemp)*4*np.pi


def template(k, func, extrap = False):
    '''Template for higher bessel integrals
    '''
    expon = numpy.exp(-0.5*k**2 * (XYlin))
    Fq = 0
    toret = 0
    for l in range(1, jn):
        Fq = expon*func(k = k, l= l)*yq**l
        ktemp, ftemp = \
                sph(qv, nu= l, q=max(0, 1.5 - l))(Fq*renorm, extrap = extrap)
        toret += 1* k**l * numpy.interp(k, ktemp, ftemp)
    return toret*4*np.pi



def integrate(func, pool, taylor = 0):
    '''Do the \mu integrals for all j's by calling the templates
    '''

    p0  = pool.starmap(template0, zip(kv, repeat(func)))
    
    pl  = pool.starmap(template, zip(kv, repeat(func)))
    
    p0, pl = np.array(p0), np.array(pl)

    toret = p0 + pl
    if taylor:
        factorial = np.arange(1, taylor + 1)
        toret *= kv**taylor / factorial.prod()

    return toret


def make_table(cl, kmin = 1e-3, kmax = 3, nk = 100, npool = 2):
    '''Make a table of different terms of P(k) between a given
    'kmin', 'kmax' and for 'nk' equally spaced values in log10 of k
    '''
    header = "k[h/Mpc]   P_Zel   P_A    P_W    P_d    P_dd     P_d^2    P_d^2d^2 \
 P_dd^2    P_s^2    P_ds^2    P_d^2s^2   P_s^2s^2   P_D2d     P_dD2d"

    pktable = numpy.zeros([nk, len(header.split())])

    global kv
    kv = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), nk)

    global qv, XYlin, sigma, yq, renorm, jn
    qv, XYlin, sigma, yq, renorm, jn = cl.qv, cl.XYlin, cl.sigma, cl.yq, cl.renorm, cl.jn

    global Xlin, Ylin, Xloop, Yloop, xi1loop, xi3loop, x10, y10, x20, y20
    Xlin, Ylin, Xloop, Yloop, xi1loop, xi3loop, x10, y10, x20, y20 = \
                cl.Xlin, cl.Ylin, cl.Xloop, cl.Yloop, cl.xi1loop, cl.xi3loop, cl.x10, cl.y10, cl.x20, cl.y20

    global u10, u11, u20, v10, v12, corr, chi, zeta
    u10, u11, u20, v10, v12, corr, chi, zeta = cl.u10, cl.u11, cl.u20, cl.v10, cl.v12, cl.corr, cl.chi, cl.zeta

    
    pool = mp.Pool(npool)

    pktable[:, 0] = kv[:]        
    pktable[:, 1] = integrate(func = za, pool = pool)
    pktable[:, 2] = integrate(aloop, taylor = 2, pool = pool)
    pktable[:, 3] = integrate(wloop, taylor = 3, pool = pool)
    pktable[:, 4] = integrate(b1, pool = pool)
    pktable[:, 5] = integrate(b1sq, pool = pool)
    pktable[:, 6] = integrate(b2, pool = pool)
    pktable[:, 7] = integrate(b2sq, pool = pool)
    pktable[:, 8] = integrate(b1b2, pool = pool)
    pktable[:, 9] = integrate(bs2, pool = pool)
    pktable[:,10] = integrate(b1bs2, pool = pool)
    pktable[:,11] = integrate(b2bs2, pool = pool)
    pktable[:,12] = integrate(bs2sq, pool = pool)
    pktable[:,13] = np.zeros_like(kv)
    pktable[:,14] = np.zeros_like(kv)

    pool.close()


    del kv

    del qv, XYlin, sigma, yq, renorm, jn

    del Xlin, Ylin, Xloop, Yloop, xi1loop, xi3loop, x10, y10, x20, y20

    del u10, u11, u20, v10, v12, corr, chi, zeta

    return pktable

    #pool.join()


###################################################################################
#Functions from table in Appendix B
def za( k, l):
    return np.ones_like(qv)

def aloop( k, l):
    return  -(Xloop + Yloop - 2*l*Yloop/Ylin/k**2.)  

def wloop( k, l):
    return  bool(l)*(6/5.*xi1loop - 6/5.*xi3loop  + 
                      6*(l-1)*xi3loop/(k**2 *Ylin))/yq /k

def b1( k, l):
    #return lambda l: -k**2 *(x10 + y10) + 2*l*y10/Ylin -2*qv* bool(l)*u10/Ylin
    #Ad-hoc factor of 2 to match ZV files
    return  -k**2 *(x10 *2. + y10 *2) + 2*l*y10*2/Ylin -2*qv* bool(l)*u10/Ylin

def b1sq(  k, l):
    return  corr - k**2 *u10**2 + 2*l*u10**2/Ylin -qv* bool(l)*u11/Ylin

def b2(  k, l):
    return   -k**2 *u10**2 + 2*l*u10**2/Ylin -qv* bool(l)*u20/Ylin

def b2sq(  k, l):
    return   corr**2 /2.

def b1b2(  k, l):
    return   -2 *bool(l)*qv*u10*corr/Ylin

def bs2(  k, l):
    return  - k**2 *(x20 + y20) + 2*l*y20/Ylin -2*qv* bool(l)*v10/Ylin

def b1bs2(  k, l):
    return   -1*qv* bool(l)*v12/Ylin

def b2bs2( k, l):
    return   chi

def bs2sq(  k, l):
    return  zeta

