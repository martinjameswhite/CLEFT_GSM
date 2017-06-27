import numpy
import numpy as np
from mcfit import SphericalBessel as sph
#mcfit multiplies by sqrt(2/pi)*x**2 to the function. 
#Divide the funciton by this to get the correct form 

from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from qfunc import Qfunc


class CLEFT():
    '''
    Class to evaluate HZ Power Spectrum given a linear power spectrum k, p.
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

        
    #################
    #Bessel Integrals for \mu

    def template0(self, k, func, extrap = False):
        '''Template for j0 integral
        '''
        expon0 = numpy.exp(-0.5*k**2 * (self.XYlin - self.sigma))
        suppress = numpy.exp(-0.5*k**2 *self.sigma)
        Fq = expon0*func(0) 
        if abs(func(0)[-1] ) > 1e-7:
            #print("Subtracting large scale constant = ", func(0)[-1], k)
            sigma2 = func(0)[-1]
            Fq -= sigma2
        ktemp, ftemp = \
            sph(self.qv, nu= 0, q=1.5)(Fq*self.renorm,extrap = extrap)
        ftemp *= suppress
        return 1* numpy.interp(k, ktemp, ftemp)*4*np.pi

    def template(self, k, func, extrap = False):
        '''Template for higher bessel integrals
        '''
        expon = numpy.exp(-0.5*k**2 * (self.XYlin))
        Fq = 0
        toret = 0
        for l in range(1, self.jn):
            Fq = expon*func(l)*self.yq**l
            ktemp, ftemp = \
                sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm, extrap = extrap)
            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)
        return toret*4*np.pi


    def integrate(self, func, taylor = 0):
        '''Do the \mu integrals for all j's by calling the templates
        '''
        p0  = np.array([self.template0(k, func(k)) for k in self.kv])
        pl  = np.array([self.template(k, func(k)) for k in self.kv])
        toret = p0 + pl
        if taylor:
            factorial = np.arange(1, taylor + 1)
            toret *= self.kv**taylor / factorial.prod()
        return toret


    def make_table(self, kmin = 1e-3, kmax = 3, nk = 100):
        '''Make a table of different terms of P(k) between a given
        'kmin', 'kmax' and for 'nk' equally spaced values in log10 of k
        '''
        header = "k[h/Mpc]   P_Zel   P_A    P_W    P_d    P_dd     P_d^2    P_d^2d^2 \
     P_dd^2    P_s^2    P_ds^2    P_d^2s^2   P_s^2s^2   P_D2d     P_dD2d"
        
        self.pktable = numpy.zeros([nk, len(header.split())])
        self.kv = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), nk)

        self.pktable[:, 0] = self.kv[:]        
        self.pktable[:, 1] = self.integrate(self.za)
        self.pktable[:, 2] = self.integrate(self.aloop, 2)
        self.pktable[:, 3] = self.integrate(self.wloop, 3)
        self.pktable[:, 4] = self.integrate(self.b1)
        self.pktable[:, 5] = self.integrate(self.b1sq)
        self.pktable[:, 6] = self.integrate(self.b2)
        self.pktable[:, 7] = self.integrate(self.b2sq)
        self.pktable[:, 8] = self.integrate(self.b1b2)
        self.pktable[:, 9] = self.integrate(self.bs2)
        self.pktable[:,10] = self.integrate(self.b1bs2)
        self.pktable[:,11] = self.integrate(self.b2bs2)
        self.pktable[:,12] = self.integrate(self.bs2sq)
        self.pktable[:,13] = np.zeros_like(self.kv)
        self.pktable[:,14] = np.zeros_like(self.kv)
        


    ###################################################################################
    #Functions from table in Appendix B
    def za(self, k):
        return lambda l: np.ones_like(self.qv)

    def aloop(self, k):
        return lambda l: -(self.Xloop + self.Yloop - 2*l*self.Yloop/self.Ylin/k**2.)  

    def wloop(self, k):
        return lambda l: bool(l)*(6/5.*self.xi1loop - 6/5.*self.xi3loop  + 
                          6*(l-1)*self.xi3loop/(k**2 *self.Ylin))/self.yq /k

    def b1(self, k):
        return lambda l: -k**2 *(self.x10 + self.y10) + 2*l*self.y10/self.Ylin -2*self.qv* bool(l)*self.u10/self.Ylin

    def b1sq(self, k):
        return lambda l: self.corr - k**2 *self.u10**2 + 2*l*self.u10**2/self.Ylin -self.qv* bool(l)*self.u11/self.Ylin

    def b2(self, k):
        return lambda l:  -k**2 *self.u10**2 + 2*l*self.u10**2/self.Ylin -self.qv* bool(l)*self.u20/self.Ylin

    def b2sq(self, k):
        return lambda l:  self.corr**2 /2.
    
    def b1b2(self, k):
        return lambda l:  -2 *bool(l)*self.qv*self.u10*self.corr/self.Ylin

    def bs2(self, k):
        return lambda l: - k**2 *(self.x20 + self.y20) + 2*l*self.y20/self.Ylin -2*self.qv* bool(l)*self.v10/self.Ylin
 
    def b1bs2(self, k):
        return lambda l:  -1*self.qv* bool(l)*self.v12/self.Ylin
        
    def b2bs2(self, k):
        return lambda l:  self.chi
    
    def bs2sq(self, k):
        return lambda l:  self.zeta
    
