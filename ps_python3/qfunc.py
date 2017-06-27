import numpy
import kernels
from mcfit import SphericalBessel as sph
#mcfit multiplies by sqrt(2/pi)*x**2 to the function. 
#Divide the funciton by this to get the correct form 

from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.misc import derivative
modpath = "/global/u1/c/chmodi/Programs/Py_codes/modules"
import sys
sys.path.append(modpath)

class Qfunc:

    def __init__(self, k, p, Qfile = None, Rfile = None):
        
        self.kp = k
        self.p = p
        self.ipk = interpolate(k, p)
        self.ilpk = self.loginterp(k, p)
        self.renorm = numpy.sqrt(numpy.pi/2.)
        self.tpi2 = 2*numpy.pi**2.
        self.kint = numpy.logspace(-5, 5, 1e4)

        if Qfile is None:
            self.kq, self.Q1, self.Q2, self.Q3, self.Q5, self.Q8, self.Qs2 = self.calc_Q()
        else:
            self.kq, self.Q1, self.Q2, self.Q3, self.Q5, self.Q8, self.Qs2  = numpy.loadtxt(Qfile, unpack=True)
        self.ilQ1 = self.loginterp(self.kq, self.Q1, rp = -5)
        self.ilQ2 = self.loginterp(self.kq, self.Q2, rp = -5)
        self.ilQ3 = self.loginterp(self.kq, self.Q3, rp = -5)
        self.ilQ5 = self.loginterp(self.kq, self.Q5, rp = -5)
        self.ilQ8 = self.loginterp(self.kq, self.Q8, rp = -5)
        self.ilQs2 = self.loginterp(self.kq, self.Qs2, rp = -5)

        if Rfile is None:
            self.kr, self.R1, self.R2 = self.calc_R()
        else:
            self.kr, self.R1, self.R2  = numpy.loadtxt(Rfile, unpack=True)
        self.ilR1 = self.loginterp(self.kr, self.R1)
        self.ilR2 = self.loginterp(self.kr, self.R2)


    def calc_Q(self):
        #k = numpy.logspace(-4, 4, 2e3)
        print('Evaluating Q integrals. Recommend saving them')
        k = numpy.logspace(-4, 4, 2e3)
        p = self.ilpk(k)
        Qk = kernels.Q(k, p, self.ilpk)
        Q1, Q2, Q3, Q5, Q8, Qs2 = numpy.zeros_like(k), numpy.zeros_like(k), numpy.zeros_like(k), numpy.zeros_like(k), numpy.zeros_like(k), numpy.zeros_like(k)
        for foo in range(k.size):
            Q1[foo] = Qk.Q_external(k[foo], 1)
            Q2[foo] = Qk.Q_external(k[foo], 2)
            Q3[foo] = Qk.Q_external(k[foo], 3)
            Q5[foo] = Qk.Q_external(k[foo], 5)
            Q8[foo] = Qk.Q_external(k[foo], 8)
            Qs2[foo] = Qk.Q_external(k[foo], -1)
        return k, Q1, Q2, Q3, Q5, Q8, Qs2

    def calc_R(self):
        print('Evaluating R integrals. Recommend saving them')
        k = numpy.logspace(-4, 4, 2e3)
        p = self.ilpk(k)
        Rk = kernels.R(k, p, self.ilpk)
        R1, R2 = numpy.zeros_like(k), numpy.zeros_like(k)
        for foo in range(k.size):
            R1[foo] = Rk.R_external(k[foo], 1)
            R2[foo] = Rk.R_external(k[foo], 2)
        return k, R1, R2


    def dosph(self, n, x, f, tilt = 1.5, extrap = True):
        f = f*self.renorm
        return sph(x, nu = n, q = tilt)(f, extrap = extrap)
        

    #correlatin function
    def corr(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (1.*self.tpi2)
        return self.dosph(0, kint, integrand, tilt = tilt)

    #Xi integrals from 1506.05264; ZV LEFT paper
    #0 lag
    def xi0lin0(self, kmin = 1e-6, kmax = 1e3):
        val = quad(self.ilpk, kmin, kmax, limit = 200)[0]/self.tpi2
        return val

    def xi0loop0(self, kmin = 1e-6, kmax = 1e3):
        kint = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), 1e4)
        qt = 9./98.*self.ilQ1(kint)
        rt = 10./21.*self.ilR1(kint)
        integrand = qt + rt 
        linterp = self.loginterp(kint, integrand)
        val = quad(linterp, kmin, kmax, limit = 200)[0]/self.tpi2
        return val

    #j0
    def xi0lin(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(0, kint, integrand, tilt = tilt)
    
    def xi0loop(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        qt = 9./98.*self.ilQ1(kint)
        rt = 10./21.*self.ilR1(kint)
        integrand = qt + rt 
        integrand /= (kint**2 *self.tpi2)
        return self.dosph(0, kint, integrand, tilt = tilt)

    def xi0eft(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (self.tpi2)        
        return self.dosph(0, kint, integrand, tilt = tilt)
    
    #j2
    def xi2lin(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(2, kint, integrand, tilt = tilt)

    def xi2loop(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        qt = 9./98.*self.ilQ1(kint)
        rt = 10./21.*self.ilR1(kint)
        integrand = qt + rt 
        integrand /= (kint**2 *self.tpi2)
        return self.dosph(2, kint, integrand, tilt = tilt)

    def xi2eft(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (self.tpi2)        
        return self.dosph(2, kint, integrand, tilt = tilt)
        
    #j1
    def xi1loop(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint        
        qt = self.ilQ1(kint) - 3*self.ilQ2(kint)
        rt = 2.*self.ilR1(kint) - 6*self.ilR2(kint)
        integrand = qt + rt 
        integrand *= (-3./7.)/(kint**3.*self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)

    def xi1eft(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand *= (-3./7.)/(kint *self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)

    #j3
    def xi3loop(self, kint = None, tilt = 1.5):
        #T112
        if kint is None:
            kint = self.kint
        qt = self.ilQ1(kint) + 2*self.ilQ2(kint)
        rt = 2.*self.ilR1(kint) + 4*self.ilR2(kint)
        integrand = qt + rt 
        integrand *= (-3./7.)/(kint**3. *self.tpi2)
        return self.dosph(3, kint, integrand, tilt = tilt)

    def xi3eft(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand *= (-3./7.)/(kint *self.tpi2)
        return self.dosph(3, kint, integrand, tilt = tilt)


    #Integrals from Apprndix B3 of 1209.0780; Carlson CLPT paper
    #V
    def v1(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilR1(kint)
        integrand *= (-3./7.)/(kint**3 *self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)

    def v3(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilQ1(kint)
        integrand *= (-3./7.)/(kint**3 *self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)


    def s1(self, kint = None, tilt = 1.5):
        ### FACTOR OF kq in denominator
        if kint is None:
            kint = self.kint
        integrand = self.ilR1(kint)
        integrand *= (3./7.)/(kint**3 *self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)

    #U
    def u1(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = -1*kint*self.ilpk(kint)
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)
        
    
    def u3(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        #U3 in Martin's code
        integrand = (-5./21)*kint*self.ilR1(kint)
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)

    def u11(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = (-6./7)* kint* (self.ilR1(kint) + self.ilR2(kint))
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)

    def u20(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = (-3./7)* kint* self.ilQ8(kint) 
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)

    #x10&y10 loop
    #These are off by factor of 2 from ZV files
    def x10(self, kint = None):
        lag = self.x100lag()
        j0 =  self.x10j0(kint = kint)
        j2 = self.x10j2(kint = kint)
        return (j0[0], lag + j0[1] + j2[1])

    def x100lag(self, kmin = 1e-6, kmax = 1e3, kint = None):
        if kint is None:
            kint = self.kint
        integrand = 1/7.*(self.ilR1(kint) - self.ilR2(kint))
        linterp = self.loginterp(kint, integrand)
        val = quad(linterp, kmin, kmax, limit = 200)[0]/self.tpi2
        return val

    def x10j0(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = (-1./14)*( 4*self.ilR2(kint) + 2*self.ilQ5(kint))
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(0, kint, integrand, tilt = tilt)

    def x10j2(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = (-1./14)*(3*self.ilR1(kint) + 4*self.ilR2(kint) + 2*self.ilQ5(kint))
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(2, kint, integrand, tilt = tilt)
    
    def y10(self, kint = None, tilt = 1.5):    
        #This is simply -3*x10j2?
        if kint is None:
            kint = self.kint
        integrand = (3./14)*(3*self.ilR1(kint) + 4*self.ilR2(kint) + 2*self.ilQ5(kint))
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(2, kint, integrand, tilt = tilt)
                          
    #T
    def t112(self, kint = None, tilt = 1.5):
        #same as xi3loop
        if kint is None:
            kint = self.kint
        qt = self.ilQ1(kint) + 2*self.ilQ2(kint)
        rt = 2.*self.ilR1(kint) + 4*self.ilR2(kint)
        integrand = qt + rt 
        integrand *= (-3./7.)/(kint**3. *self.tpi2)
        return self.dosph(3, kint, integrand, tilt = tilt)

    #SHEAR integrals
    #from 1609.02908; ZV, EC; CLEFT-GSM, Appendix D
    #Jshear
    def jshear(self, kint = None, tilt = 1.5):    
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (kint**1.*self.tpi2)
        j1 = self.dosph(1, kint, integrand, tilt = tilt)
        j3 = self.dosph(3, kint, integrand, tilt = tilt)
        q = (j1[0] + j3[0])/2.
        js2 = 2*j1[1]/15. - j3[1]/5.
        js3 = -j1[1]/5. - j3[1]/5.
        js4 = j3[1].copy()
        return (q, js2, js3, js4)


    def js2(self, kint = None, tilt = 1.5):    
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (1.*self.tpi2)
        return self.dosph(2, kint, integrand, tilt = tilt)

    def v10(self, kint = None, tilt = 1.5):    
        if kint is None:
            kint = self.kint
        integrand = (-2/7.)*self.ilQs2(kint)
        integrand /= (kint * self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)


    def zeta(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (1.*self.tpi2)
        q,l0 = self.dosph(0, kint, integrand, tilt = tilt) 
        q,l2 = self.dosph(2, kint, integrand, tilt = tilt)
        q,l4 = self.dosph(4, kint, integrand, tilt = tilt)
        toret = 8/35.*l4**2 + 8/63.*l2**2 + 4/45.*l0**2
        return (q, 2*toret)
    


    ### Enterpolate functions in log-sapce beyond the limits
    def loginterp(self, x, y, yint = None, side = "both", lorder = 15, rorder = 15, lp = 1, rp = -1, \
                  ldx = 1e-6, rdx = 1e-6):
        '''Extrapolate function by evaluating a log-index of left & right side
        '''
        if yint is None:
            yint = interpolate(x, y, k = 5)

        if side == "both":
            side = "lr"
            l =lp
            r =rp
        lneff = derivative(yint, x[l], dx = x[l]*ldx, order = lorder)*x[l]/y[l]
        rneff = derivative(yint, x[r], dx = x[r]*rdx, order = rorder)*x[r]/y[r]
        if lneff < 0:
            print('ERROR: Runaway index on left side, bad interpolation. Left index = ', lneff)
        if rneff > 0:
            print('ERROR: Runaway index on right side, bad interpolation. Reft index = ', rneff)
        #print('Log index on left & right edges are = ', lneff, rneff)

        xl = numpy.logspace(-18, numpy.log10(x[l]), 10**6.)
        xr = numpy.logspace(numpy.log10(x[r]), 10., 10**6.)
        yl = y[l]*(xl/x[l])**lneff
        yr = y[r]*(xr/x[r])**rneff

        xint = x[l+1:r].copy()
        yint = y[l+1:r].copy()
        if side.find("l") > -1:
            xint = numpy.concatenate((xl, xint))
            yint = numpy.concatenate((yl, yint))
        if side.find("r") > -1:
            xint = numpy.concatenate((xint, xr))
            yint = numpy.concatenate((yint, yr))
        yint2 = interpolate(xint, yint, k = 5)

        return yint2

