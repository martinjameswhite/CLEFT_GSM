import numpy
#import kernels
from mcfit import SphericalBessel as sph
#mcfit multiplies by sqrt(2/pi)*x**2 to the function. 
#Divide the funciton by this to get the correct form 

from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.misc import derivative
modpath = "/global/u1/c/chmodi/Programs/Py_codes/modules"
import sys
sys.path.append(modpath)

class Hzpt:
    '''
    Class to evaluate HZ Power Spectrum given a linear power spectrum k, p.
    Call hzpt.make_table() to create the table of power spectra
    The order is k, ZA, b1, b2, b1sq, b2sq, b1b1
    '''
    def __init__(self, k, p, Qfile = None, Rfile = None):
        
        self.kp = k
        self.p = p
        self.ilpk = self.loginterp(k, p)
        self.renorm = numpy.sqrt(numpy.pi/2.) #mcfit normaliztion
        self.tpi2 = 2*numpy.pi**2.
        self.qt = numpy.logspace(-5, 5, 1e4)
        self.kint = numpy.logspace(-5, 5, 1e4)
        self.jn = 10 #number of bessels to sum over
        self.setup()

    def setup(self):
        '''
        Create X_L, Y_L, xi_L, U1_L \& 0lag sigma
        '''
        self.xi0lag = self.xi0lin0() 
        self.qv, xi0v = self.xi0lin()
        xi2v = self.xi2lin()[1]
        self.corlin = self.corr()[1]
        self.Ulin = self.u10lin()[1]
        
        self.Xlin = 2/3.*(self.xi0lag - xi0v - xi2v)
        ylinv = 2*xi2v
        #Since we divide by ylin, check for zeros
        mask = (ylinv == 0)
        ylinv[mask] = interpolate(self.qv[~mask], ylinv[~mask])(self.qv[mask])
        self.Ylin = ylinv
        self.XYlin = (self.Xlin + self.Ylin)
        self.sigma = self.XYlin[-1]
        self.yq = (1*self.Ylin/self.qv)
        


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
        print('Log index on left & right edges are = ', lneff, rneff)

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

        
    #################
    #Function to do bessel integral using FFTLog for kernels
    def dosph(self, n, x, f, tilt = 1.5):
        f = f*self.renorm
        return sph(x, nu = n, q = tilt)(f, extrap = True)

    #PT kernels below
    #0 lag
    def xi0lin0(self, kmin = 1e-6, kmax = 1e3):
        val = quad(self.ilpk, kmin, kmax, limit = 200)[0]/self.tpi2
        return val


    #j0
    def xi0lin(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(0, kint, integrand, tilt = tilt)
    
    
    #j2
    def xi2lin(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(2, kint, integrand, tilt = tilt)

    #u1
    def u10lin(self, kint = None,  tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = -1*kint*self.ilpk(kint)
        integrand /= (kint**2.*self.tpi2)
        return self.dosph(1, kint, integrand, tilt = tilt)

    #correlatin function
    def corr(self, kint = None, tilt = 1.5):
        if kint is None:
            kint = self.kint
        integrand = self.ilpk(kint)
        integrand /= (1.*self.tpi2)
        return self.dosph(0, kint, integrand, tilt = tilt)
        

    #################
    #Bessel Integrals for \mu

    def template(self, k, l, func, expon, expon0, suppress, j0 = True):
        '''Generic template that is followed by mu integrals
        j0 is different since its exponent has sigma subtracted that is
        later used to suppress integral
        '''
        Fq = 0
        if l:
            Fq = expon*func*self.yq**l
        elif j0:
            Fq = expon0*func
            
        ktemp, ftemp = \
            sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
        if l==0 and j0:
            ftemp *= suppress
        return 1* k**l * numpy.interp(k, ktemp, ftemp)

    def integrals(self, k):
        '''Do the \mu integrals for various parameters for give 'k'
        '''
        expon = numpy.exp(-0.5*k**2 * self.XYlin)
        expon0 = numpy.exp(-0.5*k**2 * (self.XYlin - self.sigma))
        expon0m1 = numpy.expm1(-0.5*k**2 * (self.XYlin - self.sigma))
        suppress = numpy.exp(-0.5*k**2 *self.sigma)
        
        za, b1, b2, b1sq, b2sq, b1b2 = 0, 0, 0, 0, 0, 0
        #l indep functions
        fza = 1.
        fb2sq = (self.corlin**2/2.)
        fb1 = (-2*self.qv*self.Ulin/self.Ylin)
        fb1b2 = (-2*self.qv*self.Ulin*self.corlin/self.Ylin)
        
        for l in range(self.jn):
            #l-dep functions
            fb1sq = (self.corlin + (2*l/self.Ylin - k**2)*self.Ulin**2)
            fb2 = ((2*l/self.Ylin - k**2)*self.Ulin**2)
            #do integrals
            za += self.template(k, l, fza, expon, expon0m1, suppress, j0 = True)
            b1 += self.template(k, l, fb1, expon, expon0, suppress, j0 = False)
            b2 += self.template(k, l, fb2, expon, expon0, suppress, j0 = True)
            b1sq += self.template(k, l, fb1sq, expon, expon0, suppress, j0 = True)
            b2sq += self.template(k, l, fb2sq, expon, expon0, suppress, j0 = True)
            b1b2 += self.template(k, l, fb1b2, expon, expon0, suppress, j0 = False)
            
        return 4*numpy.pi*numpy.array([za, b1, b2, b1sq, b2sq, b1b2])


    def make_table(self, kmin = 1e-3, kmax = 3, nk = 100):
        '''Make a table of different terms of P(k) between a given
        'kmin', 'kmax' and for 'nk' equally spaced values in log10 of k
        '''
        self.pktable = numpy.zeros([nk, 7])
        kv = numpy.logspace(numpy.log10(kmin), numpy.log10(kmax), nk)
        self.pktable[:, 0] = kv[:]

        for foo in range(nk):
            self.pktable[foo, 1:] = self.integrals(kv[foo])



    
    #################
    #Bessel Integrals for \mu
    #Explicit function definitions below
    ##### CAN BE DELETED WHEN CONFIDENT ########

    def za_integral(self, k):

        Fq = 0
        toret = 0    
        expon = numpy.exp(-0.5*k**2 * self.XYlin)
        expon0 = numpy.expm1(-0.5*k**2 * (self.XYlin - self.sigma))
        suppress = numpy.exp(-0.5*k**2 *self.sigma)
        for l in range(self.jn):
            func = 1.
            if l:
                Fq = expon*func*self.yq**l
            else:
                Fq = expon0*func
            
            ktemp, ftemp = \
                sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
            if l == 0:
                ftemp *= suppress
            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)

##        Fq = 0
##        toret = 0    
##        for l in range(self.jn):
##            if l:
##                Fq = numpy.exp(-0.5*k**2 * self.XYlin)*(1*Ylin/qv)**l 
##            else :
##                Fq = numpy.expm1(-0.5*k**2 * (XYlin - sigma))
##            
##            ktemp, ftemp = \
##                sph(qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
##            if l == 0:
##                ftemp *= numpy.exp(-0.5*k**2*sigma)
##            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)
##    
        return toret


    def b1_integral(self, k):

        Fq = 0
        toret = 0    
        expon = numpy.exp(-0.5*k**2 * self.XYlin)
        expon0 = numpy.exp(-0.5*k**2 * (self.XYlin - self.sigma))
        suppress = numpy.exp(-0.5*k**2 *self.sigma)
        for l in range(self.jn):
            func = (-2*self.qv*self.Ulin/self.Ylin)
            if l:
                Fq = expon*func*self.yq**l
            #else:
            #    Fq = expon0*func
            
            ktemp, ftemp = \
                sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
            if l == 0:
                ftemp *= suppress
            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)

##        Fq = 0
##        toret = 0    
##        for l in range(self.jn):
##            if l:
##                Fq = numpy.exp(-0.5*k**2 * XYlin)*(-2*qv*ulin/Ylin)*(1*Ylin/qv)**l 
##            
##            ktemp, ftemp = \
##                sph(qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
##            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)
##    
        return toret


    def b1sq_integral(self, k):

        Fq = 0
        toret = 0    
        expon = numpy.exp(-0.5*k**2 * self.XYlin)
        expon0 = numpy.exp(-0.5*k**2 * (self.XYlin - self.sigma))
        suppress = numpy.exp(-0.5*k**2 *self.sigma)
        for l in range(self.jn):
            func = (self.corlin + (2*l/self.Ylin - k**2)*self.Ulin**2)
            if l:
                Fq = expon*func*self.yq**l
            else:
                Fq = expon0*func
            
            ktemp, ftemp = \
                sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
            if l == 0:
                ftemp *= suppress
            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)

##        Fq = 0
##        toret = 0    
##        for l in range(self.jn):
##            if l:
##                Fq = numpy.exp(-0.5*k**2 * XYlin)*(cor + (2*l/Ylin - k**2)*ulin**2)*(1*Ylin/qv)**l 
##            else:
##                Fq = numpy.exp(-0.5*k**2 *( XYlin -sigma))*(cor + (2*l/Ylin - k**2)*ulin**2)
##            ktemp, ftemp = \
##                sph(qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
##            if l == 0:
##                ftemp *= numpy.exp(-0.5*k**2*sigma)
##            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)
##    
        return toret



    def b2_integral(self, k):

        Fq = 0
        toret = 0    
        expon = numpy.exp(-0.5*k**2 * self.XYlin)
        expon0 = numpy.exp(-0.5*k**2 * (self.XYlin - self.sigma))
        suppress = numpy.exp(-0.5*k**2 *self.sigma)
        for l in range(self.jn):
            func = ((2*l/self.Ylin - k**2)*self.Ulin**2)
            if l:
                Fq = expon*func*self.yq**l
            else:
                Fq = expon0*func
            
            ktemp, ftemp = \
                sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
            if l == 0:
                ftemp *= suppress
            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)

##        for l in range(self.jn):
##            if l:
##                Fq = numpy.exp(-0.5*k**2 * XYlin)*((2*l/Ylin - k**2)*ulin**2)*(1*Ylin/qv)**l 
##            else:
##                Fq = numpy.exp(-0.5*k**2 *( XYlin -sigma))*(- ulin**2)
##            
##            ktemp, ftemp = \
##                sph(qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
##            if l == 0:
##                ftemp *= k**2 * numpy.exp(-0.5*k**2*sigma)
    
        return toret


    def b2sq_integral(self, k):

        Fq = 0
        toret = 0    
        expon = numpy.exp(-0.5*k**2 * self.XYlin)
        expon0 = numpy.exp(-0.5*k**2 * (self.XYlin - self.sigma))
        suppress = numpy.exp(-0.5*k**2 *self.sigma)
        for l in range(self.jn):
            func = (self.corlin**2/2.)
            if l:
                Fq = expon*func*self.yq**l
            else:
                Fq = expon0*func
            
            ktemp, ftemp = \
                sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
            if l == 0:
                ftemp *= suppress
            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)
            
        return toret



    def b1b2_integral(self, k):

        Fq = 0
        toret = 0    
        expon = numpy.exp(-0.5*k**2 * self.XYlin)
        for l in range(self.jn):
            func = (-2*self.qv*self.Ulin*self.corlin/self.Ylin)
            if l:
                Fq = expon*func*self.yq**l 
            
            ktemp, ftemp = \
                sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,extrap = False)
            toret += 1* k**l * numpy.interp(k, ktemp, ftemp)
        return toret
        


##
##    def setupinterpolatef(self):        
##        qt = self.qt
##        xi0lin0 = self.xi0lin0()
##        ilxi0lin = self.loginterp(*self.xi0lin())
##        ilxi2lin = self.loginterp(*self.xi2lin(tilt = 0))
##
##        Xlin = 2/3.*( ilxi0lin(qt) + ilxi2lin(qt))
##        iltXlin = self.loginterp(qt, Xlin)
##
##        ilXlin = lambda x: 2/3.*xi0lin0 - iltXlin(x)
##
##        #xlin goes negative. Check that and change extrapolate
##        #check if this is actually needed
##        qcheck = numpy.logspace(-5, -2, 1e4)
##        qstop = qcheck[::-1][numpy.where(ilXlin(qcheck[::-1]) < 0)[0][0]]
##        qforX = numpy.logspace(numpy.log10(qstop), 8, 1e5)
##        itXlin = lambda x: numpy.interp(x, qforX, ilXlin(qforX))
##
##        #ylin
##        Ylin = 2*ilxi2lin(qt)
##        ilYlin = self.loginterp(qt, Ylin)
##        #itYlin = lambda x: numpy.interp(x, qforX, ilYlin(qforX))
##
##        #ulin
##        ilulin = self.loginterp(*self.u10lin(), lp = 200)
##        #return  itXlin, ilYlin, ilulin
##        return  ilXlin, ilYlin, ilulin
##
##
##
