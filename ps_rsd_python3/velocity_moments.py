#
from __future__ import print_function,division

import numpy as np
#import kernels
from mcfit import SphericalBessel as sph
#mcfit multiplies by sqrt(2/pi)*x**2 to the function. 
#Divide the funciton by this to get the correct form 

from matplotlib import pyplot as plt

from scipy.integrate import quad
from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.misc import derivative
import sys

class VelocityMoments:
    '''
    Class to evaluate the Zeldovich power spectrum, given a linear power
    spectrum k, p [in compatible units, e.g. h/Mpc and (Mpc/h)^3].
    This can be used as the basis for both the HaloZeldovich and ZEFT
    approximations, and for auto and cross-spectra.
    See:
    Modi, White & Vlah, arXiv:1706.03173
    https://arxiv.org/abs/1706.03173
    for more information.
    To main method is make_table(), which creates the table of power spectra
    components.  The order is k, ZA, b1, b2, b1sq, b2sq, b1b1
    Convenience functions are provided for common calls.
    '''
    def __init__(self, k, p, f=1):
        '''k,p are the linear theory power spectra in compatible units,
        e.g. h/Mpc and (Mpc/h)^3.
            f is the growth-factor derivative'''
        self.kp    = k
        self.p     = p
        self.f     = f
        self.ilpk  = self.loginterp(k, p)
        self.renorm=np.sqrt(np.pi/2.) #mcfit normaliztion
        self.tpi2  = 2*np.pi**2.
        self.qt    = np.logspace(-5, 5, 1e4)
        self.kint  = np.logspace(-5, 5, 1e4)
        self.jn    = 5 #number of bessels to sum over
        
        # set up velocity table
        self.pktable = None
        self.num_power_components = 3
        
        self.vktable = None
        self.num_velocity_components = 4
        
        self.setup()
        #
    def setup(self):
        '''
        Create X_L, Y_L, xi_L, U1_L \& 0lag sigma.
        '''
        self.xi0lag = self.xi0lin0() 
        self.qv, xi0v = self.xi0lin()
        xi2v = self.xi2lin()[1]
        self.corlin = self.corr()[1]
        self.Ulin = self.u10lin()[1]
        #
        self.Xlin = 2/3.*(self.xi0lag - xi0v - xi2v)
        ylinv = 2*xi2v
        #Since we divide by ylin, check for zeros
        mask = (ylinv == 0)
        ylinv[mask] = interpolate(self.qv[~mask], ylinv[~mask])(self.qv[mask])
        self.Ylin = ylinv
        self.XYlin = (self.Xlin + self.Ylin)
        self.sigma = self.XYlin[-1]
        self.yq = (1*self.Ylin/self.qv)
    
        # time derivatives (currently only at lowest order)
        self.Udot = self.f * self.Ulin
        self.Xdot = self.f * self.Xlin
        self.Ydot = self.f * self.Ylin
    
    ### Interpolate functions in log-sapce beyond the limits
    def loginterp(self, x, y, yint = None, side = "both",\
                  lorder = 15, rorder = 15, lp = 1, rp = -1, \
                  ldx = 1e-6, rdx = 1e-6):
        '''
        Extrapolate function by evaluating a log-index of left & right side
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
        #
        xl = np.logspace(-18, np.log10(x[l]), 10**6.)
        xr = np.logspace(np.log10(x[r]), 10., 10**6.)
        yl = y[l]*(xl/x[l])**lneff
        yr = y[r]*(xr/x[r])**rneff
        #
        xint = x[l+1:r].copy()
        yint = y[l+1:r].copy()
        if side.find("l") > -1:
            xint = np.concatenate((xl, xint))
            yint = np.concatenate((yl, yint))
        if side.find("r") > -1:
            xint = np.concatenate((xint, xr))
            yint = np.concatenate((yint, yr))
        yint2 = interpolate(xint, yint, k = 5)
        #
        return yint2

    def dosph(self, n, x, f, tilt = 1.5):
        #Function to do bessel integral using FFTLog for kernels
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
    def template(self, k, l, func, expon, suppress, power=1, za = False, expon_za = 1.):
        ''' Simplified vs the original code. Beta.
            Generic template that is followed by mu integrals
        j0 is different since its exponent has sigma subtracted that is
        later used to suppress integral
        '''
        
        Fq = np.zeros_like(self.qv)
        
        if za == True and l == 0:
            Fq = expon_za * func * self.yq**l
        else:
            Fq = expon * func * self.yq**l
        
        # note that the angular integral for even powers of mu gives J_(l+1)
        ktemp, ftemp = sph(self.qv, nu= l+(power%2), q=max(0,1.5-l))(Fq*self.renorm,extrap = False)
        
        #        if l or ( power%2 == 1 ): # note that the angular integral for even powers of mu gives J_(l+1)
        #    Fq = expon*func*self.yq**l
        #elif j0:
        #    Fq = expon0*func
        #
        #if power == 1:
        #    ktemp, ftemp = \
        #        sph(self.qv, nu= l+1, q=max(0, 1.5 - l))(Fq*self.renorm,\
        #        extrap = False)
        #            #else:
        #            #ktemp, ftemp = \
        #        sph(self.qv, nu= l, q=max(0, 1.5 - l))(Fq*self.renorm,\
        #                                                 extrap = False)

        ftemp *= suppress


        return 1* k**l * np.interp(k, ktemp, ftemp)

    
    def p_integrals(self, k):
        '''Since the templates function is modified this currently contains only a subset of terms.
            
            '''
        ksq = k**2
        expon = np.exp(-0.5*ksq * (self.XYlin - self.sigma))
        exponm1 = np.expm1(-0.5*ksq * (self.XYlin - self.sigma))
        suppress = np.exp(-0.5*ksq *self.sigma)

        #
        za, b1, b2 = 0, 0, 0
        
        #l indep functions
        fza = 1.
        fb1 = -2 * k * self.Ulin
        #fb2sq = (self.corlin**2/2.)
        #fb1b2 = (-2*self.qv*self.Ulin*self.corlin/self.Ylin)
        #
        
        for l in range(self.jn):
            #l-dep functions
            fb2 = (2*l/self.Ylin/ksq - 1) * self.Ulin**2
            #fb1sq = (self.corlin + (2*l/self.Ylin - k**2)*self.Ulin**2)
            #fb2 = ((2*l/self.Ylin - k**2)*self.Ulin**2)
            
            #do integrals
            za += self.template(k,l,fza,expon,suppress,power=0,za=True,expon_za=exponm1)
            b1 += self.template(k,l,fb1,expon,suppress,power=1)
            b2 += self.template(k,l,fb2,expon,suppress,power=0)
            #Xdot += self.template(k,l,fXdot,expon,expon0,suppress,j0=True,power=0)
            #Ydot += self.template(k,l,fYdot,expon,expon0,suppress,j0=True,power=1)
        
        
        return 4*np.pi*np.array([za,b1,b2])

    
    
    def v_integrals(self, k):
        '''Do the \mu integrals for various parameters for give 'k'
            (Currently for the pairwise velocity in k space.)
            (Commented out sections come from equivalent expressions for density power spectrum.)
        '''
        ksq = k**2
        expon = np.exp(-0.5*ksq * (self.XYlin - self.sigma))
        exponm1 = np.expm1(-0.5*ksq * (self.XYlin - self.sigma))
        suppress = np.exp(-0.5*ksq *self.sigma)

        # Note: collect all constant offset terms into 'offset'
        A_const, A_ldep, b1, offset = 0, 0, 0, 0
        
        #l indep functions
        fb1 = self.Ulin
        fA_const = k * (self.Xdot + self.Ydot - self.sigma)
        foffset = k * self.sigma

        for l in range(self.jn):
            #l-dep functions
            fA_ldep = k * ( - 2*l/ksq/self.Ylin) * self.Ydot
            
            #do integrals
            b1 += self.template(k,l,fb1,expon,suppress,power=1)
            A_const += self.template(k,l,fA_const,expon,suppress,power=0)
            A_ldep += self.template(k,l,fA_ldep,expon,suppress,power=0)
            offset += self.template(k,l,foffset,expon,suppress,power=0,za=True,expon_za=exponm1)

        
        return 4*np.pi*np.array([A_const, A_ldep, b1, offset])

    def make_ptable(self, kmin = 1e-3, kmax = 3, nk = 100):
        '''Make a table of different terms of P(k) between a given
        'kmin', 'kmax' and for 'nk' equally spaced values in log10 of k
        This is the most time consuming part of the code.
        '''
        self.pktable = np.zeros([nk, self.num_power_components+1]) # one column for ks
        kv = np.logspace(np.log10(kmin), np.log10(kmax), nk)
        self.pktable[:, 0] = kv[:]
        for foo in range(nk):
            self.pktable[foo, 1:] = self.p_integrals(kv[foo])



    def make_vtable(self, kmin = 1e-3, kmax = 3, nk = 100):
        '''Make a table of different terms of P(k) between a given
        'kmin', 'kmax' and for 'nk' equally spaced values in log10 of k
        This is the most time consuming part of the code.
        '''
        self.vktable = np.zeros([nk, self.num_velocity_components+1]) # one column for ks
        kv = np.logspace(np.log10(kmin), np.log10(kmax), nk)
        self.vktable[:, 0] = kv[:]
        for foo in range(nk):
            self.vktable[foo, 1:] = self.v_integrals(kv[foo])

    

    def write_table(self,fn):
        '''Writes the table to an ascii text file, fn.'''
        if self.vktable is None:
            print("Zeldovich table not created.")
            return
        # The order is k, b1, Xdot, Ydot, ....
        hdr= ["k","b1","Xdot","Ydot"]
        ff = open(fn,"w")
        ff.write("# Components of the pairwise velocity\n")
        str = "# %14s"%hdr[0]
        for hh in hdr[1:]:
            str += " %15s"%hh
        ff.write(str+"\n")
        for i in range(self.vktable.shape[0]):
            str = ""
            for t in self.vktable[i,:]:
                str += " %15.5e"%t
            ff.write(str+"\n")
        ff.close()


def loglog_abs(x,y,color='k'):
    '''Routine to make loglog plots where negative values are indicated by dashed lines.
        Mostly for my own convenience.'''
    
    plt.loglog(x,y,color)
    plt.loglog(x,-y,color+'--')
    return None


if __name__ == '__main__':
    pks = np.loadtxt('pklin.test')
    velda = VelocityMoments(pks[:,0],pks[:,1])
    velda.make_ptable()
    velda.make_vtable()


