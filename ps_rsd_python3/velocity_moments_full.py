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

from cleftpool import CLEFT

class VelocityMoments(CLEFT):
    '''
    Class to evaluate the cleft velocity moments.
    
    '''
    def __init__(self, f=1, **kw):
        '''k,p are the linear theory power spectra in compatible units,
        e.g. h/Mpc and (Mpc/h)^3.
            f is the growth-factor derivative'''
        
        # do the X, Y, etc. integrals
        CLEFT.__init__(self,**kw)
        self.f = f 
        
        # set up velocity table
        self.pktable = None
        self.num_power_components = 5
        
        self.vktable = None
        self.num_velocity_components = 12
        
        self.sparktable = None
        self.num_spar_components = 9
        
        self.stracektable = None
        self.num_strace_components = 9
        
        self.setup_dm()
        self.setup_time_derivatives()
    
        #
    def setup_time_derivatives(self):
        '''
        Create time derivatives. The factors of n in time derivatives is annoying...
        '''
        
        qf = self.qf
        
        # time derivatives (currently only at lowest order)
        self.Udot = self.f * self.Ulin
        
        self.Xdot = self.f * self.Xlin
        self.Xdotdot = self.f**2 * self.Xlin
        
        self.Ydot = self.f * self.Ylin
        self.Ydotdot = self.f**2 * self.Ylin
    
        # add one loop
        print("Computing loops for time derivatives.")
        xi0loop013 = qf.xi0loop013()
        xi0loop13 = qf.xi0loop13()[1]
        xi2loop13 = qf.xi2loop13()[1]
        
        xi0loop022 = qf.xi0loop022()
        xi0loop22 = qf.xi0loop22()[1]
        xi2loop22 = qf.xi2loop22()[1]
        
        # important note: the X13 terms coded in cleftpool correspond to X13 + X31 = 2*X13; this is different from the Zvonomir's GSM paper, explaining the factors of two...
        
        self.Xdot1loop = self.f * 2./3 * (2 * (xi0loop013 - xi0loop13 - xi2loop13) + 2*(xi0loop022 - xi0loop22 - xi2loop22 ) )
        self.Ydot1loop = self.f * 2 * ( 2*xi2loop13 + 2*xi2loop22 )

        self.X10dot = self.f * 3./2 * self.x10
        self.Y10dot = self.f * 3./2 * self.y10
        
        self.X10ddot = self.f**2 * 2 * self.x10
        self.Y10ddot = self.f**2 * 2 * self.y10
        
        self.Udot1loop = self.f * 3 * self.u30
        
        self.U11dot = self.f * 2 * self.u11
        self.U20dot = self.f * 2 * self.u20

        self.Xddot1loop = self.f**2 * 2./3 * (3 * (xi0loop013 - xi0loop13 - xi2loop13) + 4*(xi0loop022 - xi0loop22 - xi2loop22 ) )
        self.Yddot1loop = self.f**2 * 2 * ( 3*xi2loop13 + 4*xi2loop22   )
        
        #self.v1 = qf.v1()[1] + qf.S()[1]/self.qv
        #self.v3 = qf.v3()[1] + qf.S()[1]/self.qv
        # 'alt' version has better convergence due to noncancellation of terms
        self.v1 = qf.v1alt()[1] + qf.Salt()[1]
        self.v3 = qf.v3alt()[1] + qf.Salt()[1]
        self.v  = 2 * self.v1 + self.v3
        self.T112  = qf.t112(tilt=0.5)[1]
        
        self.xi0loop013 = xi0loop013
        self.xi0loop13 = xi0loop13
        self.xi2loop13 = xi2loop13

        self.xi0loop022 = xi0loop022
        self.xi0loop22 = xi0loop22
        self.xi2loop22 = xi2loop22

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
    def template(self, k, l, func, expon, suppress, power=1, za = False, expon_za = 1.,tilt=None):
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
        
        if tilt is not None:
            q = max(0,tilt-l)
        else:
            q = max(0,1.5-l)
        
        # note that the angular integral for even powers of mu gives J_(l+1)
        ktemp, ftemp = sph(self.qv, nu= l+(power%2), q=q)(Fq*self.renorm,extrap = False)
        
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
        za, b1, b2, Aloop, offset_loops = 0, 0, 0, 0, 0
        
        #l indep functions
        fza = 1.
        fb1 = -2 * k * self.Ulin
        #fb2sq = (self.corlin**2/2.)
        #fb1b2 = (-2*self.qv*self.Ulin*self.corlin/self.Ylin)
        foffset_loops = -0.5*ksq*2./3 * (self.xi0loop013 + self.xi0loop022)
        #
        
        for l in range(self.jn):
            #l-dep functions
            fb2 = (2*l/self.Ylin/ksq - 1) * self.Ulin**2

            X13 = 2./3 * (self.xi0loop013 - self.xi0loop13 - self.xi2loop13)
            X22 = 2./3 * (self.xi0loop022 - self.xi0loop22 - self.xi2loop22)
            Y13 = 2*self.xi2loop13
            Y22 = 2*self.xi2loop22
            
            fAloop = -0.5 * ksq * ( X13 + X22 + (1 - 2*l/ksq/self.Ylin) * (Y13 + Y22) ) - foffset_loops
            
            
            #fb1sq = (self.corlin + (2*l/self.Ylin - k**2)*self.Ulin**2)
            #fb2 = ((2*l/self.Ylin - k**2)*self.Ulin**2)
            
            #do integrals
            za += self.template(k,l,fza,expon,suppress,power=0,za=True,expon_za=exponm1)
            b1 += self.template(k,l,fb1,expon,suppress,power=1)
            b2 += self.template(k,l,fb2,expon,suppress,power=0)
            Aloop += self.template(k,l,fAloop,expon,suppress,power=0)
            offset_loops += self.template(k,l,foffset_loops,expon,suppress,power=0,za=True,expon_za=exponm1)
            #Xdot += self.template(k,l,fXdot,expon,expon0,suppress,j0=True,power=0)
            #Ydot += self.template(k,l,fYdot,expon,expon0,suppress,j0=True,power=1)
        
        
        return 4*np.pi*np.array([za,b1,b2,Aloop,offset_loops])

    
    
    def v_integrals(self, k):
        '''Do the \mu integrals in the pairwise velocity for various parameters for give 'k'
            (Currently for the pairwise velocity in k space.)
            (Commented out sections come from equivalent expressions for density power spectrum.)
        '''
        ksq = k**2
        expon = np.exp(-0.5*ksq * (self.XYlin - self.sigma))
        exponm1 = np.expm1(-0.5*ksq * (self.XYlin - self.sigma))
        suppress = np.exp(-0.5*ksq *self.sigma)

        # Note: collect all constant offset terms into 'offset'
        A, Aloop, W, b1, b2, b1b2, b1sq, offset_za, offset_Aloop, offset_b1_1, offset_b1_2,offset_b1_3= 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        
        #l indep functions
        fb1b2 = 2 * self.corlin * self.Udot
        
        foffset_za   = k * self.f * self.sigma #offset for k*A term
        foffset_Aloop = k*self.Xdot1loop[-1] #technically just the zero lag term but this is safer
        #foffset_corlin = k * self.f * self.sigma * self.corlin #offset for b1^2 term
        foffset_b1 = 2 * k * self.X10dot[-1]

        for l in range(self.jn):
            #l-dep functions
            fA = k * (self.Xdot + self.Ydot - self.f * self.sigma) + k * ( - 2*l/ksq/self.Ylin) * self.Ydot
            fAloop = k * (self.Xdot1loop + self.Ydot1loop) + k * ( - 2*l/ksq/self.Ylin) * self.Ydot1loop - foffset_Aloop

            fb1sqpb2 = 2 * (1 - 2*l/ksq/self.Ylin) * k * self.Ulin * self.Udot
            fb2 = self.U20dot
            fb1sq = k * self.corlin * (self.Xdot + (1 - 2*l/ksq/self.Ylin)*self.Ydot)
            fb1sq_odd = self.U11dot
            #fb1sq = k * self.corlin * (self.Xdot + (1 - 2*l/ksq/self.Ylin)*self.Ydot - self.f * self.sigma) + fb1sqpb2

            fb1 = 2 * self.Udot - 2 * ksq * self.Ulin * (self.Xdot + (1 - 2*l/ksq/self.Ylin)*self.Ydot) + 2 * self.Udot1loop
            fb1_even = 2 * k * (self.X10dot+self.Y10dot*(1 - 2*l/ksq/self.Ylin)) - foffset_b1
            
            fW  = - 2 * self.f * ksq * (self.v + (1 - 2*l/ksq/self.Ylin)*self.T112)
            
            
            #do integrals
            b1 += self.template(k,l,fb1,expon,suppress,power=1) + self.template(k,l,fb1_even,expon,suppress,power=0) + self.template(k,l,foffset_b1,expon,suppress,power=0,za=True,expon_za=exponm1)
            A += self.template(k,l,fA,expon,suppress,power=0) + self.template(k,l,foffset_za,expon,suppress,power=0,za=True,expon_za=exponm1)
            Aloop += self.template(k,l,fAloop,expon,suppress,power=0) + self.template(k,l,foffset_Aloop,expon,suppress,power=0,za=True,expon_za=exponm1)
            b1b2 += self.template(k,l,fb1b2,expon,suppress,power=1)

            b1sqpb2 = self.template(k,l,fb1sqpb2,expon,suppress,power=0)
            b2 += self.template(k,l,fb2,expon,suppress,power=1) + b1sqpb2
            b1sq += self.template(k,l,fb1sq,expon,suppress,power=0) + self.template(k,l,fb1sq_odd,expon,suppress,power=1) + b1sqpb2

            W += self.template(k,l,fW,expon,suppress,power=1)

            # these are now all comparison terms
            offset_za += self.template(k,l,fb1sq,expon,suppress,power=0) + b1sqpb2
            offset_Aloop += self.template(k,l,fb1sq_odd,expon,suppress,power=1)
            offset_b1_1 += self.template(k,l,fb1-2*self.Udot1loop,expon,suppress,power=1)
            offset_b1_2 += self.template(k,l,2*self.Udot1loop,expon,suppress,power=1)
            offset_b1_3 += self.template(k,l,fb1_even,expon,suppress,power=0) + self.template(k,l,foffset_b1,expon,suppress,power=0,za=True,expon_za=exponm1)
        
        return 4*np.pi*np.array([A, Aloop, W, b1, b2, b1b2, b1sq, offset_za, offset_Aloop, offset_b1_1,offset_b1_2,offset_b1_3])
    
    
    
    def spar_integrals(self, k):
        '''Do the \mu integrals for various parameters for give 'k' in \sigma_{||}
            (Commented out sections come from equivalent expressions for density power spectrum.)
            '''
        ksq = k**2
        expon = np.exp(-0.5*ksq * (self.XYlin - self.sigma))
        exponm1 = np.expm1(-0.5*ksq * (self.XYlin - self.sigma))
        suppress = np.exp(-0.5*ksq *self.sigma)
        
        # Note: collect all constant offset terms into 'offset'
        A, Aloop,W, b1,b2, b1sq, offset_za, offset_Aloop, offset_b1 = 0, 0, 0, 0, 0, 0, 0, 0, 0
        
        #l indep functions

        # pick the first or second offset based on which A term below
        foffset_za   = -ksq * self.f**2 * self.sigma**2 + self.f**2 * self.sigma #offset for k*A term
        foffset_za   = self.f**2 * self.sigma
        foffset_Aloop = self.Xddot1loop[-1] -ksq * self.f**2 * self.sigma**2
        foffset_b1 = 2 * self.X10ddot[-1]
        
        for l in range(self.jn):
            #l-dep functions
            fb1 = -2 * k * ( (2*self.Udot*self.Xdot + self.Ulin * self.Xdotdot) + (2*self.Udot*self.Ydot + self.Ulin * self.Ydotdot) * (1 - 2*l/ksq/self.Ylin))
            fb1_even = 2 * (self.X10ddot + self.Y10ddot * (1-2*l/ksq/self.Ylin)) - foffset_b1
            # keeping all terms present if field quadratic
            fA  = (-ksq*self.Xdot**2+self.Xdot) + (-ksq*2*self.Xdot*self.Ydot+self.Ydotdot)*(1-2*l/ksq/self.Ylin) -ksq*self.Ydot**2 * (1-4*l/ksq/self.Ylin + 4*l*(l-1)/ksq**2/self.Ylin**2) - foffset_za
            # for keeping to order P only; pick this to agree with Zvonimir's code
            fA  = (self.Xdotdot) + (self.Ydotdot)*(1-2*l/ksq/self.Ylin) - foffset_za
            fAloop = self.Xddot1loop + self.Yddot1loop*(1-2*l/ksq/self.Ylin) -ksq*self.Xdot**2 -ksq*2*self.Xdot*self.Ydot*(1-2*l/ksq/self.Ylin) - ksq*self.Ydot**2 * (1-4*l/ksq/self.Ylin + 4*l*(l-1)/ksq**2/self.Ylin**2) - foffset_Aloop
            
            fW  = - 5 * self.f**2 * k * (self.v + (1 - 2*l/ksq/self.Ylin)*self.T112)
            
            fb1sq = self.corlin*self.Xdotdot + (self.corlin*self.Ydotdot + 2*self.Udot**2)*(1-2*l/ksq/self.Ylin)
            fb2 = 2 * self.Udot**2 * (1-2*l/ksq/self.Ylin)
            
            #do integrals
            b1 += self.template(k,l,fb1,expon,suppress,power=1) + self.template(k,l,fb1_even,expon,suppress,power=0)
            A     += self.template(k,l,fA,expon,suppress,power=0)
            Aloop += self.template(k,l,fAloop,expon,suppress,power=0)
            W += self.template(k,l,fW,expon,suppress,power=1)
            b1sq += self.template(k,l,fb1sq,expon,suppress,power=0)
            b2 += self.template(k,l,fb2,expon,suppress,power=0)
            offset_za += self.template(k,l,foffset_za,expon,suppress,power=0,za=True,expon_za=exponm1)
            offset_Aloop += self.template(k,l,foffset_Aloop,expon,suppress,power=0,za=True,expon_za=exponm1)
            offset_b1 += self.template(k,l,foffset_b1,expon,suppress,power=0,za=True,expon_za=exponm1)
        
        
        return 4*np.pi*np.array([A,Aloop,W,b1, b1sq,b2,offset_za,offset_Aloop,offset_b1])
    
    
    def strace_integrals(self, k):
        '''Do the \mu integrals for various parameters for give 'k' in Tr(\sigma)
            (Commented out sections come from equivalent expressions for density power spectrum.)
            '''
        ksq = k**2
        expon = np.exp(-0.5*ksq * (self.XYlin - self.sigma))
        exponm1 = np.expm1(-0.5*ksq * (self.XYlin - self.sigma))
        suppress = np.exp(-0.5*ksq *self.sigma)
        
        # Note: collect all constant offset terms into 'offset'
        A, Aloop,W, b1,b2, b1sq, offset_za, offset_loops, offset_b1 = 0, 0, 0, 0, 0, 0, 0, 0, 0
        
        #l indep functions
        # pick the first or second offset based on which A term below
        foffset_za   = -ksq * self.f**2 * self.sigma**2 + 3 * self.f**2 * self.sigma #offset for k*A term
        foffset_za = 3 * self.f**2 * self.sigma
        
        #note: technically the second piece is just the sum of the zero lag pieces but
        # I am not coding it in for convenience
        foffset_loops = -ksq * self.f**2 * self.sigma**2 + 3*self.Xddot1loop[-1]
        foffset_b1 = 6 * self.X10ddot[-1]
        
        
        for l in range(self.jn):
            #l-dep functions
            fb1 = -2 * k * ( 2*(self.Xdot+self.Ydot)*self.Udot + self.Ulin*(self.Xdotdot+self.Ydotdot))
            fb1_even = 2 * (3*self.X10ddot + self.Y10ddot) - foffset_b1
            
            # keeping all terms present if field quadratic
            fA  = -ksq*self.Xdot**2 + 3*self.Xdotdot + self.Ydotdot - ksq * (self.Ydot**2+2*self.Xdot*self.Ydot)*(1-2*l/ksq/self.Ylin) - foffset_za
            # for keeping to order P only; pick this to agree with Zvonimir's code
            fA  =  3*self.Xdotdot + self.Ydotdot - foffset_za #if we keep only to order P
            fAloop = 3*self.Xddot1loop + self.Yddot1loop - ksq*(self.Xdot+self.Ydot)**2 + 2*l*(self.Ydot**2+2*self.Ydot*self.Xdot)/self.Ylin - foffset_loops
            
            fb1sq = self.corlin*(3*self.Xdotdot+self.Ydotdot)+2*self.Udot**2
            fb2 = 2 * self.Udot**2
            
            fW = - self.f**2 * k * (18 * self.v1 + 7 * self.v3 + 5 * self.T112)
                            
            #do integrals
            b1 += self.template(k,l,fb1,expon,suppress,power=1) + self.template(k,l,fb1_even,expon,suppress,power=0)
            A  += self.template(k,l,fA,expon,suppress,power=0)
            Aloop += self.template(k,l,fAloop,expon,suppress,power=0)
            W += self.template(k,l,fW,expon,suppress,power=0)
            b1sq += self.template(k,l,fb1sq,expon,suppress,power=0)
            b2 += self.template(k,l,fb2,expon,suppress,power=0)
            offset_za += self.template(k,l,foffset_za,expon,suppress,power=0,za=True,expon_za=exponm1)
            offset_loops += self.template(k,l,foffset_loops,expon,suppress,power=0,za=True,expon_za=exponm1)
            offset_b1 += self.template(k,l,foffset_b1,expon,suppress,power=0,za=True,expon_za=exponm1)
        
                                    
        return 4*np.pi*np.array([A,Aloop,W, b1, b1sq,b2,offset_za,offset_loops,offset_b1])


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



    def make_vtable(self, kmin = 1e-3, kmax = 3, nk = 100, ks=None,silent=True):
        '''Make a table of different terms of P(k) between a given
        'kmin', 'kmax' and for 'nk' equally spaced values in log10 of k
        This is the most time consuming part of the code.
        '''
        if ks is not None:
            kv = ks
        else:
            kv = np.logspace(np.log10(kmin), np.log10(kmax), nk)
        
        self.vktable = np.zeros([len(kv), self.num_velocity_components+1]) # one column for ks
        self.vktable[:, 0] = kv[:]
        for foo in range(len(kv)):
            if not silent:
                print(foo)
            self.vktable[foo, 1:] = self.v_integrals(kv[foo])

    def make_spartable(self, kmin = 1e-3, kmax = 3, nk = 100, ks=None,silent=True):
        '''Make a table of different terms of P(k) between a given
            'kmin', 'kmax' and for 'nk' equally spaced values in log10 of k
            This is the most time consuming part of the code.
        '''
        if ks is not None:
            kv = ks
        else:
            kv = np.logspace(np.log10(kmin), np.log10(kmax), nk)

        self.sparktable = np.zeros([len(kv), self.num_spar_components+1]) # one column for ks
        self.sparktable[:, 0] = kv[:]
        for foo in range(len(kv)):
            if not silent:
                print(foo)
            self.sparktable[foo, 1:] = self.spar_integrals(kv[foo])

    def make_stracetable(self, kmin = 1e-3, kmax = 3, nk = 100, ks=None):
        '''Make a table of different terms of P(k) between a given
            'kmin', 'kmax' and for 'nk' equally spaced values in log10 of k
            This is the most time consuming part of the code.
        '''
        if ks is not None:
            kv = ks
        else:
            kv = np.logspace(np.log10(kmin), np.log10(kmax), nk)
        
        self.stracektable = np.zeros([len(kv), self.num_strace_components+1]) # one column for ks
        self.stracektable[:, 0] = kv[:]
        for foo in range(len(kv)):
            self.stracektable[foo, 1:] = self.strace_integrals(kv[foo])




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


def loglog_abs(x,y,color='k',label=None):
    '''Routine to make loglog plots where negative values are indicated by dashed lines.
        Mostly for my own convenience.'''
    
    plt.loglog(x,y,color, label=label)
    plt.loglog(x,-y,color+'--')
    return None


if __name__ == '__main__':
    pks = np.loadtxt('pklin.test')
    f=1

    velda = VelocityMoments(pks[:,0],pks[:,1],f=f)
    #velda.make_ptable()
    velda.make_vtable()
    #velda.make_spartable()
    #velda.make_stracetable()

    #Let's make some plots:
    #velocity moments (pardon the dirty notation for now):
    #order of variables: A, b1, b1b2, b1sqpb2, b1sq, offset_za, offset_corlin
    ks = velda.vktable[:,0]
    ZA = velda.vktable[:,1] + velda.vktable[:,6]
    b1 = velda.vktable[:,2]
    b1sq = velda.vktable[:,4] + velda.vktable[:,5] + velda.vktable[:,7]
    b2 = velda.vktable[:,4]
    b1b2 = velda.vktable[:,3]

    loglog_abs(ks,ZA,label='ZA')
    loglog_abs(ks,b1,color='b',label='b1')
    loglog_abs(ks,b1sq,color='r',label='b1sq')
    loglog_abs(ks,b2,color='g',label='b2')
    loglog_abs(ks,b1b2,color='m',label='b1b2')

    plt.legend()
    plt.show()


