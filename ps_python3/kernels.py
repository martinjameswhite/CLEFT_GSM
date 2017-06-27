import numpy
from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.integrate import trapz, simps
## Define functions

#constants

class Q:

    def __init__(self, k, p, ilpk):
        
        self.kp = k
        self.p = p
        self.ilpk = ilpk
        self.kint = numpy.logspace(-6, 3, 1e3)
        self.pkint = self.ilpk(self.kint)

        self.tol = 10**-5. 
        self.glxval, self.glwval = numpy.loadtxt("gl_128.txt", unpack = True)
        self.glxval = self.glxval.reshape(1, -1)
        self.glwval = self.glwval.reshape(1, -1)

        #self.listQapprox = [self.Q1approx, self.Q2approx, self.Q3approx, 0, self.Q5approx, 0, 0, self.Q8approx, self.Qs2approx]
        #self.listQ = [self.Q1, self.Q2, self.Q3, 0, self.Q5, 0, 0, self.Q8, self.Qs2]
        self.dictQapprox = {1:self.Q1approx, 2:self.Q2approx, 3:self.Q3approx, 5:self.Q5approx, 8:self.Q8approx, -1:self.Qs2approx}
        self.dictQ = {1:self.Q1, 2:self.Q2, 3:self.Q3, 5:self.Q5, 8: self.Q8, -1:self.Qs2}

    #Define Q:
    def Q1approx(self, r, x):
        return 0.25*(1 + x)**2
        
    def Q1(self, r, x):
        y =  1 + r**2 - 2*r*x
        return (r**2 *(1 - x**2)**2)/y**2

    def Q2approx(self, r, x):
        return 0.25*x*(1 + x)
        
    def Q2(self, r, x):
        y =  1 + r**2 - 2*r*x
        return ((1- x**2)*r*x*(1 - r*x))/y**2

    def Q3approx(self, r, x):
        return 0.25*x**2
        
    def Q3(self, r, x):
        y =  1 + r**2 - 2*r*x
        return (x**2 *(1 - r*x)**2)/y**2

    def Q5approx(self, r, x):
        y =  1 + r**2 - 2*r*x
        return x*(1 + x)/2.

    def Q5(self, r, x):
        y =  1 + r**2 - 2*r*x
        return r*x*(1 - x**2.)/y

    def Q8approx(self, r, x):
        y =  1 + r**2 - 2*r*x
        return (1 + x)/2.

    def Q8(self, r, x):
        y =  1 + r**2 - 2*r*x
        return r**2.*(1 - x**2.)/y

    def Qs2approx(self, r, x):
        y =  1 + r**2 - 2*r*x
        return (1 + x)*(1- 3*x)/4.

    def Qs2(self, r, x):
        y =  1 + r**2 - 2*r*x
        return r**2.*(x**2 - 1.)*(1 - 2*r**2 + 4*r*x - 3*x**2)/y**2


    def Q_internal(self, k, n, r, approx = False):
        
        #if approx:
        #    func = lambda r, x: self.ilpk(k*numpy.sqrt(1 + r**2 - 2*r*x))*self.listQapprox[n-1](r, x)
        #else:
        #    func = lambda r, x: self.ilpk(k*numpy.sqrt(1 + r**2 - 2*r*x))*self.listQ[n-1](r, x)
        if approx:
            func = lambda r, x: self.ilpk(k*numpy.sqrt(1 + r**2 - 2*r*x))*self.dictQapprox[n](r, x)
        else:
            func = lambda r, x: self.ilpk(k*numpy.sqrt(1 + r**2 - 2*r*x))*self.dictQ[n](r, x)
        return (self.glwval*func(r, self.glxval)).sum(axis = -1)

    def Q_external(self, k, n):
        fac = k**3 /(2*numpy.pi)**2
        r = (self.kint/k)
        absr1 = abs(r-1)
        mask = absr1 < self.tol
        y = numpy.zeros_like(r)
        y[~mask]  = self.Q_internal(k, n, r[~mask].reshape(-1, 1))
        if mask.sum():
            y[mask] = self.Q_internal(k, n,  r[mask].reshape(-1, 1), approx = True)
        y *= self.pkint
        return fac*trapz(y, r)
        #return fac*simps(y, r)


class R:

    def __init__(self, k, p, ilpk):
        
        self.kp = k
        self.p = p
        self.ilpk = ilpk
        self.kint = numpy.logspace(-6, 3, 1e3)
        self.pkint = self.ilpk(self.kint)

        self.tol = 10**-5. 
        self.glxval, self.glwval = numpy.loadtxt("gl_128.txt", unpack = True)
        self.glxval = self.glxval.reshape(1, -1)
        self.glwval = self.glwval.reshape(1, -1)

        self.listRapprox = [self.R1approx, self.R2approx]
        self.listR = [self.R1, self.R2]

        #TODO:See if this works
        #self.rvals = numpy.logspace(-6, 5, 10**6)
        #self.integral1 = self.dobefore(self.rvals, 1)
        #self.integral2 = self.dobefore(self.rvals, 2)

    def R1approx(self, r, x):
        return 0.5*(1- x)*(1 + x)**2
        
    def R1(self, r, x):
        y =  1 + r**2 - 2*r*x
        return (r**2 *(1 - x**2)**2)/y

    def R2approx(self, r, x):
        return 0.5*x*(1 + x)*(1-x)
        
    def R2(self, r, x):
        y =  1 + r**2 - 2*r*x
        return ((1- x**2)*r*x*(1 - r*x))/y


    def R_internal(self, n, r, approx = False):
        
        if approx:
            func = lambda r, x: self.listRapprox[n-1](r, x)
        else:
            func = lambda r, x: self.listR[n-1](r, x)
        return (self.glwval*func(r, self.glxval)).sum(axis = -1)

    def R_external(self, k, n):

        fac = k**3 /(2*numpy.pi)**2 *self.ilpk(k)
        r = (self.kint/k)
        absr1 = abs(r-1)
        mask = absr1 < self.tol
        y = numpy.zeros_like(r)
        y[~mask]  = self.R_internal(n, r[~mask].reshape(-1, 1))
        if mask.sum():
            y[mask] = self.R_internal(n,  r[mask].reshape(-1, 1), approx = True)
        y *= self.pkint
        return fac*trapz(y, r)

    #TODO:
#    def dobefore(self, r, n):
#        y = numpy.zeros_like(r)
#        absr1 = abs(r-1)
#        mask = absr1 < self.tol
#        y[~mask]  = self.R_internal(n, r[~mask].reshape(-1, 1))
#        if mask.sum():
#            y[mask] = self.R_internal(n,  r[mask].reshape(-1, 1), approx = True)
#        return y
#
#    def R_external2(self, k, n):
#        fac = k**3 /(2*numpy.pi)**2 *self.pk(k)
#        if n ==1:
#            integral = self.integral1
#        else: 
#            integral = self.integral2
#        #More reshping here
#        kvals = self.rvals*k
#        min = numpy.where(kvals > kmin)[0][0]
#        max = numpy.where(kvals < kmax)[0][0]
#        pvals = self.pk(kvals[min:max])
#        return fac*(pvals*integral[min:max], self.rvals[min:max])
