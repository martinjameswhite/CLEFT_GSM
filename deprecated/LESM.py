# Python wrapper class for LSM.
from __future__ import print_function,division



import commonstuff as C
import numpy       as np
import socket
import ctypes
import os



class LSM:
    """
    LSM:
    A Python class to conviently "wrap" calls to the
    Lagrangian streaming model (LSM) code.
    This uses temporary files to pass the information around.
    """
    __author__ = "Martin White"
    __version__ = "1.0"
    __email__  = "mwhite@berkeley.edu"
    #
    def __call__(self,pkfile,ff,b1,b2,bs,Aeft,Aeft1,s2FoG,Apar,Aperp):
        """    
        Runs the code and returns a NumPy array containing
        the data returned by the code.
        Note: as currently configured returns s^2 xi_ell, not  xi_ell.
        """
        Aeft2 = 0.0
        os.environ['OMP_NUM_THREADS']=str(C.Nthread)
        ret = self.mylib.call_lesm(ctypes.c_char_p(pkfile),\
          ctypes.c_double(ff),ctypes.c_double(b1),ctypes.c_double(b2),\
          ctypes.c_double(bs),ctypes.c_double(Aeft),\
          ctypes.c_double(Aeft1),ctypes.c_double(Aeft2),\
          ctypes.c_double(s2FoG),ctypes.c_double(Apar),ctypes.c_double(Aperp),\
          ctypes.c_char_p(self.tmpfn))
        if (ret==0)&(os.path.isfile(self.tmpfn)):
            dd = np.loadtxt(self.tmpfn)
            os.remove(self.tmpfn)
        else:
            outstr = "LESM call failed with: "+pkfile+","+str(ff)+","+str(b1)+\
                     ","+str(b2)+","+str(bs)+","+str(Aeft)+","+str(Aeft1)+","+\
                     str(s2FoG)+","+str(Apar)+","+str(Aperp)
            raise RuntimeError,outstr
            dd = None
        return(dd)
        #
    def __init__(self):
        """
        __init__(self):
        """
        # Basic initialization, including a temporary file
        # whose name is based on the current host, PPID and PID.
        self.tmpfn = "lesm_%s_%d_%d.txt"%\
          (socket.gethostname(),os.getppid(),os.getpid())
        self.mylib = ctypes.CDLL(C.basedir+"/sm//lesm_ctypes.so")
    #




def peak_background_bias(nu):
    """
    peak_background_bias(nu):
    Returns the Lagrangian biases, (b1,b2), given nu.
    This is helpful if we want to make our basis set f, nu and sFog.
    """
    delc = 1.686
    a    = 0.707
    p    = 0.30
    anu2 = a*nu**2
    b1   = (anu2-1+2*p/(1+anu2**p))/delc
    b2   = (anu2**2-3*anu2+2*p*(2*anu2+2*p-1)/(1+anu2**p))/delc**2
    return( (b1,b2) )
    #
