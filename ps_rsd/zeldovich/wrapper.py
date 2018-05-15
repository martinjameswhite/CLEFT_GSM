# Python wrapper class.
from __future__ import print_function,division


import numpy  as np
import socket
import ctypes
import os



class Zeldovich:
    """
    Zeldovich:
    A Python class to conviently "wrap" calls to the Zeldovich code.
    This uses temporary files to pass the information around.
    """
    __author__ = "Martin White"
    __version__ = "1.0"
    __email__  = "mwhite@berkeley.edu"
    #
    def __call__(self,pkfile,ff):
        """    
        Runs the code and returns a NumPy array containing
        the data returned by the code.
        """
        os.environ['OMP_NUM_THREADS']=str(self.Nthread)
        ret = self.mylib.call_lesm(ctypes.c_char_p(pkfile),\
          ctypes.c_double(ff),ctypes.c_char_p(self.tmpfn))
        if (ret==0)&(os.path.isfile(self.tmpfn)):
            dd = np.loadtxt(self.tmpfn)
            os.remove(self.tmpfn)
        else:
            outstr = "Call failed with: "+pkfile+","+str(ff)
            raise RuntimeError,outstr
            dd = None
        return(dd)
        #
    def __init__(self,Nthread=1):
        """
        Basic initialization.
        """
        # Basic initialization, including a temporary file
        # whose name is based on the current host, PPID and PID.
        self.tmpfn = "zeld_%s_%d_%d.txt"%\
          (socket.gethostname(),os.getppid(),os.getpid())
        self.mylib = ctypes.CDLL(os.getcwd()+"/libzeld.so")
        self.Nthread = Nthread
    #



