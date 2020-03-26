import stk.util as util
import numpy as np

class PowerSpectrum:

    def __init__(self, infile=None):
        self.infile_name = infile
        self.k = []
        self.pk = []
        self.interp = ""
        if infile is not None:
            self.init_from_file(infile)

    def init_from_file(self, infile):
        self.infile_name = infile
        kt,pkt = np.genfromtxt(infile, comments="#", dtype="float64", unpack=True, usecols=(0,1))
        self.k = kt
        self.pk = pkt
        self.interp = util.log_interpolate(self.k, self.pk)

    def init_from_arrays(self, in_k, in_pk):
        self.k = in_k
        self.pk = in_pk
        self.interp = util.log_interpolate(self.k, self.pk)

class TransferFunction:

    def __init__(self, infile, s8, ns, renorm = True):
        self.infile_name = infile
        k, tm = np.genfromtxt(infile, comments="#", dtype="float64", unpack=True, usecols=(0,6))
        self.k = k
        self.tm=tm
        self.s2 = util.ComputeSigma2(self.k, self.tm, ns)
        snorm = s8**2/self.s2
        if(renorm):
            self.power = self.compute_pk(snorm, ns) 
        else:
            self.power = self.compute_pk(1.0, ns)
    def compute_pk(self, snorm, ns):
        pk = snorm*self.k**ns*self.tm**2
        temp = PowerSpectrum()
        temp.init_from_arrays(self.k, pk)
        return temp 
