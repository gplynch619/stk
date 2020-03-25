import stk.util as util
import numpy as np

class PowerSpectrum:

    def __init__(self, infile):
        self.infile_name = infile
        kt,pkt = np.genfromtxt(infile, comments="#", dtype="float64", unpack=True, usecols=(0,1))
        self.k = kt
        self.pk = pkt
        self.interp = util.log_interpolate(self.k, self.pk)
