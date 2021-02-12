import os
import numpy as np
from os import listdir
from os.path import isfile, join

import reading as read
import plotting as pl

# The directories to be read can be stored here by updating the dictionary dirs
# with any new runs that want to be stored
# The runs id can be chosen by the user to identify the run, and the directory name of the run
# has to be given to define the pointer
def read_dirs():

    dirs = {}
    dirs.update({'ini1':'M1152e_exp6k4_M4b'})
    dirs.update({'ini2':'M1152e_exp6k4'})
    dirs.update({'ini3':'M1152e_exp6k4_k60b'})
    dirs.update({'hel1':'F1152d2_sig1_t11_M2c_double'})
    dirs.update({'hel2':'F1152a_sig1_t11d_double'})
    dirs.update({'hel3':'F1152a_sig1'})
    dirs.update({'hel4':'F1152a_k10_sig1'})
    dirs.update({'noh1':'F1152b_sig0_t11_M4'})
    dirs.update({'noh2':'F1152a_sig0_t11b'})
    dirs.update({'ac1':'E1152e_t11_M4d_double'})
    dirs.update({'ac2':'E1152e_t11_M4a_double'})
    dirs.update({'ac3':'E1152e_t11_M4e_double'})

    return dirs


# Define the tuple run that will contain the spectra and time series (ts) dictionaries with the corresponding
# spectral values (and wave numbers and times) and the time evolution of the averaged values of the fields.
# The class run contains the list spectra_avail to list the read spectra
# The class when initialized, reads the spectra and the time series
class run():

    def __init__(self, name_run, dir0, quiet=True):
        
        if quiet:
            np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
        else:
            np.warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)
            
        print('Reading run ' + name_run + '\n')

        self.name_run = name_run
        dirs = read_dirs()
        self.dir_run = dirs.get(name_run)    
        self.spectra = read.read_spectra_runs(self.dir_run, dir0)
        keys = self.spectra.keys()
        print('Spectra computed: ')
        self.spectra_avail = [s for s in self.spectra.keys() if not s=="k"]
        self.spectra_avail = [s for s in self.spectra_avail if not s=="t"]
        self.spectra_avail = [s for s in self.spectra_avail if not 't_' in s]
        print(self.spectra_avail)
        print('\n')
        self.ts = read.read_ts(self.dir_run, dir0)

    def characterize_run(self, min_col=-2, max_col=0):

        # initial energy density of magnetic field
        if 'ac' in self.name_run:
            indmax = np.argmax(self.ts.get('EEK'))
            self.OmMmax = self.ts.get('EEK')[indmax]
        else: 
            indmax = np.argmax(self.ts.get('EEM'))
            self.OmMmax = self.ts.get('EEM')[indmax]
        
        self.tini = self.ts.get('t')[indmax]
   	
	# assign a color as a function of the value of \OmM 
        
        self.color = pl.pseudocolor(np.log10(self.OmMmax), max_col, min_col)
        self.compute_max_spectra()
        # initial position of spectral peak
        if 'ac' in self.name_run:
            kf2 = np.interp(self.ts.get('t'), self.spectra.get('t_kin'),
                     self.spectra.get('kin_kpeak'))
        else:
            kf2 = np.interp(self.ts.get('t'), self.spectra.get('t_mag'),
                     self.spectra.get('mag_kpeak'))
            
        self.kf = kf2[indmax]
        # Alfven speed of initial magnetic field
        self.vA = np.sqrt(2*self.OmMmax)
        # Eddy turnover time of initial magnetic field
        self.te = 1/self.kf/self.vA

    def compute_max_spectra(self):
    
        k = self.spectra.get('k')
    
        for m in self.spectra_avail:
        
            t = self.spectra.get('t_'+ m)
            Nt = len(t)
            kpeak = np.zeros(Nt)
            Emax = np.zeros(Nt)
            E = self.spectra.get(m)
            for i in range(0, Nt):
                kpeak[i], Emax[i] = compute_kpeak(k, E[i,:], quite=True)
            
            self.spectra.update({m + '_kpeak':kpeak})
            self.spectra.update({m + '_Emax':Emax})

def compute_kpeak(k, E, quite):
    
    max1 = np.argmax(E)
    # In case the spectrum is flat, compute also max of k*E (e.g. EGW)
    max2 = np.argmax(k*E)
    
    indmax = max1
    if E[max2] > E[max1]:
        indmax = max2
        
    Emax = E[indmax]
    kpeak = k[indmax]
    
    if not quite:
        print('The maximum value of the spectrum is', Emax, 'and the peak k is', kpeak)
    
    return kpeak, Emax