
import os
import numpy as np
from os import listdir
from os.path import isfile, join

# Function that reads all the spectra files stored in the run directory
# It reads all files that start with power_ and powerhel_ (with the exeption of power_krms)
# It stores the spectra in the dictionary spectra
# To restore the spectra values, the wavenumber and the time arrays from the dictionary one should use:
# ------------------------------------
# - spectra.get('#name_spectrum'), for #name_spectrum referring to the different spectra read.
# Some examples are:
#      - GWs (spectrum from the time derivatives of the strains)
#      - GWh (spectrum from the strains)
#      - mag (spectrum from the magnetic field)
#      - kin (spectrum from the velocity field)
#      - Tpq (spectrum from the unprojected stress, computed from velocity and magnetic fields)
#      - SCL (spectrum from the scalar mode of the stress tensor)
#      - VCT (spectrum from the vector mode of the stress tensor)
#      - Str (spectrum from the projected stress TT)
#       Note that this list may vary depending on the run and could not include all of them or it could add some
#          new spectra computed, to print the resulting spectra one can print:
#          [s for s in spectra.keys()]
# - spectra.get('hel#name_spectrum'), for the helical spectra of the same fields
# - spectra.get('k') to get the wave numbers
# - spectra.get('t') to get the times of the computed spectra

def read_spectra_runs(dir_run, dir0):

    os.chdir(dir0 + dir_run + '/data/')

    # defines the list of power spectra to be read in matching and matchinghel (for helical spectra)
    onlyfiles = [f for f in listdir() if isfile(join(f))]
    matching = [s for s in onlyfiles if "power" in s]
    matching = [s for s in matching if not "krms" in s]
    matchinghel = [s for s in matching if "hel" in s]
    matchinghel = [s for s in matchinghel if not "swp" in s]
    matching = [s for s in matching if not "hel" in s]
    matching = [s for s in matching if not "swp" in s]
    
    # Read the wave number from power_krms.dat and normalize using the size of the box length L
    # (assuming a cubic domain)
    k = read_k(dir_run, dir0)
    L = read_L(dir_run, dir0)
    k = k*2*np.pi/L
    
    spectra = {}                # initialize the dictionary spectra              
    spectra.update({'k':k})     # add the wavenumber to the dictionary
    
    # read and add to the dictionary spectra all the spectra to be read from matching list
    for i in matching:
        aux = i.replace('power_', '')
        aux = aux.replace('.dat', '')
        times, sps = read_spectrum(dir_run, aux,dir0,  hel=False)
        sps = np.asarray(sps, dtype=object)
        spectra.update({aux:sps})
    # add the times array to the dictionary from any of the spectra read
        spectra.update({'t_' + aux:times})

    # read and add to the dictionary spectra all the helical spectra to be read from matchinghel list
    for i in matchinghel:
        aux = i.replace('powerhel_', '')
        aux = aux.replace('.dat', '')
        times, sps = read_spectrum(dir_run, aux, dir0, hel=True)
        sps = np.asarray(sps, dtype=object)
        spectra.update({'hel' + aux:sps})
        spectra.update({'t_hel' + aux:times})

    os.chdir(dir0)
    
    return spectra

# read time data series, which contains the averaged values of the fields as a function of time
def read_ts(dir_run, dir0):
    
    os.chdir(dir0 + dir_run + '/data/')
    
    # Read the file from time_series.dat
    af = np.loadtxt('time_series.dat')
    
    # legend.dat contains the values of the fields that are stored in time_series.dat
    # To change this one should change the print.in file before executing the run, to decide which fields
    # should be stored in the time series
    with open('legend.dat') as fp:
        leg = fp.readline()
        leg = leg.split('-')
        leg = [s for s in leg if s.isalpha()]
    
    # define the ts (time series) dictionary and update with the values read from the time series file
    Nts = len(leg)
    ts = {}
    for i in range(0, Nts):
        #mag = np.array(af[:,i], dtype='float')
        ts.update({leg[i]:(af[:,i])})
        
    os.chdir(dir0)
    
    return ts

# Reads the spectrum file in dir_runs named power_"spectrum" or powerhel_"spectrum" depending if hel (helical)
# field is True or False
def read_spectrum(dir_run, spectrum, dir0, hel=False):
    
    os.chdir(dir0 + dir_run + '/data/') # run directory storing the spectrum file
    
    # Decide if reading power or powerhel for symmetric or antisymmetric spectrum, respectively
    # and define the file to read
    # hel=False is the default
    power = 'power_'
    if hel:
        power = 'powerhel_'
    file = power + spectrum + '.dat'

    # reading the file of the spectrum data and store the values of the spectrum (specs) as a function of k,
    # for every value of time (stored in variable times)
    with open(file) as fp:
        line = fp.readline()
        times = []
        sp = [[]]
        cnt = 1
        content = line.strip()
        times.append(content)
        len_st = len(content)
        specs=[]
        while line:
            #print("Line {}: {}".format(cnt, line.strip()))
            line = fp.readline()
            content = line.strip()
            if len(content) == len_st:
                times.append(content)
                sp.append(specs)
                specs = []
            else:
                spec = content.split()
                specs.append(spec)
            cnt += 1
    sp.append(specs)

    # define the length of the times values read, and compare with the size of the 2D array spec
    nt = np.shape(sp)[0] - 1
    nt2 = len(times)
    if nt != nt2:
        print('The number of points in time does not coincide with the number of spectra functions')

    # read the array read for spec and rewrite it as a 2D array as a function of time (first index)
    # and as a function of k (second index).
    # Note that previously spec had the format of the data file (chunks of values at every time)
    sps = []
    test = False
    for l in range(0, min(nt, nt2)):
        a = np.shape(sp[l + 1])
        if len(a) == 1:
            a = np.shape(sp[l + 1])[0] - 1
            b = np.shape(sp[l + 1][0])[0]
        else:
            a, b = a
        sp0 = np.zeros(a*b)
        cnt = 0
        for i in range(0, a):
            for j in range(0, b):
                sp0[cnt] = sp[l + 1][i][j]
                cnt += 1
        sps.append(np.array(sp0, dtype='double'))
        
    # redefine the times and sps arrays as numpy arrays to return from the function
    times = np.array(times, dtype='double')
    sps = np.array(sps)

    os.chdir(dir0)
    
    return times, sps

# Read the values of wave numbers that correspond to the spectra files
# Note that the return ks array is normalized such that the smallest wave number is 1
# (independently of the size of the box)
def read_k(dir_run, dir0):
    
    os.chdir(dir0 + dir_run + '/data/')

    # The values of the wave numbers are stored in power_krms.dat
    # Read file and store in numpy array ks
    ak = np.loadtxt('power_krms.dat')
    a, b = np.shape(ak)
    ks = np.zeros(a*b)
    cnt = 0
    for i in range(0, a):
        for j in range(0, b):
            ks[cnt]=ak[i][j]
            cnt+=1        
    ks = np.array(ks, dtype='double')
            
    os.chdir(dir0)
    
    return ks

# Read the length of the domain to compute the actual wave numbers, since the values
# read in previous function are normalized
def read_L(dir_runs, dir0):
    
    os.chdir(dir0 + dir_runs + '/data/')
    
    # The length of the size domain is stored in the file param.nml (where many parameters of the
    # run are stored)
    # Note that we assume cubic domain such that the volum is Lx^3
    with open('param.nml') as fp:
        content = fp.readlines()
        content = [x.strip() for x in content] 

    matching = [s for s in content if "LXYZ" in s]
    LXYZ = matching[0].split()[1]
    L = float(LXYZ.split('*')[1])
    
    os.chdir(dir0)
    
    return L


# This function is used to write the directory (s2) and name (s1) of a new run in the 'run.py' file, by
# generating a copy of run.py in 'fp'
def add_run(s1, s2, sf):
    
    fp2 = open(sf, 'w')
    with open('run.py') as fp:
        line = fp.readline()
        fp2.write(line)
        i = 0
        while line:
            line = fp.readline()
            fp2.write(line)
            content = line.split()
            if content:
                if content[0] == 'def':
                    if content[1] == 'read_dirs():':
                        fp2.write(fp.readline())
                        fp2.write(fp.readline())
                        fp2.write('    dirs.update({\'%s'%s1)
                        fp2.write('\':\'%s\'})'%s2)
                        fp2.write('\n')

