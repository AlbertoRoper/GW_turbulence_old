# This file contains plotting routines

from matplotlib.lines import Line2D
from colorsys import hsv_to_rgb

import matplotlib.pyplot as plt
import numpy as np

# For general use (assign colors to every run as a function of maximum value of energy)
def colors_runs(runs, diff=True):
   
    # define lines and labels for the legend
    custom_lines = []
    legs = []
    if diff:
        custom_lines_mag = []
        custom_lines_kin = []
        legs_mag = []
        legs_kin = []
    
    for i in runs:
        
        run = runs.get(i)
        custom_lines.append(Line2D([0], [0], color=run.color, lw=1))

        legs.append(run.name_run)
        
        if diff:
            if 'ac' in run.name_run:
                legs_kin.append(run.name_run)
                custom_lines_kin.append(Line2D([0], [0], color=run.color, lw=1))
            else:
                legs_mag.append(run.name_run)
                custom_lines_mag.append(Line2D([0], [0], color=run.color, lw=1))
                
    if diff:
        return legs, legs_mag, legs_kin, custom_lines, custom_lines_mag, custom_lines_kin
   
    return legs, custom_lines

def pseudocolor(val, minval, maxval):
    """ Convert val in range minval..maxval to the range 0..120 degrees which
        correspond to the colors Red and Green in the HSV colorspace.
    """
    h = (float(val-minval) / (maxval-minval)) * 120

    # Convert hsv color (h,1,1) to its rgb equivalent.
    # Note: hsv_to_rgb() function expects h to be in the range 0..1 not 0..360
    r, g, b = hsv_to_rgb(h/360, 1.,  1.)
    return r, g, b

# This routine is used to plot the averaged field values as a function of time (from the time series)
# for GW energy (time derivative of t), magnetic and velocity fields

# It can be used to also plot hc, EKmax, EMmax, rho, and total energy, if the corresponding
# logicals are set to True when calling the function

def palette_colors(min_col=-2, max_col=-1.5, min_Om=2.5e-3, max_Om=.2):
    
    plt.figure(figsize=(5,3))
    
    Om = np.logspace(np.log10(min_Om), min(-1.5, np.log10(max_Om)), 100)
    Om2 = np.logspace(min(-1.5, np.log10(max_Om)), np.log10(max_Om), 10)
    
    cols = []
    for i in range(0, len(Om)):
        col = pseudocolor(np.log10(Om[i]), max_col, min_col)
        aux = np.linspace(0, 1, 10)
        plt.plot(Om[i]*aux**0, aux, color = col, lw = 3.5)
    for j in range(0, len(Om2)):
        plt.plot(Om2[j]*aux**0, aux, color = 'red', lw = 30)
    plt.xlim([Om[0], Om2[-1]])
    plt.ylim([0, 1])
    plt.xscale('log')
    plt.xlabel('E (M or K)')
    plt.yticks([])

def EEM_EEK_EEGW_vs_t(run, lhc=False, lEKmax=False, lEMmax=False,
                     save=False, figsize=(8,5), EGW0=1., EGW1=1., EK0=1., EK1=1., t0=1., t1=1.,
                     hc0=1., hc1=1., ET0=1., ET1=1.):
    
    # The default routine plots EGW, EM, EK for the runs
    EEGW = run.ts.get('EEGW')
    EEM = run.ts.get('EEM')
    EEK = run.ts.get('EEK')
    t = run.ts.get('t')
    
    # Plot EM
    if not 'ac' in run.name_run:
        plt.figure(1, figsize=figsize)
        plt.ylabel('EM')
        plt.plot(t - 1, EEM, label = run.name_run, color=run.color)
        aux = np.logspace(np.log10(run.OmMmax/1.5), np.log10(1.5*run.OmMmax), 5)
        # aux2 = np.logspace(np.log10(t[0] + run.tini - 2), np.log10(run.te + run.tini - 1), 5)
        aux2 = np.logspace(np.log10(t[1] + run.tini - 2), np.log10(run.te + run.tini - .8), 5)
        plt.plot((run.tini + .25*run.te - 1)*aux**0, aux, ls = 'dashed', color = run.color)
        plt.plot(aux2, run.OmMmax*aux2**0, ls = 'dashed', color = run.color)
    if t0 != 1. or t1 != 1.:
        plt.xlim([t0, t1])
    if save:
        plt.savefig('EM_vs_t_runs.pdf')
        
    # Plot EK
    if 'ac' in run.name_run:
        plt.figure(2, figsize=figsize)
        plt.ylabel('EK')
        plt.plot(t[1:] - 1, EEK[1:], label = run.name_run, color = run.color)
        if EK0 != 1 or EK1 != 1:
            plt.ylim([EK0, EK1])
        if t0 != 1. or t1 != 1.:
            plt.xlim([t0, t1])
        if save:
            plt.savefig('EK_vs_t_runs.pdf')
        
    # Plot EGW
    plt.figure(3, figsize=figsize)
    plt.ylabel('EGW')
    plt.plot(t[1:] - 1, EEGW[1:], label = run.name_run, color=run.color)
    if EGW0 == 1 and EGW1 == 1:
        ymin, ymax = plt.ylim()
    else:
        ymin = EGW0
        ymax = EGW1
    aux = np.logspace(np.log10(ymin), np.log10(ymax), 10)
    plt.plot(run.tini - 1 + 1/run.kf*aux**0, aux, color = 'black', ls = 'dashed')
    if t0 != 1. or t1 != 1.:
        plt.xlim([t0, t1])
    if save:
        plt.savefig('EGW_vs_t_runs.pdf')
    plt.ylim([ymin, ymax])
    
    i = 4
    if lhc:
        # compute and plot hc
        hc  = run.ts.get('hrms')
        plt.figure(i, figsize=figsize)
        plt.ylabel('hc')
        plt.plot(t[1:] - 1, hc[1:], label = run.name_run, color = run.color)
        if hc0 != 1 or hc1 == 1:
            plt.ylim([hc0, hc1])
        if t0 != 1. or t1 != 1.:
            plt.xlim([t0, t1])
        if save:
            plt.savefig('hc_vs_t_runs.pdf')
        i += 1
        
    if lEMmax:
        if not 'ac' in run.name_run:
            # compute max value of EEM and plot
            EEMmax = run.ts.get('bmax')**2/2
            plt.figure(i, figsize=figsize)
            plt.ylabel('max EM')
            plt.plot(t[1:] - 1, EEMmax[1:], label = run.name_run, color = run.color)
            if t0 != 1. or t1 != 1.:
                plt.xlim([t0, t1])
            if save:
                plt.savefig('EMmax_vs_t_runs.pdf')
        i += 1
        
    if lEKmax:
        if 'ac' in run.name_run:
            # compute max value of EEK and plot
            EEKmax = run.ts.get('umax')**2/2
            plt.figure(i, figsize=figsize)
            plt.ylabel('max EK')
            plt.plot(t[1:] - 1, EEKmax[1:], label = run.name_run, color = run.color)
            if t0 != 1. or t1 != 1.:
                plt.xlim([t0, t1])
            if save:
                plt.savefig('EKmax_vs_t_runs.pdf')
        i += 1

    i -= 1
    
    for j in range(1, i + 1):
        plt.figure(j, figsize=figsize)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('dt = t - tini')
        
    return i
