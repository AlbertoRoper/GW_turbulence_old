# This file contains plotting routines

from matplotlib.lines import Line2D
from colorsys import hsv_to_rgb

import matplotlib.pyplot as plt
import numpy as np

# function that defines the legends for the specific runs of the paper
# from https://arxiv.org/pdf/2011.05556.pdf
def legends_PRR(runs):
    
    # define lines and labels for the legend
    custom_lines = []
    legs = []
    custom_lines_mag = []
    legs_mag = []
    custom_lines_kin = []
    legs_kin = []
    
    for i in runs:
        
        run = runs.get(i)
        legs.append(run.name_run)
        
        if 'K' in run.name_run:
            run.color = (run.color[0], 0, run.color[2])
            legs_kin.append(run.name_run)
            custom_lines.append(Line2D([0], [0], color=run.color, lw=1))
            custom_lines_kin.append(Line2D([0], [0], color=run.color, lw=1))
        
        else:
            run.color = (0, run.color[1], run.color[2])
            legs_mag.append(run.name_run)
            custom_lines.append(Line2D([0], [0], color=run.color, lw=1))
            custom_lines_mag.append(Line2D([0], [0], color=run.color, lw=1))
   
    return legs, custom_lines, legs_kin, legs_mag, custom_lines_kin, custom_lines_mag

# function that defines the legends for the specific runs of the paper
# from https://arxiv.org/pdf/1903.08585.pdf
def legends_PRD(runs):
   
    # define lines and labels for the legend
    custom_lines = []
    legs = []
    custom_lines_mag = []
    custom_lines_kin = []
    legs_mag = []
    legs_kin = []
    
    for i in runs:
        
        run = runs.get(i)
        custom_lines.append(Line2D([0], [0], color=run.color, lw=1))

        legs.append(run.name_run)
        
        if 'ac' in run.name_run:
            legs_kin.append(run.name_run)
            custom_lines_kin.append(Line2D([0], [0], color=run.color, lw=1))
        else:
            legs_mag.append(run.name_run)
            custom_lines_mag.append(Line2D([0], [0], color=run.color, lw=1))
            
    return legs, legs_mag, legs_kin, custom_lines, custom_lines_mag, custom_lines_kin

# return legends for general use
def legends(runs):
    
    # define lines and labels for the legend
    custom_lines = []
    legs = []
    
    for i in runs:
        
        run = runs.get(i)
        custom_lines.append(Line2D([0], [0], color=run.color, lw=1))
        legs.append(run.name_run)
   
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

# This routine is used to plot the averaged field values as a function of time (from the time series)
# for GW energy (time derivative of t), magnetic and velocity fields

# It can be used to also plot hc, EKmax, EMmax, rho, and total energy, if the corresponding
# logicals are set to True when calling the function   
 
import numpy as np

def EEM_EEK_EEGW_vs_t(run, lhc=False, lEKmax=False, lEMmax=False,
                     save=False, figsize=(8,5), EGW0=0, EGW1=0, EK0=0, EK1=0, t0=0, t1=0,
                     hc0=0, hc1=0, ET0=0, ET1=0, EM0=0, EM1=0,
                     EMmax0=0, EMmax1=0, EKmax0=0, EKmax1=0, lM=True, lK=True, lGW=True):
    
    # The default routine plots EGW, EM, EK for the runs
    t = run.ts.get('t')
    
    i = 1
    # Plot EM
    if lM:
        EEM = run.ts.get('EEM')
        plt.figure(2, figsize=figsize)
        plt.ylabel('EM')
        plt.plot(t - 1, EEM, label = run.name_run, color=run.color)
        aux = np.logspace(np.log10(run.OmMmax/1.5), np.log10(1.5*run.OmMmax), 5)
        aux2 = np.logspace(np.log10(t[1] + run.tini - 2), np.log10(run.te + run.tini - .8), 5)
        plt.plot((run.tini + .25*run.te - 1)*aux**0, aux, ls = 'dashed', color = run.color)
        plt.plot(aux2, run.OmMmax*aux2**0, ls = 'dashed', color = run.color)
        adjust_ax_lims(EM0, EM1, t0, t1)
        if save == True: plt.savefig('EM_vs_t_runs.pdf')
        i += 1
        
    # Plot EK
    if lK:
        EEK = run.ts.get('EEK')
        plt.figure(3, figsize=figsize)
        plt.ylabel('EK')
        plt.plot(t[1:] - 1, EEK[1:], label = run.name_run, color = run.color)
        adjust_ax_lims(EK0, EK1, t0, t1)
        if save == True: plt.savefig('EK_vs_t_runs.pdf')
        i += 1
        
    # Plot EGW
    if lGW:
        EEGW = run.ts.get('EEGW')
        plt.figure(1, figsize=figsize)
        plt.ylabel('EGW')
        plt.plot(t[1:] - 1, EEGW[1:], label = run.name_run, color=run.color)
        adjust_ax_lims(EGW0, EGW1, t0, t1)
        ymin, ymax = plt.ylim()
        aux = np.logspace(np.log10(ymin), np.log10(ymax), 10)
        plt.plot(run.tini - 1 + 1/run.kf*aux**0, aux, color = 'black', ls = 'dashed')
        if save == True: plt.savefig('EGW_vs_t_runs.pdf')
        i += 1

    if lhc:
        if lGW:
        # compute and plot hc
            hc  = run.ts.get('hrms')
            plt.figure(4, figsize=figsize)
            plt.ylabel('hc')
            plt.plot(t[1:] - 1, hc[1:], label = run.name_run, color = run.color)
            adjust_ax_lims(hc0, hc1, t0, t1)
            if save == True: plt.savefig('hc_vs_t_runs.pdf')
            i += 1
        
    if lEMmax:
        if lM:
            # compute max value of EEM and plot
            EEMmax = run.ts.get('bmax')**2/2
            plt.figure(5, figsize=figsize)
            plt.ylabel('max EM')
            plt.plot(t[1:] - 1, EEMmax[1:], label = run.name_run, color = run.color)
            adjust_ax_lims(EMmax0, EMmax1, t0, t1)
            if save == True: plt.savefig('EMmax_vs_t_runs.pdf')
            i += 1
        
    if lEKmax:
        if lK:
            # compute max value of EEK and plot
            EEKmax = run.ts.get('umax')**2/2
            plt.figure(6, figsize=figsize)
            plt.ylabel('max EK')
            plt.plot(t[1:] - 1, EEKmax[1:], label = run.name_run, color = run.color)
            adjust_ax_lims(EKmax0, EKmax1, t0, t1)
            if save == True: plt.savefig('EMmax_vs_t_runs.pdf')
            i += 1
    i += 1
    
    for j in range(1, i + 1):
        plt.figure(j, figsize=figsize)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('dt = t - tini')
        
    return i

def adjust_ax_lims(E0, E1, t0, t1):

    ymin, ymax = plt.ylim()
    if E0 != 0: ymin = E0
    if E1 != 0: ymax = E1
    plt.ylim([ymin, ymax])
    tmin, tmax = plt.xlim()
    if t0 != 0: tmin = t0
    if t1 != 0: tmax = t1
    t0, t1 = plt.xlim([tmin, tmax])
    