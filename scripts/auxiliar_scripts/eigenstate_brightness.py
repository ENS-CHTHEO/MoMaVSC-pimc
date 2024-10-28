#!/bin/python3
import numpy as np
import scipy.linalg as lin
import argparse as ap
import os

parser=ap.ArgumentParser()
parser.add_argument('--nbin',help='Number of bins for the histogram',type=int)
args=parser.parse_args()
nbin=args.nbin

data=np.loadtxt('dd_equi_freq_ips_entr_cav_dipole.data')
ir_spectr_hist,ir_spectr_bins=np.histogram(data[:,0],weights=data[:,4],bins=nbin,range=(2800,4000))


ir_spectr_hist=ir_spectr_hist/np.size(data[1:,0]) #Normalize by the number of matter modes

forsaving=np.vstack((ir_spectr_bins[1:],ir_spectr_hist)).T



darkid=np.argwhere(ir_spectr_hist>0)[1:-1] #Drop the indices of the polaritons, first and last nonzero bins

np.savetxt('dark_state_brightnesses.data',forsaving)
average_light_dark=np.average(ir_spectr_hist[darkid])


with open('av_light_in_dark.data','w') as f:
    f.write(str(average_light_dark))

