#!/bin/python3

import numpy as np
import argparse as ap

#Get from input
parser=ap.ArgumentParser()
parser.add_argument('--up_cut',help='upper gaussian cutoff',type=float)
parser.add_argument('--low_cut',help='lower gaussian cutoff',type=float)
parser.add_argument('--avg',help='center of gaussian',type=float)
parser.add_argument('--desvest',help='std dev of the gaussian ',type=float)
parser.add_argument('--nosc',help='number of oscillators, cavity incl.',type=int)
args=parser.parse_args()

nosc=args.nosc
avg=args.avg
desvest=args.desvest
low_cut=args.low_cut
up_cut=args.up_cut


distribution=np.random.normal(loc=avg,scale=desvest,size=(nosc))

for i in range(nosc):
    if (distribution[i]>up_cut) or (distribution[i]<low_cut):
        distribution[i]=np.random.normal()*desvest+avg #Redraw another, go again

out_name='inp_freqs_'+str(nosc)+'_'+str(int(desvest))+'.data'

np.savetxt(out_name,distribution)
