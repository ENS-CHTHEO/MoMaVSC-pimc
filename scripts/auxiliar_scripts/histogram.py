#!/bin/python3
import numpy as np
import argparse as ap


#Get from input
parser=ap.ArgumentParser()
parser.add_argument('--nbin',help='Input bin number. If not, 100 bins are used.',type=int)
parser.add_argument('--minfreq',help='Name of the data to histogram',type=float)
parser.add_argument('--maxfreq',help='Name of the data to histogram',type=float)
args=parser.parse_args()
nbins=args.nbin
maxf=args.maxfreq
minf=args.minfreq

#FOR FF 
#Read freqs
data=np.loadtxt('eigenvalues_ff.data')
#Generate histogram data.
values,bin_index=np.histogram(data,bins=nbins,range=(minf,maxf))

#Normalize values dividing by the oscillator number.
nosc=np.size(data)
values=values/nosc

#Write the data in the histogram file
output=open('dos_ff.data','w')
bin_index=bin_index[:-1]
for i in range(nbins):
    text=str(bin_index[i])+'  '+str(values[i])+'\n'
    output.write(text)
output.close()




#Generate now weighted spectrum by the transition dipole momentum.
analyzed_corf=np.loadtxt('freq_ips_entr_cav_dipole.data')
dipole=analyzed_corf[:,4] #Load the data corresponding to the dipole

val_spectr,bin_spectr=np.histogram(data,bins=nbins,range=(minf,maxf),weights=dipole)
val_spectr=val_spectr/nosc #Normalize.

output=open('ff_spectr.data','w')
bin_index=bin_index[:-1]
for i in range(nbins):
    text=str(bin_spectr[i])+'  '+str(val_spectr[i])+'\n'
    output.write(text)
output.close()



#FOR THE DD CASE WITH FX. ALL COMMENTED OUT BECAUSE OF POSSIBLE PROBLEMS WITH BIDIRECTIONAL DISPLACEMENTS
data=np.loadtxt('eigenvalues_dd.data')
#Generate histogram data.
values,bin_index=np.histogram(data,bins=nbins,range=(minf,maxf))

#Normalize values dividing by the oscillator number.
nosc=np.size(data)
values=values/nosc

#Write the data in the histogram file
output=open('dos_dd.data','w')
bin_index=bin_index[:-1]
for i in range(nbins):
    text=str(bin_index[i])+'  '+str(values[i])+'\n'
    output.write(text)
output.close()




#Generate now weighted spectrum by the transition dipole momentum.
#For the dd case
analyzed_corf=np.loadtxt('dd_freq_ips_entr_cav_dipole.data')
dipole=analyzed_corf[:,4] #Load the data corresponding to the dipole

val_spectr,bin_spectr=np.histogram(data,bins=nbins,range=(minf,maxf),weights=dipole)
val_spectr=val_spectr/nosc #Normalize.

output=open('dd_spectr.data','w')
bin_index=bin_index[:-1]
for i in range(nbins):
    text=str(bin_spectr[i])+'  '+str(val_spectr[i])+'\n'
    output.write(text)
output.close()



#FOR THE DD CASE WITH EQUIPARTITION
data=np.loadtxt('eigenvalues_dd_equi.data')
#Generate histogram data.
values,bin_index=np.histogram(data,bins=nbins,range=(minf,maxf))

#Normalize values dividing by the oscillator number.
nosc=np.size(data)
values=values/nosc

#Write the data in the histogram file
output=open('dos_dd_equi.data','w')
bin_index=bin_index[:-1]
for i in range(nbins):
    text=str(bin_index[i])+'  '+str(values[i])+'\n'
    output.write(text)
output.close()




#Generate now weighted spectrum by the transition dipole momentum.
#For the dd case
analyzed_corf=np.loadtxt('dd_equi_freq_ips_entr_cav_dipole.data')
dipole=analyzed_corf[:,4] #Load the data corresponding to the dipole

val_spectr,bin_spectr=np.histogram(data,bins=nbins,range=(minf,maxf),weights=dipole)
val_spectr=val_spectr/nosc #Normalize.

output=open('dd_equi_spectr.data','w')
bin_index=bin_index[:-1]
for i in range(nbins):
    text=str(bin_spectr[i])+'  '+str(val_spectr[i])+'\n'
    output.write(text)
output.close()


#Generate the IPR histogram as I do with the IR spectrum
ipr=analyzed_corf[:,1]
val_ipr,bin_ipr=np.histogram(analyzed_corf[:,0],bins=nbins,range=(minf,maxf),weights=ipr)
val_ipr=val_ipr/nosc

output=open('dd_equi_ipr.data','w')
bin_index=bin_index[:-1]
for i in range(nbins):
    text=str(bin_ipr[i])+'  '+str(val_ipr[i])+'\n'
    output.write(text)
output.close()
