#!/bin/python3
import numpy as np
import scipy.linalg as lin
import argparse as ap
import os

#Get from input
parser=ap.ArgumentParser()
parser.add_argument('--run_id',help='Input number of runs.',type=int)
parser.add_argument('--nosc',help='Number of oscillators',type=int)
parser.add_argument('--redmass',help='Reduced mass for diagonalization',type=float)
parser.add_argument('--temp',help='temperature',type=float)
args=parser.parse_args()
nosc=args.nosc
index=args.run_id
mass=args.redmass
temperature=args.temp

ff_corr=np.zeros([nosc,nosc])
fx_corr=np.zeros([nosc,nosc])
dd_corr=np.zeros([nosc,nosc])

#Average the correlation functions.
for idx in range(1,index+1):

    name='run_'+str(idx)
    os.chdir(name)
    print(os.getcwd())
    load_ff=np.loadtxt('force_corr.data')
    load_fx=np.loadtxt('force_displacement_corr.data')
    load_dd=np.loadtxt('disp_disp_corr.data')
    ff_corr=ff_corr+load_ff
    fx_corr=fx_corr+load_fx
    dd_corr=dd_corr+load_dd
    os.chdir('..')


ff_corr=ff_corr/index
fx_corr=fx_corr/index
dd_corr=dd_corr/index


#Make diagonal mass matrix.
mass=mass*1.60054*(10**(-27))
mass_matrix=np.zeros([nosc,nosc])
for i in range(nosc):
    mass_matrix[i,i]=mass



#Average done, scale by beta.
k_boltz=1.38064852*(10**(-23))
ff_corr=ff_corr/(k_boltz*temperature)


#Average done, diagonalize averages for ff.
eigval_ff,eigvec_ff=lin.eigh(ff_corr,mass_matrix,eigvals_only=False)


#Eigenvalue problem ff and fx




#Eigenvalue problem for dd (with fx)
dd_inv=lin.inv(dd_corr)
mass_inv=lin.inv(mass_matrix)
b=mass_inv*fx_corr
b=np.matmul(mass_inv,fx_corr)
b_inv=lin.inv(b)
eigval_dd,eigvec_dd=lin.eigh(dd_inv,b_inv,eigvals_only=False)


#Eigenvalue problem for dd without fx, using equipartion to approximate the momentum correlation.
dd_inv=dd_inv*k_boltz*temperature
eigval_dd2,eigvec_dd2=lin.eigh(dd_inv,mass_matrix,eigvals_only=False)



#Convert to wavenumbers
c_light=299792458.0 #meters/second
eigval_ff=np.sqrt(eigval_ff)/(c_light*100*2*np.pi)
eigval_dd=np.sqrt(eigval_dd)/(c_light*100*2*np.pi)
eigval_dd2=np.sqrt(eigval_dd2)/(c_light*100*2*np.pi)



#Print eigenvalue list
np.savetxt('eigenvalues_ff.data',eigval_ff)
np.savetxt('eigenvalues_dd.data',eigval_dd)
np.savetxt('eigenvalues_dd_equi.data',eigval_dd2)


#Calculate rabi
rabi_ff=eigval_ff[-1]-eigval_ff[0]
rabi_dd=eigval_dd[-1]-eigval_dd[0]
rabi_dd2=eigval_dd2[-1]-eigval_dd2[0]

save_ff=open('rabi_ff.data','w')
save_ff.write(str(rabi_ff))
save_ff.close()

save_dd=open('rabi_dd.data','w')
save_dd.write(str(rabi_dd))
save_dd.close()


save_dd2=open('rabi_dd_equi.data','w')
save_dd2.write(str(rabi_dd2))
save_dd2.close()


#Print eigenvectors, normalized
norm_evecs=0
for i in range(nosc):
    norm_evecs=norm_evecs+eigvec_ff[i,0]**2

eigvec_ff=eigvec_ff/np.sqrt(norm_evecs)
np.savetxt('eigenvectors_ff.data',eigvec_ff)


norm_evecs=0
for i in range(nosc):
    norm_evecs=norm_evecs+eigvec_dd[i,0]**2

eigvec_dd=eigvec_dd/np.sqrt(norm_evecs)
np.savetxt('eigenvectors_dd.data',eigvec_dd)
cav_in_pols_dd=np.array([eigvec_dd[0,0]**2,eigvec_dd[0,-1]**2])

#Save the contributions of the cavity to the polaritons.
output=open('cav_in_pols_dd.data','w')
text=str(cav_in_pols_dd[0])+' '+str(cav_in_pols_dd[1])
output.write(text)
output.close()

norm_evecs=0
for i in range(nosc):
    norm_evecs=norm_evecs+eigvec_dd2[i,0]**2

eigvec_dd2=eigvec_dd2/np.sqrt(norm_evecs)
cav_in_pols_dd2=np.array([eigvec_dd2[0,0]**2,eigvec_dd2[0,-1]**2])

#Save the contributions of the cavity to the polaritons.
output=open('cav_in_pols_dd_equi.data','w')
text=str(cav_in_pols_dd2[0])+' '+str(cav_in_pols_dd2[1])
output.write(text)
output.close()



np.savetxt('eigenvectors_dd_equi.data',eigvec_dd2)


#NOW CALCULATE LOCALIZATION, ENTROPY, CAVITY AND SHINE, STORE ALL IN ARRAY AND ENJOY.

#Renormalize eigenvectors for the matter
c_mat=eigvec_ff[1:,:]
c_mat2=c_mat**2

norm_coef=np.sum(c_mat2,axis=0)
norm_coef=np.sqrt(norm_coef)
c_mat=c_mat/norm_coef

#Once it is renormalized, calculate IPS as sum of 4th powers
c_4=c_mat**4
ips_ren=np.sum(c_4,axis=0)


#Calculate shannon entropy: S=sum(p_i*ln(p_i))
#p_i: the probability of each different state(aka, my oscillators)
entropy_matrix=(c_mat**2)*np.log(c_mat**2)
entropy=np.sum(entropy_matrix,axis=0)

#Calculate square of cavity contribution
cavity=eigvec_ff[0,:]**2

#Calculate transition dipolar moment.
dipole=(np.sum(c_mat,axis=0))**2


#Save everything into single matrix
results=np.transpose(eigval_ff)
results=np.vstack([results,ips_ren])
results=np.vstack([results,entropy])
results=np.vstack([results,cavity])
results=np.vstack([results,dipole])

results_print=np.transpose(results)

np.savetxt('freq_ips_entr_cav_dipole.data',results_print)


#NOW for DD with fx (no equipartition).
#ALL COMMENTED OUT DUE TO PROBLEMS WITH BIDIRECTIONAL HARMONIC OSC.
#Renormalize eigenvectors for the matter
c_mat=eigvec_dd[1:,:]
c_mat2=c_mat**2

norm_coef=np.sum(c_mat2,axis=0)
norm_coef=np.sqrt(norm_coef)
c_mat=c_mat/norm_coef

#Once it is renormalized, calculate IPS as sum of 4th powers
c_4=c_mat**4
ips_ren=np.sum(c_4,axis=0)


#Calculate shannon entropy: S=sum(p_i*ln(p_i))
#p_i: the probability of each different state(aka, my oscillators)
entropy_matrix=(c_mat**2)*np.log(c_mat**2)
entropy=np.sum(entropy_matrix,axis=0)

#Calculate square of cavity contribution
cavity=eigvec_dd[0,:]**2

#Calculate transition dipolar moment.
dipole=(np.sum(c_mat,axis=0))**2


#Save everything into single matrix
results=np.transpose(eigval_dd)
results=np.vstack([results,ips_ren])
results=np.vstack([results,entropy])
results=np.vstack([results,cavity])
results=np.vstack([results,dipole])

results_print=np.transpose(results)

np.savetxt('dd_freq_ips_entr_cav_dipole.data',results_print)


#NOW for DD (equipartition).

#Renormalize eigenvectors for the matter
c_mat=eigvec_dd2[1:,:]
c_mat2=c_mat**2

norm_coef=np.sum(c_mat2,axis=0)
norm_coef=np.sqrt(norm_coef)
c_mat=c_mat/norm_coef

#Once it is renormalized, calculate IPS as sum of 4th powers
c_4=c_mat**4
ips_ren=np.sum(c_4,axis=0)


#Calculate shannon entropy: S=sum(p_i*ln(p_i))
#p_i: the probability of each different state(aka, my oscillators)
entropy_matrix=(c_mat**2)*np.log(c_mat**2)
entropy=np.sum(entropy_matrix,axis=0)

#Calculate square of cavity contribution
cavity=eigvec_dd2[0,:]**2

#Calculate transition dipolar moment.
dipole=(np.sum(eigvec_dd2[1:,:],axis=0))**2


#Save everything into single matrix
results=np.transpose(eigval_dd2)
results=np.vstack([results,ips_ren])
results=np.vstack([results,entropy])
results=np.vstack([results,cavity])
results=np.vstack([results,dipole])

results_print=np.transpose(results)

np.savetxt('dd_equi_freq_ips_entr_cav_dipole.data',results_print)


