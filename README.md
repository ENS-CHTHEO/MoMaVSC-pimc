# MoMaVSC-pimc
Path Integral Monte Carlo Code for VSC



### Written by
* Jaime de la Fuente
* Riccardo Spezia
* Rodolphe Vuilleumier

### Summary:
This code serves to perform Path Integral Monte Carlo simulations on a simple model of matter under Vibrational Strong Coupling. The model under consideration is a set of one-dimensional oscillators that can be described either by Morse potentials or by harmonic potentials. A linear Hamiltonian coupling scheme of a QED-consistent Pauli-Fierz Hamiltonian coupling can be selected. The code outputs force-force, force-displacement and displacement-displacement correlation functions, from which the analysis scripts are able to compute vibrational spectra and density of states.

This foundations of this code have been presented at: J. Chem. Phys. 161, 184114 (2024), for which it has produced most of the results.

This code has also been used to provide results in:


### Files:
The src directory cotains the files that compose the software. The scripts file includes bash scripts that make it easier to put the code into production as well as to analyze the results.
#### Src directory
1) Main.f90: it is the main file of the program. It defines all physical constants, reads the input file and prints progress reports.
2) mc_mod.f90: it is the module containing the monte carlo related subroutines used to compute the energy of a configuration, the action of a configuration, and to generate and accept/reject new configurations.
3) pot_mod.f90: functions to compute potential and gradient of matter and cavity coupled by either linear of Pauli\--Fierz Hamiltonian (for a description of the Hamiltonians see the paper referenced at the start). Different descriptions of matter are given by different potential names. In particular:
   - har: harmonic potential for matter. Linear coupling.
   - mor: morse potential for matter. Linear coupling.
   - fhr: harmonic potential for matter. PF coupling.
   - fmr: morse potential for matter. PF coupling.
   - fh2: harmonic potential for matter. PF coupling. Random orientation of molecules.
   - fm2: morse potential for matter. PF coupling. Random orientation of molecules.
#### Scripts directory
There are two scripts and a folder containing small snippets of Python. The code is much easier to use if it is paired with the scripts. 
1) runs_to_avg_v5.sh: It launches multiple replicas of the calculations using a slurm queue system. It has to be edited to be used. The top part is the only one that should be edited. One can specify the number of runs with different random seeds to be done, the number of monte carlo steps of each run, the cavity frequency and the targeted Rabi Splitting. The cavity can be sequentially detuned by defining a detuning range and step. Also, one can generate a frequency distribution inside the code for matter frequencies, which will be a trimmed gaussian and thus its average and standard deviation can also be defined. Alternatively, one can load a frequency distribution by setting the parameter "generate_random" to true or false respectively.

The script creates a set of directories, at the end of which the datafiles will be found. The structure of the directories is always, starting at the launching directory: "offres_$offres/$temperature K/nosc_$numberofmolecules/gauss_$stdevofmatter/nrun_$indexofrepication/"

"offres" indicates the detuning of the cavity from its central value, and is given by $offres, the temperature is also indicated. Then $numberofmolecules is the number of 1 dimensional modes representing matter, $stdevofmatter is the standard deviation of the matter distribution, which must be always indicated, even if the distribution is loaded from the outside, and $indexofreplication is the index of the run replica, if multiple different MC runs that differ in their random seed have been launched.

2) runs_to_avg_bipartition.sh is a version of the previous script that can automatically generate two gaussian distributions of matter at different frequencies, with the aim of studying vibrational strong co-coupling and the overlap between different spectral bands under VSC.

3) result_processor_nonres.sh: analyzes the results. When it is launched, it indicates the arguments that it needs to work. The analysis consists in the averaging of correlation functions from different runs, and the solving of the generalized eigenvalue problem associated with them. It must be launched from the same directory where the calculation is launched. At the gauss_$stdev directory, the eigenvalues and eigenvectors calculated in different ways can be found. The most commonly used are eigenvalues_dd_equi.data and eigenvectors_dd_equi.data 
