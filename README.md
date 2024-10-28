# MoMaVSC-pimc
Path Integral Monte Carlo Code for VSC



### Written by
* Jaime de la Fuente
* Riccardo Spezia
* Rodolphe Vuilleumier

### Summary:
This code serves to perform Path Integral Monte Carlo simulations on a simple model of matter under Vibrational Strong Coupling. The model under consideration is a set of one-dimensional oscillators that can be described either by Morse potentials or by harmonic potentials. A linear Hamiltonian coupling scheme of a QED-consistent Pauli-Fierz Hamiltonian coupling can be selected. The code outputs force-force, force-displacement and displacement-displacement correlation functions, from which the analysis scripts are able to compute vibrational spectra and density of states.

### Physical Model

### Files:
The src directory cotains the files that compose the software. The scripts file includes bash scripts that make it easier to put the code into production as well as to analyze the results.
#### Src directory
1) Main.f90: it is the main file of the program. It defines all physical constants, reads the input file and prints progress reports.

#### Scripts directory
