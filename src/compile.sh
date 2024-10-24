#!/bin/bash
module load scalapack
module load openblas
module load gcc
export 
####gfortran -O4 -o ../bin/code.exe pot_mod.f90 mc_mod.f90 mod_conv_print.f90 main.f90 -L$OPENBLAS_PATH/lib -lopenblas -fallow-argument-mismatch
gfortran -O4 -o ../bin/code.exe pot_mod.f90 mc_mod.f90 mod_conv_print.f90 main.f90 -L$OPENBLAS_PATH/lib -lopenblas -Wno-argument-mismatch

