#!/bin/bash
	# This script is paired with the PIMC-VSC code.
	# It launches multiple calculations trying different oscillator numbers, disorder sizes, detunings and temperatures. It can also 
	#Module that python3 needs
	module load anaconda
#	module load slurm #New cluster needs it
#	module load munge

# CONTROL AREA. MOST THINGS CAN BE DONE WITHOUT DIGGING INTO THE CODE
	#Montecarlo parameters.
	n_runs=25
	ncycles=100000000

	#Offresonance parameters.
	step=200
	max_offres=000
	cav_freq=1671
	osc_freq=1671

	
	#System parameters.
	nbeads=64
	rabi_target=215 #215
	start_gauss=11
	gauss_step=100
	final_gauss=11
	potential='fh2' # Check readme for the keywords
	#Gaussian generation controls.
	generate_random='true'
	upcut=0.25 #With respect to half rabi, dist from P+  to center line to bound the gaussian and cut the tails.
	lowcut=0.25

	#Calculate the deviations from the polaritons to give margin to the freq. distribution.
	upcut=`awk -v cutoff=$upcut -v center=$osc_freq -v rabi=$rabi_target 'BEGIN {print (center+(rabi/2)*(1-cutoff))}'`
	lowcut=`awk -v cutoff=$lowcut -v center=$osc_freq -v rabi=$rabi_target 'BEGIN {print (center-(rabi/2)*(1-cutoff))}'`

	start_dir=`pwd`
	if [ $generate_random == 'true' ]
	then
		rm -rf saved_gaussians
		mkdir saved_gaussians
	fi

	for nosc in `seq $1 $3 $2`
	do

	nosc=$(($nosc**2+1))


	for offres in `seq $((-$max_offres)) $step $max_offres`
	do
		w_cav=`awk -v c_freq=$cav_freq -v offset=$offres 'BEGIN {print c_freq+offset}'`
		echo $w_cav
		mkdir offres_$offres
		cd offres_$offres

		for t in `seq 2 2 2`
		do
			temperature=300
			#temperature=$((3*10**$t))
			mkdir "$temperature"K
			cd "$temperature"K


			#Set individual coupling
			coup=`awk -v sqr_n=$nosc -v omeg_t=$rabi_target 'BEGIN {print omeg_t/sqrt(sqr_n-1)}'` #Static gamma
			# coup=`awk -v sqr_n=$nosc -v w_c=$w_cav -v w_o=$osc_freq -v omeg_t=$rabi_target 'BEGIN {print (w_c/w_o)^(3/2)*(omeg_t/sqr_n)}'` #Dynamic gamma for the off resonance

		
			mkdir nosc_$nosc
			cd nosc_$nosc


			for gauss in `seq $start_gauss $gauss_step $final_gauss`
			do
				mkdir gauss_$gauss
				cd gauss_$gauss
				if [ $generate_random == 'true' ]
				then
					#Get the generator for the corrected gaussian.
						cp ~/custom_commands/normal_corrected.py .

					#Use the generator
						python3 normal_corrected.py --up_cut=$upcut --low_cut=$lowcut --avg=$osc_freq --desvest=$gauss --nosc=$nosc

					#Set results in nice name for fortran90
						cp inp_freqs_"$nosc"_"$gauss".data starting_freqs.data
						mv inp_freqs_"$nosc"_"$gauss".data inp_freqs_nosc"$nosc"_gauss"$gauss"_temp"$temperature"_off"$offres".data

					#Move to our safe folder of distributions
						mv inp_freqs_nosc"$nosc"_gauss"$gauss"_temp"$temperature"_off"$offres".data $start_dir/saved_gaussians/.

					#Delete generator for cleanness
						rm normal_corrected.py
				else
					#If we don't generate, we read from saved_gaussians. There must be one for our precise case
						cp "$start_dir"/saved_gaussians/inp_freqs_nosc"$nosc"_gauss"$gauss"_temp"$temperature"_off"$offres".data starting_freqs.data

				fi

				for run_id in `seq 1 1 $n_runs`
				do
					mkdir run_$run_id
					cd run_$run_id
					cp ../starting_freqs.data .
					#Generate slurm launcher file

					cat <<-EOF >launcher
						#!/bin/bash
						#SBATCH -J nosc_$nosc
						#SBATCH --ntasks=1
						#SBATCH --time=10000:00:00
						#SBATCH --output coup_point_$nosc.out
						#SBATCH --error coup_point_$nosc.err
						#SBATCH --partition=main 

						#########################################################

						#########module purge
						#########module load gcc


						mycode=/home/jdelafuentediez/momavsc-harm/bin/code.exe #Route for where your compiled code is stored
						time \$mycode < input.in > output.$nosc.data


					EOF

					#Generate input file
					cat <<-EOF >input.in
						&seetings
						model=$potential
						rabi_split=$coup
						nosc=$nosc
						nbeads=$nbeads
						temp=$temperature
						w_oh=$osc_freq
						mass_1=16.d0
						mass_2=1.d0
						cycles=$ncycles
						step=0.2d0
						start_sampling=0.3d0
						estimator_interval=10
						freq_coll=0.2d0
						verbose=.False.
						width_freq=$gauss.d0
						w_cav=$w_cav
						de=500.d0
						outlier=.True.
						outlier_freq=1794.d0
						load_rand=.True.
						&end

					EOF
					sbatch launcher
					cd ..

				done
				cd ..
			done
			cd ..
			cd ..
		done
		cd ..
	done
done
