#!/bin/bash

module load anaconda
#Demand parameters
echo 'POSITIONAL PARAMETERS: START_N END_N STEP_N BINS RUNS START_GAUSS END_GAUSS STEP_GAUSS'
temperature_exponent=1
start_gauss=$6
final_gauss=$7
step_gauss=$8
nruns=$5
bin_n=$4
minrange=0
maxrange=4900
red_mass=0.94118 #Reduced mass for an OH oscillator


#Offresonance parameters. This are not received from cline, since they are less common
step=100
max_offres=000
cav_freq=3290
scriptroute="/home/jdelafuentediez/scripts_pimc/scripts_safe_postanalysis" # Route for the script in your setup
	rm rabi*
	rm histogram*
	rm *nonres_*
	rm cavity*
	rm *data
	for offres in `seq $((-$max_offres)) $step $max_offres`
        do
        cd offres_$offres

 
	for t in `seq 1 1 $temperature_exponent`
	do
		#temperature=$((100*$t))
		temperature=300
		cd "$temperature"K
		rm *rabi_vs_nosc_gauss*
		for nosc in `seq $1 $3 $2`
		do
			nosc=$(($nosc**2+1))
			cd nosc_$nosc
				for gauss in `seq $start_gauss $step_gauss $final_gauss`
				do
					cd gauss_$gauss
					rm assembled_energies.data #In order not to overappend if executed 2 times simultaneously.

					for run_id in `seq 1 1 $nruns`
					do
						cd run_$run_id
							rm output*
							#We run here new python algorithm to identify polaritons under disorder insted of the previous approach
							energy_virial=$(grep 'Virial' summary.data | cut -f 11 -d ' ')
							echo $energy_virial>>../assembled_energies.data
						cd ..
					done


					# Bring the scripts to work dir
					cp $scriptroute/corf_analyzer.py .
					cp $scriptroute/histogram.py .
					cp $scriptroute/average_file.py .
					cp $scriptroute/eigenstate_brightness.py .
					
					# Average the corr functions from different equivalent runs. Diagonalize the GEV problem
					python3 corf_analyzer.py --run_id $nruns --nosc $nosc --redmass $red_mass --temp $temperature #All corr functions are generated here

					# Build histograms directly from the output of corf_analyzer.py. Useful if no graphical terminal from cluster to do a quick check
					python3 histogram.py --nbin $bin_n --minfreq $minrange --maxfreq $maxrange #All histograms are generated here

					# Compute the IR brightness of each eigenstate and print on separate file, again useful for quick check
					python3 eigenstate_brightness.py --nbin $bin_n

					# Average the energies of the runs to check how spread are the results of different runs.
					python3 average_file.py -name assembled_energies.data > average_vir_energy.data
					rm assembled_energies.data


					#We store the rabi vs nosc in the temperature folder
#					rabi_ff_newav=$(cat rabi_ff_"$nosc"_"$gauss".data)
					rabi_ff_newav=$(cat rabi_ff.data)
					echo $rabi_ff_newav
#					rabi_dd2=$(cat rabi_dd_equi_"$nosc"_"$gauss".data)
					rabi_dd2=$(cat rabi_dd_equi.data)
					av_light_dark=$(cat av_light_in_dark.data)
					cat av_light_in_dark.data
					echo $rabi_dd2
					echo $nosc $rabi_ff_newav >> ../../"newav_ff_rabi_vs_nosc_gauss_$gauss.data"
					echo $nosc $rabi_dd2 >> ../../"dd_eq_rabi_vs_nosc_gauss_$gauss.data"

					echo $offres 'This is the off resonance index' # In case there is detuning done, tell on which detuning we are
					awk -v w=$offres '{print w" "$1" "$2}' dos_ff.data >> ../../../../histograms_nonres_ff_"$nosc"_"$gauss".data
					awk -v w=$offres '{print w" "$1" "$2}' dos_dd_equi.data >> ../../../../histograms_nonres_dd_eq_"$nosc"_"$gauss".data
					awk -v w=$offres '{print w" "$1" "$2}' cav_in_pols_dd_equi.data >> ../../../../cavity_pols_nonres_dd_eq_"$nosc"_"$gauss".data
					awk -v w=$offres '{print w" "$1" "$2}' cav_in_pols_dd.data >> ../../../../cavity_pols_nonres_dd_"$nosc"_"$gauss".data
					rm cav_in_pols_dd_equi.data #Remove the partial file with the squares.

					#Now for the ir spectrum.
					awk -v w=$offres '{print w" "$1" "$2}' ff_spectr.data >> ../../../../spectra_nonres_ff_"$nosc"_"$gauss".data
					awk -v w=$offres '{print w" "$1" "$2}' dd_equi_spectr.data >> ../../../../spectra_nonres_dd_eq_"$nosc"_"$gauss".data
					
					echo $offres $rabi_ff_newav >> ../../../../rabi_ff_off_"$nosc"_"$gauss".data
					echo $offres $rabi_dd2 >> ../../../../rabi_dd_eq_off_"$nosc"_"$gauss".data
					echo $gauss $av_light_dark >> ../../../../darklight_vs_gauss_nosc_$nosc.data

					cd ..
				done
				cd ..
			done
			cd ..
		done
		cd ..
done

