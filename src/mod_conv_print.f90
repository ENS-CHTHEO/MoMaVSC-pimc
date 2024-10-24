module mod_conv
        implicit none
        !Phys. constants
	real*8,parameter::pi=DACOS(-1.d0)
	real*8,parameter::hbar=1.054571628D-34 !h/2pi
	real*8,parameter::avog=6.02214179D23 !Avogadro's number
	real*8,parameter::h=6.62606896D-34 !Planck constant
	real*8,parameter::c_vac=2.99792458D8 !Light speed vacuum

        !Conversion factors (hardcoded)
	real*8,parameter::atmass_kg=1.60054D-27
	real*8,parameter::angs_m=1.0D-10
	real*8,parameter::wnum_hz=100.d0*c_vac !Wavenumber(cm-1) to hertz.
	real*8,parameter::wnum_e=h*c_vac*100.d0 !Wavenumber(cm-1) to energy(Joule).
	real*8,parameter::jul_hartree=2.293710449D17 !Joule to hartree

        public::conversion,heading
        contains


        !Receives all physical constants and converts them all to SI units for internal work or vice versa
        subroutine conversion(in_or_out,mass,freq,nbead,nosc,pos,step,coup_const,e1,e2,req,de)
                implicit none

                !Entry variables
                integer,intent(in)::in_or_out,nosc,nbead
		real*8,intent(inout),dimension(nosc)::freq
		real*8,intent(inout)::req
		real*8,intent(inout)::mass,step,coup_const,de
		real*8,intent(inout),dimension(nosc,nbead)::pos
		real*8,intent(inout)::e1,e2 !Energies to be converted




                if (in_or_out==1) then !entry units are turned into SI units for iternal use
                        mass=mass*atmass_kg !Set mass to kg
                        freq(:)=freq(:)*wnum_hz !Set frequency to hertz
                        coup_const=coup_const*wnum_hz !Coupling constant from cm-1 to hertz
                        pos(:,:)=pos(:,:)*angs_m !Coordinates in meters
                        step=step*angs_m !Stepsize in meters
                        req=req*angs_m
                        de=de*1000/avog

                else if (in_or_out==2) then !Units reconverted to entry units for printout
                        pos(:,:)=pos(:,:)/angs_m !Coordinates in angstrom again
                        e1=e1*jul_hartree !enegies in joules go out in hartree
                        e2=e2*jul_hartree
                else
                        write(6,*) 'Mistake on the conversion mode. Please use 1 or 2'
                endif
        end subroutine conversion


        subroutine heading(nosc,nbeads,cycles,estimator_interval,temp,w_oh,red_mass,start_sampling,step,freq_coll,coup_const)
                implicit none
                integer,intent(in)::nosc,nbeads,estimator_interval
                integer(kind=8),intent(in)::cycles
		real*8,intent(in)::temp,w_oh,red_mass,start_sampling,step,freq_coll,coup_const
		real::useful_moves

                open(33,file='summary.data')
		write(33,*) '           ***********************************************************************'
		write(33,*) '             -------PATH INTEGRAL MONTECARLO COUPLED HARMONIC OSCILLATORS-------'
		write(33,*) '                        (Looking for any other more artistic name)    '
		write(33,*) '                                         v0.1                          '
		write(33,*) '             ***By: Jaime de la Fuente, Rodolphe Vuilleumier & Riccardo Spezia***'
		write(33,*) '                    ***École Normale Superiéure, Sorbonne Université***'
		write(33,*) '           ***********************************************************************'

		write(33,*) ''
		write(33,*) ''
		write(33,*) ''
		write(33,*) 'Data read successfully from input file!'
		write(33,*) 'Seetings for the current simulation: '
		write(33,*) ''
		write(33,*) ''
		write(33,*) 'Number of harmonic oscillators:',nosc
		write(33,*) ''
		write(33,*) 'Oscillators frequency (cm^-1):',w_oh
		write(33,*) ''
		write(33,*) 'Coupling constant (cm-1)', coup_const
		write(33,*) ''
		write(33,*) 'Reduced mass of the oscillator (atomic units):', red_mass
		write(33,*) ''
		write(33,*) 'Number of beads considered (i time slices):',nbeads
		write(33,*) ''
		write(33,*) 'Fixed temperature for the system (Kelvin)',temp
		write(33,*) ''
		write(33,*) 'Total number of Monte Carlo trial moves:',cycles
		write(33,*) ''
		useful_moves=(1.d0-start_sampling)*100
		write(33,*) 'Proportion of MonteCarlo moves acceptable for production (%):',useful_moves
		write(33,*) ''
		write(33,*) 'Sampling interval of acceptable moves (each how many moves energy is updated):',estimator_interval
		write(33,*) ''
		write(33,*) 'Maximum montecarlo movement size:',step
		write(33,*) ''
		write(33,*) 'Frequency at which collective bead movements are tried instead of single ones (%):',freq_coll*100
		write(33,*) ''
		write(33,*) 'Launching MonteCarlo...'
		write(33,*) '----------------------------------------------------------------------------------------------------'
                close(33)
	end subroutine heading



	subroutine results(av_energy,av_virial,av_analytical,accept_count,sampled_confs,nosc,sigma_e,sigma_vir)
		implicit none
		real*8,intent(in)::av_energy,av_virial,av_analytical,sigma_e,sigma_vir
		integer,intent(in)::nosc
		integer(kind=8),intent(in)::accept_count,sampled_confs

                open(33,file='summary.data')
		write(33,*) '----------------------------------------------------------------------------------------------------'
		write(33,*) ''
		write(33,*) ''
		write(33,*) 'Movement generation ended successfully!!!'
		write(33,*) ''
		write(33,*) '**Results summary:'
		write(33,*) '->Average analytic harmonic oscillator energy (joule) for:',av_analytical
		write(33,*) ''
		write(33,*) '->Average energy from MC simulation (average of observable):',av_energy
		write(33,*) '->Std deviation for this 1st estimator',sigma_e
		write(33,*) ''
		write(33,*) '->Average energy from MC simulation (Virial theorem):',av_virial
		write(33,*) '->Std deviation for this 2nd estimator:',sigma_vir
		write(33,*) ''
		write(33,*) '->Configurations sampled for obtaining estimators:',sampled_confs
		write(33,*) ''

		write(33,*) '->Total accepted MC configurations:',accept_count
		write(33,*) ''

		close(33)
	end subroutine results


end module mod_conv
