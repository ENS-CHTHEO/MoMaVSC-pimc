program path_montecarlo
        use mc_mod
        use mod_conv
        implicit none

	real*8,parameter::pi2=dacos(-1.d0)
        !Simulation parameters
        integer::nbeads=16
        integer(kind=8)::cycles=100000 !Number of "beads" (i time slices) and MC cycles
	real*8::coup_const=0.d0,rabi_split=0.d0
        integer::nosc=1 !Number of oscillators
        integer::nblocks=10,blocksize
        character(len=3)::model='har'
	real*8::mass_1=16.d0 !Mass in atomic units of O16 (no isotopic averaging)
	real*8::mass_2=1.d0 !Proton mass (no isotopic averaging)
	real*8::mass_red
	real*8::w_oh=3404.d0 !Resonant to which cavity is resonant. cm-1
	real*8::w_cav=3404.d0 !Resonant to which cavity is resonant. cm-1
	real*8,dimension(:),allocatable::w_range !Matrix for molecule frequencies (around w_oh by width_freq)
        real*8::width_freq=0.d0 !Width of frequency distribution of molecules.cm-1
	real*8::temp=300.d0 !Temperature in Kelvin
        real*8::de=500!Dissoc energy of oh in water in kj/mol
        real*8::req=0.d0 !In angstrom
	real*8::step=0.1d0
	real*8::pol_action
	real*8::energy_vir
	real*8::energy_clas
	real*8::conf_energy
        integer::estimator_interval=10 !Each how many configurations the energies are updated
	real*8::start_sampling=0.3d0 !After which % of cycles sampling starts. 30% in this case
	real*8::production_cycles
        integer(kind=8)::sampled_confs
	real*8::av_energy,av_virial,av_analytical,die_const,rabisplit_fd
	real*8::freq_coll=0.5d0
	real*8::sigma_e,sigma_vir,av_e_block,av_vir_block
        integer::nblock1=4,nblock2=16,nblock3=32
	logical::verbose=.False. !Controls if detailed data files are written or not.
        logical::do_blocking=.False.  !Turns on of off blocking analysis for energies.
        logical::outlier=.False. !Allows to introduce an outlier frequency associated to the 2nd oscillatos (first matter one)
        real*8::outlier_freq
        logical::configs=.False. !Default for configuration printing

        !Simulation variables
	real*8::b
	real*8,allocatable,dimension(:,:)::pos,av_xx
	real*8::r2
        real*8::upper_pos,lower_pos,binsize
        real*8::binsize_spectr
        integer(kind=8)::accept_count
        integer::nbins=200
        integer::nbins_spectr
        integer::nmoves=2
        integer,dimension(:),allocatable::accepted_moves
        real*8,allocatable,dimension(:)::histogram
	real*8,allocatable,dimension(:)::grad,av_x
	real*8,allocatable,dimension(:,:)::blocking
	real*8,allocatable,dimension(:,:)::force_corr,mass_diag
	real*8,allocatable,dimension(:)::eigen
        real*8::avfieldsq,var_avfieldsq,acum_pos
        real*8,dimension(:,:),allocatable::f_pos_corr
        real*8,dimension(:,:),allocatable::disp_disp_corr


        !Variables for force force diagonalization
        real*8,dimension(:),allocatable::eigen_ff
        real*8,dimension(:,:),allocatable::f_corr_saved
        integer*8::lwork_ff,inform_ff
        real*8::optimum_ff
        real*8,dimension(:),allocatable::work_ff


        !Auxiliaries
        integer(kind=8)::i,j,k,aux_index,inform,lwork
	real*8::optimum,coordinate,norm_evecs
        real*8,dimension(:),allocatable::var_al1,var_al2
        integer,parameter::into=1,outof=2
	real*8,dimension(:),allocatable::history,work !Stores action values for each iteration
        real*8::aux
        logical::load_rand=.True.
        integer::print_it
        !Namelist
        namelist /seetings/ model,rabi_split,nosc,nbeads,temp,mass_1,mass_2,w_oh,cycles,step,start_sampling,estimator_interval,&
                freq_coll,nblock1,nblock2,nblock3,verbose,width_freq,w_cav,de,outlier,outlier_freq,load_rand,configs
        read(5,seetings)

        mass_red=mass_1*mass_2/(mass_1+mass_2)
        production_cycles=cycles*(1.d0-start_sampling)
        sampled_confs=int(production_cycles/estimator_interval)
        blocksize=int(sampled_confs/nblocks)
        call heading(nosc,nbeads,cycles,estimator_interval,temp,w_oh,mass_red,start_sampling,step,freq_coll,coup_const)


!--------------------PROGRAM START--------------------------
        call random_seed()
        b=1.d0/(kb*temp)
        allocate(pos(nosc,nbeads),grad(nosc),force_corr(nosc,nosc),eigen(nosc),mass_diag(nosc,nosc),f_pos_corr(nosc,nosc), &
                accepted_moves(nmoves),disp_disp_corr(nosc,nosc),av_xx(nosc,nosc),av_x(nosc))
        if (verbose.eqv..True.) then
                allocate(history(cycles))
        endif

        !Initialize to 0
        pos(:,:)=0.d0
        energy_vir=0.d0
        energy_clas=0.d0
        accept_count=0


        !Generation of Gaussian distribution of frequencies around w_oh.
        allocate(var_al1(nosc),var_al2(nosc),w_range(nosc))

        !Load the starting random numbers common to all replicas so that all have the same start distribution.

        if (load_rand.eqv..True.) then
                open(11,file='starting_freqs.data')
                read(11,*) w_range(:)
                close(11)
                deallocate(var_al1,var_al2) !No longer needed if we get the distribution previously generated.

        else 
                call random_number(var_al1(:))
                call random_number(var_al2(:))
                w_range(:)=dsqrt(-2.d0*dlog(var_al1(:)))*dcos(2.d0*pi2*var_al2(:)) !Box muller transform.
                w_range(:)=0.d0+w_range(:)*width_freq!Now, set the variance on the normal distribution. Average is 0.
                w_range(:)=w_oh+w_range(:) !Apply normal displacements around the average.
        endif

        w_range(1)=w_cav !Set frequency for cavity to exactly itsef.

        write(6,*) w_range(:)
        open(11,file='used_freqs.data') !Overwrite with the actual freq distribution used
                write(11,*) w_range(:)
        close(11)




        !Introduce outlier (might break the ordering of the frequencies from the previous loop.
                if (outlier) then
                        w_range(2)=outlier_freq
                endif



        !Convert now the frequencies to SI
        call conversion(into,mass_red,w_range,nbeads,nosc,pos,step,rabi_split,energy_vir,energy_clas,req,de)


        !Calculate coupling const.
        coup_const=mass_red*w_range(1)*rabi_split*4.d0*pi2**2 !Coupling constant depends on the cavity not the oscillators!!.


        !Action for the starting point
        call calc_action(model,nbeads,nosc,pos,b,mass_red,w_range,coup_const,pol_action,de,req)

        !Sampling
        av_energy=0.d0 !Through 1st energy estimator
        av_virial=0.d0 !Through 2nd energy estimator (virial formula)
        av_analytical=h*w_oh*(0.5d0+1.d0/(dexp(h*w_oh*b)-1.d0))*nosc !Average energy for a set of uncoupled harmonic oscillators
        force_corr(:,:)=0.d0
        f_pos_corr(:,:)=0.d0 
        av_xx=0.d0
        av_x=0.d0
        aux_index=1

        !Histogram for the electric field
        allocate(histogram(nbins))
        histogram(:)=0

        upper_pos=3D-11
        lower_pos=-3D-11
        binsize=(upper_pos-lower_pos)/nbins
        avfieldsq=0.d0
        die_const=0.d0

        if (do_blocking.eqv..True.) then
                allocate(blocking(sampled_confs,3))
        endif

        accepted_moves=0
        open(17,file='pos_nosc_bead.data')
        do i=1,cycles

                call mc_update(model,nbeads,nosc,pos,step,pol_action,accept_count,b,w_range,coup_const,mass_red,freq_coll,de,&
                        req,accepted_moves,nmoves)

                if ((mod(i,estimator_interval)==0).and.(i>(cycles-production_cycles))) then !Calculate estimator each n cycles
                call energies(model,nbeads,nosc,pos,grad,b,mass_red,w_range,coup_const,conf_energy,energy_vir,force_corr, & 
                        f_pos_corr,de,req,av_xx,av_x)

                        av_energy=av_energy+conf_energy
                        av_virial=av_virial+energy_vir

                        if (do_blocking.eqv..True.) then
                        !Saving into data array for blocking analysis
                                blocking(aux_index,1)=pol_action/b
                                blocking(aux_index,2)=conf_energy
                                blocking(aux_index,3)=energy_vir
                                aux_index=aux_index+1
                        endif

                        !Calculating distribution of the electric field
                        do j=1,nbeads
                                k=floor((pos(1,j)-lower_pos)/binsize)
                                if ((k<=nbins).and.(k>0)) then
                                        histogram(k)=histogram(k)+1
                                endif
                                avfieldsq=avfieldsq+pos(1,j)**2
                        enddo

                        acum_pos=0.d0
                        do j=1,nosc
                                do k=1,nbeads
                                        acum_pos=pos(j,k)+acum_pos
                                enddo
                        enddo

                        die_const=(acum_pos**2)+die_const

                	if ((configs.eqv..True.) .and. (mod(i,1000)==0) ) then
                        	write(17,*) 'MC Step: ',i
                        	do print_it=1,nosc
                                	write(17,*) pos(print_it,:)
                        	enddo
                	endif

                endif


                if (verbose.eqv..True.) then
                        write(6,*) 'iteration',i,'exponent',pol_action,'action',pol_action/b
                        history(i)=pol_action !register action value
                endif

        enddo
        close(17)

        !Write the acceptances
        open(11,file='acceptations.data')
        do i=1,nmoves
                write(11,*) i,accepted_moves(i)
        enddo
        write(11,*) '#Total accept count',accept_count
        write(11,*) '#Proportion of accepted move 1 vs attempted move 1',(accepted_moves(1)/(cycles*(1-freq_coll))*100)
        write(11,*) '#Proportion of accepted move 2 vs attempted move 2',(accepted_moves(2)/(cycles*freq_coll)*100)
        close(11)

        ! Calculate average E.
        av_energy=av_energy/sampled_confs
        av_virial=av_virial/sampled_confs

        !Average correlation functions along the ensemble.
        force_corr=force_corr/sampled_confs !Actually turning force corr in potential corr
        f_pos_corr=f_pos_corr/sampled_confs
        av_xx=av_xx/sampled_confs
        av_x=av_x/sampled_confs

        !Calculate displacement displacement corr
        do i=1,nosc
                do j=1,nosc
                        disp_disp_corr(i,j)=av_xx(i,j)-av_x(i)*av_x(j)
                enddo
        enddo

        !Normalize histogram
        histogram=histogram/(sampled_confs*nbeads)
        avfieldsq=avfieldsq/(sampled_confs*nbeads) !calculating average field^2


        !Calculate dielectric const.
        die_const=b*die_const/((nbeads**2)*sampled_confs)

        !Get variance:
        var_avfieldsq=0.d0
        do i=1,nbins
                var_avfieldsq=var_avfieldsq+((i-1)*binsize+lower_pos)**2*histogram(i)
        enddo

        var_avfieldsq=var_avfieldsq/nbins
        var_avfieldsq=var_avfieldsq-avfieldsq

        open(22,file='force_corr.data')
                do i=1,nosc
                        write(22,*) force_corr(i,:)
                enddo
        close(22)



        !Print the force-displacement correlation

        open(31,file='force_displacement_corr.data')
        do i=1,nosc
                write(31,*) f_pos_corr(i,:)
        enddo
        close(31)


        !Print the displacement displacement correlation

        open(31,file='disp_disp_corr.data')
        do i=1,nosc
                write(31,*) disp_disp_corr(i,:)
        enddo
        close(31)




!        !Lapack diagonalization (dsygv) of force-force eigenvalue problem

        mass_diag=0.d0 !Diagonal matrix with the masses
        do i=1,nosc
                mass_diag(i,i)=mass_red
        enddo
        allocate(eigen_ff(nosc),f_corr_saved(nosc,nosc))     
        f_corr_saved=b*force_corr
   
        lwork_ff=-1
       
        call dsygv(1,'N','U',nosc,f_corr_saved,nosc,mass_diag,nosc,eigen_ff,optimum_ff,lwork_ff,inform_ff)

        lwork_ff=floor(optimum_ff)
        allocate(work_ff(lwork_ff))
        call dsygv(1,'V','U',nosc,f_corr_saved,nosc,mass_diag,nosc,eigen_ff,work_ff,lwork_ff,inform) !The eigenvalues come out in hertz.


        eigen_ff(:)=dsqrt(eigen_ff(:))/(c_vac*100.d0*2.d0*pi2) !Turning the frequencies from hertz to cm^-1
        open(18,file='eigenvalues_ff.data')
        do i=1,nosc
                write(18,*) eigen_ff(i)
        enddo
        close(18)


        norm_evecs=0.d0
        open(19,file='eigenvectors_ff.data')
        do i=1,nosc
                norm_evecs=norm_evecs+f_corr_saved(i,1)**2
        enddo
        f_corr_saved(:,:)=f_corr_saved(:,:)/dsqrt(norm_evecs)

        do i=1,nosc !Print the matrix with eigenstates.
                write(19,*) f_corr_saved(i,:)
        enddo
        close(19)




        !Lapack diagonalization with dsygv of force-displacement correlation eigenvalue problem.
        mass_diag=0.d0 !Diagonal matrix with the masses
        do i=1,nosc
                mass_diag(i,i)=mass_red
        enddo

        lwork=-1
        f_pos_corr=matmul(mass_diag,f_pos_corr)
 
        call dsygv(1,'N','U',nosc,force_corr,nosc,f_pos_corr,nosc,eigen,optimum,lwork,inform)
        lwork=floor(optimum)
        allocate(work(lwork))
        call dsygv(1,'V','U',nosc,force_corr,nosc,f_pos_corr,nosc,eigen,work,lwork,inform) !The eigenvalues come out in hertz.

        !Normalization of the force-displacement eigenvectors
        norm_evecs=0.d0
        open(19,file='eigenvectors_fx.data')
        do i=1,nosc
                norm_evecs=norm_evecs+force_corr(i,1)**2
        enddo
        force_corr(:,:)=force_corr(:,:)/dsqrt(norm_evecs)

        do i=1,nosc !Print the matrix with eigenstates.
                write(19,*) force_corr(i,:)
        enddo
        close(19)



        !Conversion of the eigenvalues
        eigen(:)=dsqrt(eigen(:))/(c_vac*100.d0*2.d0*pi2) !Turning the frequencies from hertz to cm^-1

        open(33,file='eigenvalues_fx.data')
                do i=1,nosc
                        write(33,*) eigen(i) !Write a colon of eigenvalues
                enddo 
        close(33)





        !WRITING RESULTS IN FILES
        !Write the electric field as histogram.
        open(33,file='e_field.data')
                write(33,*) "#",avfieldsq
                write(33,*) "#",var_avfieldsq
                do i=1,nbins
                        coordinate=(i-1)*binsize+lower_pos
                        write(33,*) coordinate, histogram(i)
                enddo

        close(33)



        !Writing data
        if (verbose.eqv..True.) then
                open(11,file='mc_history_action.data')
                do i=1,cycles
                        write(11,*) i,history(i)
                enddo
                close(11)
        endif



        !Write dielectric constant
        open(26,file='dielectric_constant.data')
                write(26,*) die_const
        close(26)




        !STATISTICS ON THE RESULTS
        !Blocking analysis of the results
        if (verbose.eqv..True.) then
                open(12,file='blocking.data') !Save sampled data to file
                do i=1,sampled_confs
                        write(12,*) i, blocking(i,:)
                enddo
         write(12,*) 'final', sum(blocking(:,1))/sampled_confs,sum(blocking(:,2))/sampled_confs,sum(blocking(:,3))/sampled_confs
                close(12)
        endif


        !Distribute data in blocks
        sigma_e=0.d0
        sigma_vir=0.d0
        if (do_blocking.eqv..True.) then
                blocksize=int(sampled_confs/nblock1)
                open(14,file='block_averages1.data')
                do i=0,nblock1-1
                        av_e_block=0.d0
                        av_vir_block=0.d0

                        do j=1,blocksize
                                av_e_block=av_e_block+blocking(i*blocksize+j,2)
                                av_vir_block=av_vir_block+blocking(i*blocksize+j,3)
                        enddo
                        av_e_block=av_e_block/blocksize
                        av_vir_block=av_vir_block/blocksize
                        write(14,*) i,av_e_block,av_vir_block


                        sigma_e=sigma_e+(av_e_block-av_energy)**2
                        sigma_vir=sigma_vir+(av_vir_block-av_virial)**2
                enddo
                close(14)

        sigma_e=dsqrt(1.d0/(nblock1-1)*sigma_e)
        sigma_vir=dsqrt(1.d0/(nblock1-1)*sigma_vir)

        call results(av_energy,av_virial,av_analytical,accept_count,sampled_confs,nosc,sigma_e,sigma_vir)

        open(13,file='block_results1.data')
                write(13,*) blocksize,sigma_e,sigma_vir
        close(13)


        sigma_e=0.d0
        sigma_vir=0.d0
        blocksize=int(sampled_confs/nblock2)
        open(14,file='block_averages2.data')
        do i=0,nblock2-1
                av_e_block=0.d0
                av_vir_block=0.d0

                do j=1,blocksize
                        av_e_block=av_e_block+blocking(i*blocksize+j,2)
                        av_vir_block=av_vir_block+blocking(i*blocksize+j,3)
                enddo
                        av_e_block=av_e_block/blocksize
                        av_vir_block=av_vir_block/blocksize
                        write(14,*) i,av_e_block,av_vir_block


                        sigma_e=sigma_e+(av_e_block-av_energy)**2
                        sigma_vir=sigma_vir+(av_vir_block-av_virial)**2
        enddo
        close(14)

        sigma_e=dsqrt(1.d0/(nblock2-1)*sigma_e)
        sigma_vir=dsqrt(1.d0/(nblock2-1)*sigma_vir)

        open(13,file='block_results2.data')
                write(13,*) blocksize,sigma_e,sigma_vir
        close(13)


        sigma_e=0.d0
        sigma_vir=0.d0
        blocksize=int(sampled_confs/nblock3)
        open(14,file='block_averages3.data')
        do i=0,nblock3-1
                av_e_block=0.d0
                av_vir_block=0.d0

                do j=1,blocksize
                        av_e_block=av_e_block+blocking(i*blocksize+j,2)
                        av_vir_block=av_vir_block+blocking(i*blocksize+j,3)
                enddo
                        av_e_block=av_e_block/blocksize
                        av_vir_block=av_vir_block/blocksize
                        write(14,*) i,av_e_block,av_vir_block


                        sigma_e=sigma_e+(av_e_block-av_energy)**2
                        sigma_vir=sigma_vir+(av_vir_block-av_virial)**2
        enddo
        close(14)

        sigma_e=dsqrt(1.d0/(nblock3-1)*sigma_e)
        sigma_vir=dsqrt(1.d0/(nblock3-1)*sigma_vir)

        open(13,file='block_results3.data')
                write(13,*) blocksize,sigma_e,sigma_vir
        close(13)
        
endif


call results(av_energy,av_virial,av_analytical,accept_count,sampled_confs,nosc,sigma_e,sigma_vir)

write(6,*) "I'm done!!"

end program path_montecarlo
