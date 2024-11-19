module mc_mod
        use pot_mod
        implicit none
	real*8,parameter::pi=dacos(-1.d0)
	real*8,parameter::hbar=1.054571628D-34 !h/2pi
	real*8,parameter::kb=1.3806504D-23 !Boltzmann constant
        public::energies,calc_action,mc_update
        contains


        !Both whith the direct average value and the virial approach
        subroutine energies(model,nbeads,nosc,pos,grad,b,mass,freq_range,coup_const,conf_energy,energy_vir,force_corr, &
                        f_pos_corr,de,req,av_xx,av_x)
                implicit none
                character(len=3)::model
                integer,intent(in)::nbeads,nosc
		real*8,intent(in)::b,de,req
		real*8,intent(in),dimension(nosc-1)::coup_const
		real*8,dimension(nosc,nbeads),intent(in)::pos
		real*8,intent(in)::mass
		real*8,allocatable,dimension(:),intent(in)::freq_range
		real*8,intent(inout)::energy_vir,conf_energy
		real*8::r2
                real*8,dimension(nosc)::kforce
		real*8,dimension(nosc)::deltar
		real*8,dimension(nosc)::gradients,displacements
		real*8,dimension(nosc),intent(inout)::grad
		real*8,dimension(nosc,nosc),intent(inout)::force_corr,f_pos_corr
                real*8,dimension(nosc,nosc)::f_pos_corr1,f_pos_corr2,av_xx_part
                real*8,dimension(nosc,nosc),intent(inout)::av_xx
                real*8,dimension(nosc),intent(inout)::av_x
                integer::i,j,k

        !Calculate force constant (classical equation):
                kforce(:)=(4.d0*pi**2)*mass*freq_range(:)**2!kf=4*pi^2*m*w^2

        !1st estimator current configuration contribution
                conf_energy=nbeads/(2.d0*b)*nosc
                do i=1,nbeads
                        if (i==nbeads) then
                                deltar(:)=pos(:,i)-pos(:,1)
                        else
                                deltar(:)=pos(:,i)-pos(:,i+1)
                        endif

                        r2=sum(deltar(:)**2)

                        conf_energy=conf_energy-mass*nbeads/(2.d0*(hbar**2)*(b**2))*r2
                enddo

                do i=1,nbeads
                        conf_energy=conf_energy+calc_pot(model,nosc,pos(:,i),kforce,coup_const,de,req)/nbeads
                enddo


        !Calculating energy through virial
                energy_vir=0.d0
                do i=1,nbeads
                                call calc_grad(model,nosc,pos(:,i),kforce,coup_const,de,req,grad)
                                do j=1,nosc
                                        grad(j)=grad(j)*pos(j,i)
                                enddo

                                energy_vir=energy_vir+0.5d0*sum(grad(:))+calc_pot(model,nosc,pos(:,i),kforce,coup_const,de,req)
                enddo
                energy_vir=energy_vir/nbeads


        !Calculating force force correlations
        gradients(:)=0.d0
        do i=1,nbeads
                call calc_grad(model,nosc,pos(:,i),kforce,coup_const,de,req,grad)
                gradients(:)=gradients(:)+grad(:) !Accumulates all the gradients at every oscillator and bead
        enddo
        gradients=gradients/nbeads

        !Consider particle pairs
        do i=1,nosc
                do j=i,nosc
                        force_corr(i,j)=force_corr(i,j)+gradients(i)*gradients(j) !Update force correlation matrix element
                        force_corr(j,i)=force_corr(i,j) !Symmetry of matrix
                enddo
        enddo

        !Calculating force_position correlation
        gradients(:)=0.d0
        displacements(:)=0.d0

        !Calculate average displacements and average gradients on the oscillators.
        do i=1,nbeads
                call calc_grad(model,nosc,pos(:,i),kforce,coup_const,de,req,grad)
                gradients(:)=gradients(:)+grad(:)
                displacements(:)=displacements(:)+pos(:,i)
        enddo

        f_pos_corr1(:,:)=0.d0
        f_pos_corr2(:,:)=0.d0
        gradients=gradients/nbeads
        displacements=displacements/nbeads
        do j=1,nosc
                do k=1,nosc
                        f_pos_corr1(j,k)=f_pos_corr1(j,k)+gradients(j)*displacements(k) !This 2 are not symmetric in principle
                        f_pos_corr2(j,k)=f_pos_corr2(j,k)+gradients(k)*displacements(j)
                enddo
        enddo

        f_pos_corr=f_pos_corr+(f_pos_corr1+f_pos_corr2)/2 !So we average them








                !First, accumulate the bead average of the xixj product
                av_xx_part(:,:)=0.d0 !Initialize to 0.
                do i=1,nosc
                        do j=1,nosc
                        
                                av_xx_part(i,j)=av_xx_part(i,j)+displacements(i)*displacements(j)
                               
                        enddo
                enddo
                av_xx=av_xx+av_xx_part !Accumulate on the output average xx matrix

                !the averaged individual positions were already calculated as "displacements"
                av_x=av_x+displacements


        end subroutine energies



        subroutine calc_action(model,nbeads,nosc,pos,b,mass,freq_range,coup_const,pol_action,de,req)
                implicit none
                character(len=3),intent(in)::model
                integer,intent(in)::nbeads,nosc
		real*8,intent(in)::b,de,req
		real*8,intent(in),dimension(nosc-1)::coup_const
		real*8,dimension(nosc,nbeads),intent(in)::pos
		real*8,intent(in)::mass
		real*8,dimension(nosc),intent(in)::freq_range
		real*8,intent(inout)::pol_action
		real*8::r2
		real*8,dimension(nosc)::kforce
		real*8,dimension(nosc)::deltar
                integer::i,j

        !Calculate force constant (classical equation):
                kforce(:)=(4.d0*pi**2)*mass*freq_range(:)**2  !kf=4*pi^2*m*w^2

        !Calculating the energy through the direct average value
                pol_action=0.d0
                do i=1,nbeads
                        if (i==nbeads) then
                                deltar(:)=pos(:,i)-pos(:,1)
                        else
                                deltar(:)=pos(:,i)-pos(:,i+1)
                        endif

                        r2=sum(deltar(:)**2)
                        pol_action=pol_action+mass*nbeads/(2.d0*(hbar**2)*b)*r2
                enddo


                do i=1,nbeads
                        pol_action=pol_action+b*calc_pot(model,nosc,pos(:,i),kforce,coup_const,de,req)/nbeads
                enddo

        end subroutine calc_action



        subroutine mc_update(model,nbeads,nosc,pos,step,pol_action,accept_count,b,freq_range,coup_const,mass,freq_coll,de, &
                        req,accepted_moves,nmoves)
                implicit none
                character(len=3),intent(in)::model
                integer,intent(in)::nbeads,nosc
		real*8,dimension(nosc,nbeads),intent(inout)::pos
		real*8,dimension(nosc,nbeads)::postrial
		real*8::trial_action
		real*8,intent(in)::mass,freq_coll,de,req
		real*8,intent(in),dimension(nosc-1)::coup_const
                real*8,dimension(nosc),intent(in)::freq_range
                integer(kind=8),intent(inout)::accept_count
		real*8,intent(inout)::pol_action
		real*8::exp_accept
		real*8,intent(in)::step,b
		real*8,dimension(nosc)::randomnum
		real*8::randomnum2,randomnum3
		real*8,dimension(nosc)::trialmov
		real*8::trial
                integer::rand_bead,rand_osc,i
                integer,dimension(nmoves)::accepted_moves
                integer::move_id,nmoves

                call RANDOM_NUMBER(randomnum3)
                if (randomnum3>(1.d0-freq_coll)) then
                        !Collective movement of beads
                        call RANDOM_NUMBER(randomnum)
                        do i=1,nosc
                                call RANDOM_NUMBER(randomnum2)
                                if (randomnum2>0.5) then
                                        trialmov(i)=randomnum(i)*step
                                else

                                trialmov(i)=-randomnum(i)*step
                                move_id=2
                                endif
                        enddo

                        postrial=pos
                        call RANDOM_NUMBER(randomnum3)
                        rand_osc=int(randomnum3*nosc)+1 !Pick random oscillator
                        do i=1,nbeads
                                postrial(rand_osc,i)=pos(rand_osc,i)+trialmov(rand_osc)
                        enddo

                else
                        !Individual movement of beads
                        call RANDOM_NUMBER(randomnum2)
                        rand_bead=int(randomnum2*nbeads)+1 !Pick random bead

                        call RANDOM_NUMBER(randomnum3)
                        rand_osc=int(randomnum3*nosc)+1 !Pick random oscillator


                        call RANDOM_NUMBER(randomnum2) !Decides step size
                        call RANDOM_NUMBER(randomnum3) !Decides step sign

                        if (randomnum2>0.5d0) then
                                trial=step*randomnum3
                                !Trial will need to be the array containing the path between 2 random beads i,j.
                        else
                                trial=-step*randomnum3
                        endif
                        postrial=pos
                        postrial(rand_osc,rand_bead)=postrial(rand_osc,rand_bead)+trial
                        move_id=1
                endif

                call calc_action(model,nbeads,nosc,postrial,b,mass,freq_range,coup_const,trial_action,de,req)

                if (trial_action<pol_action) then

                        !Accept move
                        pos=postrial
                        pol_action=trial_action
                        accept_count=accept_count+1
                        accepted_moves(move_id)=accepted_moves(move_id)+1


                else
                        call RANDOM_NUMBER(randomnum)
                        exp_accept=exp(-(trial_action-pol_action))

                        if (randomnum(1)<exp_accept) then
                        !Accepts move as not too bad
                                pos=postrial
                                pol_action=trial_action
                                accept_count=accept_count+1
                                accepted_moves(move_id)=accepted_moves(move_id)+1
                        endif
                endif

        end subroutine mc_update

end module mc_mod
