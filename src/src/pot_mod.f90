module pot_mod
        implicit none

        contains

	real*8 function calc_pot(model,nosc,pos_osc,kforce,coup_const,de,req)
                implicit none
                character(len=3),intent(in)::model
                integer,intent(in)::nosc
		real*8,dimension(nosc)::pos_osc
                real*8,intent(in)::de,req !Dissoc energy and eq distance for morse potential.
		real*8,intent(in),dimension(nosc-1)::coup_const
                real*8,dimension(nosc),intent(in)::kforce
		real*8::energy,aux,a
                integer::i

                if (model=='har') then
                        energy=0.5d0*kforce(1)*pos_osc(1)**2
                        do i=2,nosc
                                energy=energy+0.5d0*kforce(i)*pos_osc(i)**2+coup_const(i-1)*pos_osc(i)*pos_osc(1) !Simplified interaction hamiltonian
                        enddo


                else if (model=='fhr') then!Harmonic oscillator including full interaction in the hamiltonian
                        energy=0.d0
                        aux=0.d0
                        do i=2,nosc
                                energy=energy+coup_const(i-1)/kforce(1)*pos_osc(i)
                                aux=aux+0.5d0*kforce(i)*pos_osc(i)**2 !This is the dip moment term
                        enddo
                        energy=0.5d0*kforce(1)*(pos_osc(1)+energy)**2+aux



                else if (model=='fh2') then!Harmonic oscillator including full interaction in the hamiltonian
                        energy=0.d0
                        aux=0.d0
                        do i=2,nosc
                                energy=energy+coup_const(i-1)/kforce(1)*pos_osc(i) !This is the dip moment squared parenthesis.
                                aux=aux+0.5d0*kforce(i)*(pos_osc(i)*(-1)**i)**2
                        enddo
                        energy=0.5d0*kforce(1)*(pos_osc(1)+energy)**2+aux



                else if (model=='mor') then
                        energy=0.5d0*kforce(1)*pos_osc(1)**2 !Cavity energy contribution. Always harmonic
                        
!                        write(6,*) 'de',de,'alpha',a,'k',kforce,'E',energy
                        do i=2,nosc
                                a=dsqrt(kforce(i)/(2.d0*de))
                                energy=energy+de*(1.d0-dexp(-a*(pos_osc(i)-req)))**2 !Adding morse part
!                                       write(6,*) 'energy subroutine',i,energy
                                energy=energy+coup_const(i-1)*pos_osc(i)*pos_osc(1) !Adding interaction with the cavity
                        enddo

                else if (model=='fmr') then
                        aux=0.d0
                        energy=0.d0
                        do i=2,nosc
                                a=dsqrt(kforce(i)/(2.d0*de))
                                aux=aux+de*(1.d0-dexp(-a*(pos_osc(i)-req)))**2
                                energy=energy+coup_const(i-1)/kforce(1)*pos_osc(i)
                        enddo
                        energy=0.5d0*kforce(1)*(pos_osc(1)+energy)**2+aux


                else if (model=='fm2') then
                        aux=0.d0
                        energy=0.d0
                        do i=2,nosc
                                a=dsqrt(kforce(i)/(2.d0*de))
                                aux=aux+de*(1.d0-dexp(-a*(pos_osc(i)-req)*(-1)**i))**2
                                energy=energy+coup_const(i-1)/kforce(1)*pos_osc(i)
                        enddo
                        energy=0.5d0*kforce(1)*(pos_osc(1)+energy)**2+aux



                else if (model=='dmr') then !double well by concatenating morses
                        aux=0.d0
                        energy=0.d0
                        do i=2,nosc
                                a=dsqrt(kforce(i)/(2.d0*de))
                                if (pos_osc(i)>0.d0) then
                                        aux=aux+de*(1.d0-dexp(a*(pos_osc(i)-req)))**2
                                        energy=energy+coup_const(i-1)/kforce(i)*pos_osc(i)
                                else

                                        aux=aux+de*(1.d0-dexp(-a*(pos_osc(i)+req)))**2
                                        energy=energy+coup_const(i-1)/kforce(i)*pos_osc(i)
                                endif
                        enddo
                        energy=0.5d0*kforce(1)*(pos_osc(1)+energy)**2+aux



                else
                        write(6,*) 'Unknown model specification. Killing execution'
                        STOP

                endif

                calc_pot=energy
                return
        end function calc_pot



        subroutine calc_grad(model,nosc,pos_osc,kforce,coup_const,de,req,grad)
                implicit none
                character(len=3),intent(in)::model
                integer,intent(in)::nosc
		real*8,dimension(nosc),intent(in)::pos_osc
		real*8,intent(in)::de,req
		real*8,intent(in),dimension(nosc-1)::coup_const
		real*8,dimension(nosc),intent(in)::kforce
		real*8,dimension(nosc),intent(inout)::grad
		real*8::sumpos,a
                integer::i

                if (model=='har') then
                                grad(:)=kforce(:)*pos_osc(:) !Harmonic part

                                !Coupling part
                                sumpos=0.d0
                                do i=2,nosc
                                        grad(i)=grad(i)+coup_const(i-1)*pos_osc(1)
                                        sumpos=coup_const(i-1)*pos_osc(i)
                                enddo
                                grad(1)=grad(1)+sumpos


                else if (model=='fhr') then !Harmonic oscillator including full interaction in the hamiltonian
                        grad(:)=kforce(:)*pos_osc(:) !Harmonic part
                        sumpos=0.d0
                        do i=2,nosc
                                grad(i)=grad(i)+coup_const(i-1)*pos_osc(1)
                                sumpos=sumpos+pos_osc(i)*coup_const(i-1)
                        enddo
                        grad(2:)=grad(2:)+coup_const(:)/kforce(1)*sumpos
                        grad(1)=grad(1)+sumpos


                else if (model=='fh2') then !Harmonic oscillator including full interaction in the hamiltonian
                        grad(:)=0.d0
                        grad(1)=kforce(1)*pos_osc(1) !Harmonic potential part
                        sumpos=0.d0
                        do i=2,nosc
                                grad(i)=coup_const(i-1)*pos_osc(1)+kforce(i)*pos_osc(i)
                                sumpos=sumpos+pos_osc(i)*coup_const(i-1)
                        enddo
                        grad(2:)=grad(2:)+coup_const(:)/kforce(1)*sumpos
                        grad(1)=grad(1)+sumpos




                else if (model=='mor') then
                        grad(1)=kforce(1)*pos_osc(1)
                        do i=2,nosc
                                a=dsqrt(kforce(i)/(2.d0*de))
                                grad(i)=2.d0*de*a*dexp(-a*(pos_osc(i)-req))*(1.d0-dexp(-a*(pos_osc(i)-req)))  !Morse part
                        enddo
                        sumpos=0.d0
                        do i=2,nosc ! Add interaction to cavity
                                grad(i)=grad(i)+coup_const(i-1)*pos_osc(1)
                                sumpos=sumpos+pos_osc(i)*coup_const(i-1)
                        enddo

                        grad(1)=grad(1)+sumpos

                else if (model=='fmr') then
                        grad(:)=0.d0
                        sumpos=0.d0
                        grad(1)=kforce(1)*pos_osc(1)

                        do i=2,nosc
                                a=dsqrt(kforce(i)/(2.d0*de))
                                grad(i)=2.d0*a*de*dexp(-a*(pos_osc(i)-req))*(1-dexp(-a*(pos_osc(i)-req)))
                        enddo

                        do i=2,nosc
                                sumpos=sumpos+pos_osc(i)*coup_const(i-1)
                                grad(i)=grad(i)+coup_const(i-1)*pos_osc(1)
                        enddo

                        do i=2,nosc
                                grad(i)=grad(i)+(coup_const(i-1)/kforce(1))*sumpos
                        enddo

                        grad(1)=grad(1)+sumpos



                else if (model=='fm2') then
                        grad(:)=0.d0
                        sumpos=0.d0
                        grad(1)=kforce(1)*pos_osc(1)

                        do i=2,nosc
                                a=dsqrt(kforce(i)/(2.d0*de))
                                grad(i)=2.d0*((-1)**i)*a*de*dexp(-a*(pos_osc(i)-req)*(-1)**i)*(1-dexp(-a*(pos_osc(i)-req)*(-1)**i))
                        enddo

                        do i=2,nosc
                                sumpos=sumpos+pos_osc(i)*coup_const(i-1) !This may make the total dipole moment equal to 0.
                                grad(i)=grad(i)+coup_const(i-1)*pos_osc(1)
                        enddo


                        do i=2,nosc
                                grad(i)=grad(i)+(coup_const(i-1)/kforce(1))*sumpos
                        enddo

                        grad(1)=grad(1)+sumpos
                        






                else if (model=='dmr') then
                        grad(:)=0.d0
                        sumpos=0.d0
                        grad(1)=kforce(1)*pos_osc(1)

                        do i=2,nosc
                                if (pos_osc(i)>0.d0) then
                                        a=dsqrt(kforce(i)/(2.d0*de))
                                        grad(i)=2.d0*a*de*dexp(a*(pos_osc(i)-req))*(1-dexp(a*(pos_osc(i)-req)))
                                else
                                        a=dsqrt(kforce(i)/(2.d0*de))
                                        grad(i)=2.d0*a*de*dexp(-a*(pos_osc(i)+req))*(1-dexp(-a*(pos_osc(i)+req)))
                                endif
                        enddo


                        do i=2,nosc
                                sumpos=sumpos+pos_osc(i)*coup_const(i-1)
                                grad(i)=grad(i)+coup_const(i-1)*pos_osc(1)
                        enddo

                        do i=2,nosc
                                grad(i)=grad(i)+(coup_const(i-1)/kforce(i))*sumpos
                        enddo

                        grad(1)=grad(1)+sumpos



                else
                        write(6,*) 'Unknown model specification. Killing execution'
                        STOP

                endif
                return
        end subroutine calc_grad



end module pot_mod
