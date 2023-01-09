Module lltb_background_aofht
use lltb_params
use lltb_complex_tools
use lltb_prefuncs
use lltb_precision
use lltb_background_eqs
use lltb_background_ht
use lltb_background_errors
use lltb_background_aofht_funcs
use lltb_background_ht_funcs
use lltb_background_ht_ddok
implicit none

contains


	subroutine lltbbackgroundaofHt(vpi,a,t,sigma)
		implicit none
		type(voidparams) :: vpi, myvpi
		real(dl), optional, intent(out) :: a
		real(dl), optional, intent(in) :: t
		real(dl), optional, intent(out) :: sigma
		real(dl) :: mya, myt, realt, err_est, tBB, aBB, mysigma
		integer :: lastdist, distrepetition, i, signrepetition, lastsign
		logical :: found, beyondturning
		integer, parameter :: imax = 1000
		real(dl) :: atop, abot, dadt, dt, tturn, thist, thisa, lasta, lastt, secondlasta, secondlastt, try_a, atopfix, thisabot, thisatop
!		real(dl), parameter :: changetoNT = 1.
		real(dl), parameter :: changetoNT = 2., switchoffNT = 0.5
!		real(dl), parameter :: accuracy_downscale = 4._dl ! downscale is amount of loss we allow for if i > imax
                ! Moved downscale to lltb_precision.f90
                character(len=512) :: error_msg, error_msg_nr
    integer :: n_fullroots
    real(dl) :: aturn_fullroots


		if(present(t))then
			myt=t
			realt=t
		else
			stop 'lltbbackgroundaofHt called without specifying t.'
		end if


                vpi%crossedbranchcut=.false.
		call void_makeworkxyz(vpi,1._dl)
		call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz))

		call lltb_set_error(vpi)

                if(t<0.)call lltb_error('aofht called with t<0. Error occurred somehwere in tofa before calling this routine.') 


                ! Need to have tturn, and we store it in vpi:
                if(vpi%posrealroots>0)then

                   atop = min(vpi%smallestrealroot,vpi%smallestlinroot)
          aturn_fullroots = atop
          n_fullroots=vpi%posrealroots
                   ! There is a caveat here: previously, the roots
                   ! are calculated not knowing a, i.e. the roots are
                   ! the full roots, not perturbative.
                   ! Now that we put in a, Htofa will recalculate
                   ! whether it should do a perturbative expansion or not
                   ! for the given a. Under this a, the expansion case may change
                   ! and therefore the linear roots might have a 
                   ! turning point BEFORE the turning point that the full
                   ! roots give. Then HTofa will return an error since we ask
                   ! for t(a) with a > amax.
                   ! Prevent trouble by doing this:
                   ! Begin trouble:
                   call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz,atop))
                   atop = min(vpi%smallestrealroot,vpi%smallestlinroot)
                   ! End trouble.
                   ! HOWEVER, this solution should not hide the possible full roots:
                   if(n_fullroots>0 .and. vpi%posrealroots==0)then
                     vpi%posrealroots=n_fullroots
                     atop=aturn_fullroots
                   end if


                   call lltbbackgroundHtofa(vpi,t=tturn,a=atop)
                   vpi%tturn=tturn
                   if(isnan(vpi%tturn))then
                    write(0,*)'received a NaN time from lltbbackgroundHtofa for tturn as a=',atop
                    write(0,*)'for expansion case:',GetOmegaCase(vpi%xyz,atop)
                    write(0,*)'two options for atop:',vpi%smallestrealroot,vpi%smallestlinroot
                    stop
                  end if


                else
                   atop=1._dl
                   vpi%tturn=infinity !1.e30_dl
                   tturn = vpi%tturn

                end if
		call lltb_set_error(vpi)

                ! If tturn < infinity, calculate effective t:
                beyondturning=.false.
		if(vpi%posrealroots>0)then
                   if(myt>tturn)then
                      beyondturning=.true.
                      myt=2._dl*tturn-myt
                      if(myt<0._dl)then
                         write(*,*)'tturn:',tturn
                         write(*,*)'tcrunch:',2._dl*tturn
                         write(*,*)'t:',t
                         write(*,*)'t_eff:',myt
                         call lltb_error('xx t beyond big crunch in lltbbackgroundaofHt')
                      end if
                   end if
                end if


		! First try exact inversion
		try_a = 0._dl ! for later use

		call lltbbackgroundaofHt_TryExact(vpi,mya,myt,mysigma, try_a)

		if(present(a))then
			if(present(Sigma))sigma=mysigma
			a=mya
		else
			stop 'lltbbackgroundaofHT called without giving a: where do you want the answer?'
		end if


		if(mya>0._dl)return
		! if TryExact didn't find a satisfying accuracy, coninue
		! but its returned try_a may still be a fast starting point for the numerical inversion


		abot=0._dl
!		beyondturning=.false.
		if(vpi%posrealroots>0)then
                ! Sometimes the roots are not the same:
                !if(min(vpi%smallestrealroot,vpi%smallestlinroot)<1.e29_dl)then
!			atop = min(vpi%smallestrealroot,vpi%smallestlinroot)
!			call lltbbackgroundHtofa(vpi,t=tturn,a=atop)
!			write(*,*)'tturn:',tturn
!			if(myt>tturn)then
!				beyondturning=.true.
!				myt=2._dl*tturn-myt
!                                if(myt<0._dl)call lltb_error('t beyond big crunch in lltbbackgroundaofHt')
!				write(*,*)realt,myt,tturn
!                                !write(*,*)'Going beyond turning point'
!				call lltb_error('Going beyond turning point. Need to test this.')
!			end if
                        thist = myt
		else
			! set atop now
!			tturn = -1._dl
!			atop = 1._dl
			! test atop:
			thist = -1._dl
			i=0
			do while(thist<myt)
                           !write(*,*)thist,myt,i
				i=i+1
				call lltbbackgroundHtofa(vpi,t=thist,a=atop)
				if(thist<0)then
					if(present(a))then
						a=-1
					end if
					call lltb_set_error(vpi)
					write(*,*)'i =',i
                                        write(*,*)atop
                                        write(*,*)thist
					call lltb_error('t = -1, while trying to set atop.')
					return
				end if
				if (thist < myt)then
					atop = atop * 1.01_dl
					!write(*,*)thist, atop, realt
					!write(*,*) 'need to test this loop: test atop.'
				end if
			end do
		end if

		if((try_a>0._dl) .and. (try_a < atop))then ! Then TryExact found an a which was close, but not close enough, to the answer
			thisa=try_a
!			write(*,*)'Expansion just not good enough: but taking try_a.'
		else	! then TryExact didn't even get near the right result
			thisa = 0.5_dl*(abot+atop)
		end if

                atopfix=atop
		distrepetition = 0
		lastdist = 30
		signrepetition = 0
		lastsign = -1
		i=0
		found = .false.
		secondlasta = 0.
                lastt=thist
                lasta=thisa
		secondlastt= 0.
		do while (.not. found)
			i=i+1
			secondlastt = lastt
			lastt = thist
			call lltbbackgroundHtofa(vpi,t=thist,a=thisa)
	
			dadt = dadt_overH0(vpi%xyz,thisa)

! Comment following when not debugging:
!		err_est = abs(2._dl*(thist - myt)/(thist+myt))
!		writE(*,'(20ES23.14)')atopfix,atop,abot,thisa, thist, myt, dadt,err_est!, 'delta:' , dadt*(myt-thist), int(log10(err_est)), lastdist, vpi%ExpCase,atop-thisa,abot-thisa


			if((dadt /= 0._dl) .and. (thist>0._dl) )then
				secondlasta = lasta
				lasta = thisa
				thisa = thisa + dadt*(myt-thist)
				if(thisa>lasta)then
					abot = lasta
!					if(lastsign == 1)then
!						signrepetition = signrepetition + 1
!					else
!						signrepetition = 0
!						lastsign = 1
!					end if

					! Now force to jump to middle of domain in case of too slow improvement:

					if( ( abs( log10(thist) - log10(myt) ) > changetoNT) .or. (i>int(imax* switchoffNT )) ) thisa = atop + 1

!					if(signrepetition > 10) then
!						thisa = atop + 1
!						signrepetition = 0
!					end if
				else
					atop = lasta
					if( ( abs( log10(thist) - log10(myt) ) > changetoNT) .or. (i>int(imax* switchoffNT )) ) thisa = abot - 1
!					if(lastsign == -1)then
!						signrepetition = signrepetition + 1
!					else
!						signrepetition = 0
!						lastsign = -1
!					end if
!					if(signrepetition > 10) then
!						thisa = abot - 1
!						signrepetition = 0
!					end if
				end if
				if(thisa<abot)then
!					atop = lasta
					thisa = 0.5_dl*(lasta+abot)
					if(abs((thisa-abot)/thisa)<0*accuracy(dl)*acceptance_pure)then ! unless we're caught in our bottom
						thisa = thisa + dadt*(myt-thist)
                                                if(thisa<0._dl)then
                                                   thisa=0._dl
                                                   abot=0._dl
                                                else
                                                   abot = thisa
                                                end if
					end if
				else if (thisa>atop) then
!					abot = lasta
					thisa = 0.5_dl*(lasta + atop)
					if(abs((thisa-atop)/thisa)<0*accuracy(dl)*acceptance_pure)then ! unless we're caught in our bottom
						thisa = thisa + dadt*(myt-thist)
!						atop = thisa
                                                if(thisa>atopfix)then
                                                   thisa=atopfix
                                                   atop=atopfix
                                                else
                                                   atop=thisa
                                                end if
                                                
					end if
				end if
			else
!				call lltb_error('error in aofht')
                           ! if dadt == 0, we are exactly in the turning point amax:
                           abot = thisa
                           thisa = (thisa+atop)*0.5_dl
			end if
			err_est = abs(2._dl*(thist - myt)/(thist+myt))
			if (  err_est < accuracy(dl) * acceptance_pure ) found = .true.
			! or if the precision doesn't increase anymore:
			if ( int(log10(err_est)) == lastdist ) then
				distrepetition = distrepetition + 1
				if(distrepetition>4 .and. err_est<accuracy(dl)*acceptance)found=.true.
!				if(distrepetition>int(imax*0.5))found=.true. ! some cases are too inaccurate.
			else
				distrepetition = 0
			end if
                        if(i>10)then
                           if((lastdist < int(log10(err_est))) .or.  (i>imax-2))then ! If it takes too long, check if error is increasing again
! Why did I take secondlast?
                              if (abs(2._dl*(secondlastt - myt)/(secondlastt+myt))<accuracy(dl)*acceptance*10**(dble(i)/dble(imax)*accuracy_downscale)) then
                                 thisa=secondlasta
                                 found=.true.
! Take both now. Forgot why last was not ok.
                              else if (abs(2._dl*(lastt - myt)/(lastt+myt))<accuracy(dl)*acceptance*10**(dble(i)/dble(imax)*accuracy_downscale)) then
                                 thisa=lasta
                                 found = .true.
                              else if (abs(err_est)<accuracy(dl)*acceptance*10**(dble(i)/dble(imax)*accuracy_downscale)) then ! Or maybe anyway we were somewhat close enough
                                 thisa=lasta
                                 found = .true.
                              end if
                           end if
                        end if

                        ! sometimes we get stuck with wrong bounds due to poor initial step
                        ! <min(accuracy(dl)*acceptance,max(accuracy(dl)*acceptance_pure,abs(thisa/thisatop-1.d0)*0.1d0)) describes a million exceptions..
                        ! - if we are very close to one side of the domain, then we move the boundary by a little bit.
                        ! - but: only if total domain width > 10*distance from nearest boundary
                        ! - unless: total domain is narrower than accuracy(dl)*acceptance_pure = 1.e-14, then we allow change of boundary anyway.
                        if(abs(err_est) > accuracy(dl) * acceptance_pure)then
                           thisatop=atop
                           thisabot=abot
                           if(abs(thisa/thisabot-1.d0)<min(accuracy(dl)*acceptance,max(accuracy(dl)*acceptance_pure,abs(thisa/thisatop-1.d0)*0.1d0)))then
                              abot=abot*0.9d0
                           end if
                           if(abs(thisa/thisatop-1.d0)<min(accuracy(dl)*acceptance,max(accuracy(dl)*acceptance_pure,abs(thisa/thisabot-1.d0)*0.1d0)))then
                              atop=min(1.1d0*atop,atopfix) ! limit from above by atopfix. Don't want H^2<0.
                           end if
                           
                        end if

!if(i>int(dble(imax)*0.9))writE(*,'(20ES23.14)')atopfix,atop,abot,thisa, thist, myt, dadt, dadt*(myt-thist), real(int(log10(err_est))), dble(lastdist), err_est,accuracy(dl)*acceptance*10**(dble(i)/dble(imax)*accuracy_downscale)

			lastdist = int(log10(err_est))


			if(i>imax)then
                           if(all((/abs(atop/atopfix-1._dl),abs(thisa/atopfix-1._dl)/)<1.e-30_dl))then

                              exit
                           end if

                           writE(*,'(20ES23.14)')atopfix,atop,abot,thisa, thist, myt, dadt, dadt*(myt-thist), real(int(log10(err_est))), dble(lastdist), err_est
!writE(*,*)atop,abot,thisa, thist, myt, dadt, 'delta:' , dadt*(myt-thist), int(log10(err_est)), lastdist
				writE(*,'(A,3E30.21)')'xyz:',vpi%xyz
				writE(*,'(A,3E30.21)')'t:',myt
                                error_msg='lltbbackgroundaofHT is hanging'//new_line('A')
                                error_msg=trim(error_msg)//' !x'//new_line('A')
                                write(error_msg_nr,'(3E30.21)')vpi%xyz(1)
                                error_msg=trim(error_msg)//'inline(1) = '//trim(adjustl(error_msg_nr))//'_dl'//new_line('A')
                                error_msg=trim(error_msg)//' !y'//new_line('A')
                                write(error_msg_nr,'(3E30.21)')vpi%xyz(2)
                                error_msg=trim(error_msg)//'inline(2) = '//trim(adjustl(error_msg_nr))//'_dl'//new_line('A')
                                error_msg=trim(error_msg)//' !z'//new_line('A')
                                write(error_msg_nr,'(3E30.21)')vpi%xyz(3)
                                error_msg=trim(error_msg)//'inline(3) = '//trim(adjustl(error_msg_nr))//'_dl'//new_line('A')
                                error_msg=trim(error_msg)//' !a'//new_line('A')
                                write(error_msg_nr,'(3E30.21)')0.
                                error_msg=trim(error_msg)//'inline(4) = '//trim(adjustl(error_msg_nr))//'_dl'//new_line('A')
                                error_msg=trim(error_msg)//' !t'//new_line('A')
                                write(error_msg_nr,'(3E30.21)')myt
                                error_msg=trim(error_msg)//'inline(5) = '//trim(adjustl(error_msg_nr))//'_dl'//new_line('A')
        call lltb_set_error(vpi)                        
!				call lltb_error('lltbbackgroundaofHT is hanging')
				call lltb_error(trim(error_msg))
			end if
                        
		end do

		mya = thisa
		mysigma = myt*err_est*dadt
		if(present(sigma))sigma=mysigma
		if(lfeedback>2)write(*,'(A,I4,A,ES14.5,A,ES14.5)')'a in ',i,' steps, with error ',mysigma, '; a=',mya



		if(present(a))then
			a=mya
		else
			stop 'lltbbackgroundaofHt called without specifying a.'
		end if

!		if(beyondturning)then
!			write(*,*)'myt:',myt
!			write(*,*)'realt:',realt
!			write(*,*)'tturn:',tturn
!			write(*,*)'thist:',thist
!			write(*,*)'mya:',mya
!			stop 'beyond turningpoint.'
!		end if
	end subroutine lltbbackgroundaofHt

	subroutine lltbbackgroundaofHt_TryExact(vpi,a,t,sigma, try_a)
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(out) :: a,sigma,try_a
		real(dl), intent(in) :: t
		!--------------------!
		real(dl) :: epsilons(3), relsize(2), acclim(2),mya,myt
		integer :: i
		character(len=3) :: TSCstr
		real(dl), parameter :: weaken=0.2_dl ! parameter to weaken accuracy tests, in order to use try_a also in less accurate cases.

		a=-1._dl
		try_a=-1._dl
		sigma=1.e30_dl
		myt=t

		call GetSingleEpsilons(vpi,mya,myt,epsilons,weaken)
                ! Restore roots:
                call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz))


		vpi%ExpCaseInv="LLL" ! ExpCaseInv is actually not so relevant..
		TSCstr=""
		Do i = 1, 3
			if(vpi%xyz(i)==0._Dl)then
				vpi%ExpCaseInv(i:i)="Z"
			else if (abs(epsilons(i)/mya) < (accuracy(dl) * rel_acceptance(i))**weaken ) then
				vpi%ExpCaseInv(i:i)="S"
			end if
		end do
		TSCstr = vpi%ExpCaseInv

		! now TSCstr contains the expasion case at hand.
		if(.not.any(TSCstr==(/"LLS","LSS","LSZ","LZS"/)))then
			epsilons=1.e30_Dl
			return
		end if

		relsize(1:2)=abs(epsilons(2:3)/mya)
		acclim(1:2) = (rel_acceptance(2:3)*accuracy(dl))

		if(all(  relsize< acclim**(weaken)) ) then
			try_a = mya - epsilons(2) - epsilons(3)
			if(all(relsize< acclim))then
				a = try_a
				! DO TEST:
				mya=a
				call lltbbackgroundHtofa(vpi,mya,myt)
!				if(abs(myt/t-1._dl)>accuracy(dl)*acceptance_pure)then
				if(abs(myt/t-1._dl)>accuracy(dl))then
					a=-1
				else
!					write(*,*)a,t,myt
					continue ! Success!
				end if
			end if
		end if


	end subroutine lltbbackgroundaofHt_TryExact

	! GetEpsilons: epsilons is change in a, a=a0-epsilon.
	! epst is change in t: t=t0-epst. a(t)=a0(t-epst)=a0(t)-epst dadt (sort of..)
	subroutine GetSingleEpsilons(vpi,mya,myt,epsilons,weaken)
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(out) :: mya
		real(dl), intent(in) :: myt, weaken
		!--------------------!
		real(dl) :: epsilons(3), da0dt, xyz(3), a_r
		complex(dl) :: aoft_c_LSS, aoft_c, ddeltda0, delt_a0, epst
		! The only cases we know how to invert so far, are LLS, LSS, LSZ LZS
		epsilons(1)=1.e30_dl

		! First: assume LSS and get epsilons(2) or conclude 2=L.
		call aoft_m(vpi, aoft_c, myt)
		aoft_c_LSS = aoft_c ! although this one is always real, as omega_m>0 strictly.
		a_r = mreal(aoft_c)
		mya=a_r
                if(mya>vpi%smallestlinroot)then
                   mya=-1._dl
                   epsilons=1.e30_dl
                   return
                end if
		! get these: da0dt, ddeltda0, delta0, all for LSZ
		xyz=vpi%xyz
		xyz(2)=0._dl
                
                call lltb_getallroots(vpi, GetOmegaCase(xyz,a_r))
                if(mya>vpi%smallestlinroot)then
                   mya=-1._dl
                   epsilons=1.e30_dl
                   return
                end if

		da0dt = dadt_overH0(xyz,a_r)

		call tofa_pert_m_k(vpi,a=a_r, tc=delt_a0)

		ddeltda0 = aoft_ddeltda0_m_k(vpi,a=a_r)

		epst = delt_a0 / (1._dl + da0dt * ddeltda0)

		if( abs(maimag(epst))<accuracy(dl))epst= killimag(epst)
		if( abs(mreal(epst))<accuracy(dl))epst= killreal(epst)

		if(  (abs(maimag(epst)/mreal(epst) ) > accuracy(dl)*acceptance_pure)  )then
			epsilons(2)=1.e30_dl
		else
			epsilons(2) = mreal(epst) * da0dt
		end if

!		write(*,*)'epst:',epst
!		call lltb_error('..')
!		stop 'test this.'

		! Second: assume LSS and get epsilons(3), already using the result from epsilons(2)
		if(abs(epsilons(2))<(accuracy(dl)*rel_acceptance(2))**weaken)then ! LSS
			! LSS: can re-use some of previous quantities: a_r and da0dt

	! re-used call aoft_m(vpi, aoft_c, myt)
	! re-used aoft_c_LSS = aoft_c ! although this one is always real, as omega_m>0 strictly.
	! re-used a_r = mreal(aoft_c)

			! get these: da0dt, ddeltda0, delta0, all for LSZ
			xyz=vpi%xyz
			xyz(3)=0._dl

	! re-used da0dt = dadt_overH0(xyz,a_r)

                        call lltb_getallroots(vpi, GetOmegaCase(xyz,a_r))
                        if(a_r>vpi%smallestlinroot)then
                           a_r=-1._dl
                           epsilons=1.e30_dl
                           return
                        end if

			call tofa_pert_m_l(vpi,a=a_r, tc=delt_a0)

			ddeltda0 = aoft_ddeltda0_m_l(vpi,a=a_r)

			epst = delt_a0 / (1._dl + da0dt * ddeltda0)
			if( abs(maimag(epst))<accuracy(dl))epst= killimag(epst)
			if( abs(mreal(epst))<accuracy(dl))epst= killreal(epst)

			if(  (abs(maimag(epst)/mreal(epst) ) > accuracy(dl)*acceptance_pure)  )then
				epsilons(3)=1.e30_dl
			else
				epsilons(3) = mreal(epst) * da0dt
			end if

		else ! LLS -> no inversion I know of.

			epsilons(3)=1.e30_dl

		end if
                
                mya=a_r

	end subroutine GetSingleEpsilons

End Module lltb_background_aofht
