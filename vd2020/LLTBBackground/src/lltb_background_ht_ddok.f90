module lltb_background_ht_ddok
use lltb_background_errors
use lltb_params
use lltb_complex_tools
use lltb_prefuncs
use lltb_precision
use lltb_background_ht_ddok_funcs
implicit none


contains

        ! ddok_tofa() returns dt/dy of the integral \int_0^A \sqrt(a / (X + Y a + Z a^3)) da
        ! assuming always the positive root. Beware! This means, assuming always expansion, 
        ! not contraction. So, when calling this routing, apply your own minus sign to the result
        ! when you are considering a contracting phase.
        ! 
        ! input: vpi, a
        ! output: dt/dy is stored in tc (complex) and t (real).
	subroutine ddok_tofa(vpi,a, tc, t)
		! Pretty much the same as lltbbackgroundhtofa
		implicit none
		type(voidparams) :: vpi
!		real(dl), intent(in) :: xyz(3)
!		complex(dl), intent(in) :: uvw(3)
		real(dl), intent(in) :: a
		real(dl), intent(out), optional :: t
		complex(dl), intent(out), optional :: tc
		real(dl) :: mya,myt!, xyz(3)
		complex(dl) :: mytc
		character(len=1) :: ThisCase, ThisSubCase(3)
		character(len=3) :: TSCstr
		integer :: i
                real(dl) :: signmem

		call lltb_set_error(vpi)

                mya = a


! all cases:
! -A- {omm, omk, oml} /= 0
! -B- omk = 0, {omm, oml} /= 0
! -C- oml = 0, {omm, omk} /= 0
! -D- omm = 0, {omk, oml} /= 0
! -E- omk /= 0, {omm, oml} = 0
! -F- oml /= 0, {omm, omk} = 0
! -G- omm /= 0, {omk, oml} = 0
! -H- {omm, omk, oml} = 0, duh..

!		ThisSubCase = GetOmegaCase(vpi%workxyz,mya)
!		TSCstr = ""
!		do i = 1, 3
!			TSCstr = trim(adjustl(TSCstr))//ThisSubCase(i)
!		end do
!		vpi%ExpCase = TSCstr
		TSCstr = vpi%ExpCase
		do i = 1, 3
			ThisSubCase(i) = TSCstr(i:i)
		end do

	if( all(ThisSubCase=="L")) then
			ThisCase = 'A'
		elseif( (ThisSubCase(2)/="L") .and. (ThisSubCase(1)=="L") .and. (ThisSubCase(3)=="L") ) then
			ThisCase = 'B'
		elseif( (ThisSubCase(3)/="L") .and. all(ThisSubCase(1:2)=="L") ) then
			ThisCase = 'C'
		elseif( (ThisSubCase(1)/="L") .and. all(ThisSubCase(2:3)=="L") ) then
			ThisCase = 'D'
		elseif( (ThisSubCase(2)=="L") .and. (ThisSubCase(1)/="L") .and. (ThisSubCase(3)/="L")  ) then
			ThisCase = 'E'
		elseif( (ThisSubCase(3)=="L") .and. all(ThisSubCase(1:2)/="L") ) then
			ThisCase = 'F'
		elseif( (ThisSubCase(1)=="L") .and. all(ThisSubCase(2:3)/="L") ) then
			ThisCase = 'G'
		elseif( all(ThisSubCase/="L") ) then
			call lltb_error( 'everything zero? Impossible.. Case L in lltb_background_ht.' )
		end if

!		call lltb_getallroots(vpi, ThisSubCase)
!		call lltb_set_error(vpi)
if(lfeedback>3)then
   write(*,*)ThisCase
   write(*,*)ThisSubCase
   write(*,*)vpi%ExpCase
end if

		if(ThisCase=='A') then
			! later add a'(r,t) here.

			! This case works for A&B (also for k=0, \int dt is elliptic and encoded in void_uvw).
!			call void_uvw(vpi%workxyz,vpi%uvw, vpi%posrealroots, vpi%smallestrealroot)
!			call lltb_set_error(vpi)
			call ddok_tofa_mkl(vpi,mya, tc=mytc, t=myt)

		else if (ThisCase=='B') then

			call ddok_tofa_mkl(vpi,mya, tc=MYtc, t=myt) ! since mkl uses uvw_linear.

		else if (ThisCase=='C') then
			! this case has one root: a= -x/y, which is positive for -y/x>0, i.e. omm>1.

			call ddok_tofa_mk(vpi,mya,tc=mytc,t=myt)

		else if (ThisCase=='D') then

			call lltb_error('Case D in ddok_tofa: should not occur.')
!			call ddok_tofa_kl(vpi,mya,tc=mytc,t=myt)

		else if (ThisCase=='E') then

			! this case has no roots
			call lltb_error('Case E in ddok_tofa: should not occur.')
!                        call ddok_tofa_k(vpi,mya,tc=mytc,t=myt)

		else if (ThisCase=='F') then
			! this case has no big bang
			mytc=0.
			myt=0.
		else if (ThisCase=='G') then
			! this case has no roots
			call ddok_tofa_m(vpi,mya,tc=mytc,t=myt)
		else
!			stop 'not coded yet'
			if(lfeedback > 2) write(*,*) 'no minkowsky!'
			mytc=dcmplx(-1.,0.)
			myt=-1._dl
		end if

		! roots are set now: reset error with these roots
!		call lltb_set_error(vpi)

		! C: ddok_tofa_km
		! D: ddok_tofa_kl
		! E: ddok_tofa_k
		! F: ddok_tofa_l
		! G: ddok_tofa_m
		! H: stop it, moron.


		! and now the expansions

		if (any(ThisSubCase=="S"))then
			if(TSCstr=="LSL")then ! B
                           !call lltb_error('Asking for ddok_pert_ml_k, but do not have it.')
                           ! proposed solution to this difficulty, seems ok:
                           ! vpi%uvw_linear = vpi%uvw
                           call CopyUVWtoLinear(vpi)
                           
                           ! ugly sign trouble hack:
                           signmem = myt/abs(myt)
                           
                           call ddok_tofa_mkl(vpi,mya, tc=mytc, t=myt)
                           
                           ! ugly sign trouble hack:
                           if(myt/abs(myt)*signmem<0)then
                              mytc=mytc*signmem
                              myt=myt*signmem
                           end if
                           
!--				call ddok_pert_ml_k(vpi,mya, tc=mytc, t=myt)
			elseif(TSCstr=="LLS")then ! C
				if(ThisCase/="C")call lltb_error( 'Wrong cases for LLS' )
        call ddok_pert_mk_l(vpi,mya, tc=mytc, t=myt)
				if ( real(mytc) < -0.99_dl) then ! Expansion blew up: real elliptic integral is more accurate.
!					call void_uvw(vpi%workxyz,vpi%uvw, vpi%posrealroots,vpi%smallestrealroot)
                                        !vpi%uvw_linear = vpi%uvw
                                        call CopyUVWtoLinear(vpi)
					call ddok_tofa_mkl(vpi,mya, tc=mytc, t=myt)
				end if
			elseif(TSCstr=="SLL")then ! D
				if(ThisCase/="D")call lltb_error( 'Wrong cases for SLL')
				! SLL is non-perturbative in omega_m. Fortunately this case should not occur in physical configurations
				! -- because ZLL has no big bang
				if(lfeedback>2)write(*,*)'Non-perturbative case: SLL. Should not occur in viable cosmologies.'
				myt=-1._dl
			elseif(TSCstr=="SLS")then ! E
				if(ThisCase/="E")call lltb_error( 'Wrong cases for SLS')
!				stop 'not coded yet: SLS'
				if(lfeedback>2)write(*,*)'Non-perturbative case: SLS. Should not occur in viable cosmologies.'
				! SLS is non-perturbative. Fortunately this case should not occur in physical configurations
				! -- because ZLZ has no big bang
				myt=-1._dl
			elseif(TSCstr=="SSL")then ! F
				if(ThisCase/="F")call lltb_error( 'Wrong cases for SSL' )
				if(lfeedback>2)write(*,*)'Non-perturbative case: SSL. Should not occur in viable cosmologies.'
				! SLL is non-perturbative. Fortunately this case should not occur in physical configurations
				! -- because ZZL has no big bang
				myt=-1._dl
			elseif(TSCstr=="LSS")then ! G
				if(ThisCase/="G")call lltb_error( 'Wrong cases for LSS' )
!				stop 'not coded yet: LSS'
!				myt=-1._dl
				call ddok_pert_m_l(vpi,mya, tc=mytc, t=myt)
				call ddok_pert_m_k(vpi,mya, tc=mytc, t=myt)

			! 6 SZ cases: LSZ, ZLS, SZL, LZS, SLZ, ZSL

			elseif(TSCstr=="LSZ")then ! G
				if(ThisCase/="G")call lltb_error( 'Wrong cases for LSZ' )
!				stop 'not coded yet: LSS'
!				myt=-1._dl
				call ddok_pert_m_k(vpi,mya, tc=mytc, t=myt)
			elseif(TSCstr=="LZS")then ! G
				if(ThisCase/="G")call lltb_error( 'Wrong cases for LZS' )
!				stop 'not coded yet: LSS'
!				myt=-1._dl
				call ddok_pert_m_l(vpi,mya, tc=mytc, t=myt)

! and the next four are only hypothetically: we never haver omm< 1.d-...
			elseif(TSCstr=="SLZ")then ! G
				if(ThisCase/="E")call lltb_error( 'Wrong cases for SLZ' )
				call lltb_error( 'not coded yet: SLZ' )
!				call ddok_pert_k_m(vpi,mya, tc=mytc, t=myt)
			elseif(TSCstr=="ZLS")then ! G
				if(ThisCase/="E")call lltb_error( 'Wrong cases for ZLS' )
				call lltb_error( 'not coded yet: ZLS' )
!				call ddok_pert_k_l(vpi,mya, tc=mytc, t=myt)

			elseif(TSCstr=="SZL")then ! G
				if(ThisCase/="F")call lltb_error( 'Wrong cases for SZL')
				call lltb_error( 'not coded yet: SLZ' )
!				call ddok_pert_l_m(vpi,mya, tc=mytc, t=myt)
			elseif(TSCstr=="ZSL")then ! G
				if(ThisCase/="F")call lltb_error( 'Wrong cases for ZSL')
				call lltb_error( 'not coded yet: ZLS' )
!				call ddok_pert_l_k(vpi,mya, tc=mytc, t=myt)

			else
				write(*,*)'This unknown case: ', TSCstr
				call lltb_error( 'Do not know what to do with this case.' )
			end if

		end if



		mytc=mytc*sqrt(vpi%workxyzfac)
		myt=myt*sqrt(vpi%workxyzfac)


		if(present(tc))tc=mytc
		if(present(t))t=myt


	return

	end subroutine ddok_tofa

end module lltb_background_ht_ddok
