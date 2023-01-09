module lltb_background_ht
use lltb_background_errors
use lltb_params
use lltb_complex_tools
use lltb_prefuncs
use lltb_precision
use lltb_background_ht_funcs
use lltb_background_eqs
implicit none

contains
	subroutine lltbbackgroundHtofa(vpi,a,t) ! more quantities should come in the future
		! calling without t only returns roots
		! I guess eventually this routine will become backgroundoft, calling the inversion routines here.
		implicit none
		type(voidparams) :: vpi, myvpi
		real(dl), intent(in) :: a
		real(dl), optional, intent(out) :: t
		real(dl) :: mya,myt!, xyz(3)
		complex(dl) :: mytc
		character(len=1) :: ThisCase, ThisSubCase(3)
		character(len=3) :: TSCstr
		integer :: i
		integer, parameter :: localfeedback = 0

		call lltb_set_error(vpi)


		if(lfeedback>3)then
			if((vpi%xyz(3)>1.5_dl).and. ( abs(vpi%xyz(2))>1.5_dl ))&
			write(*,*)'From tests it was found that for such large lambda, the code somethimes loses accuracy from 1.e-12 back to 1.e-5'
                end if
		mya = a
		call void_makeworkxyz(vpi,mya)
		vpi%crossedbranchcut = .false.

		mya = vpi%workrescalefac*mya

! all cases:
! -A- {omm, omk, oml} /= 0
! -B- omk = 0, {omm, oml} /= 0
! -C- oml = 0, {omm, omk} /= 0
! -D- omm = 0, {omk, oml} /= 0
! -E- omk /= 0, {omm, oml} = 0
! -F- oml /= 0, {omm, omk} = 0
! -G- omm /= 0, {omk, oml} = 0
! -H- {omm, omk, oml} = 0, duh..

		ThisSubCase = GetOmegaCase(vpi%workxyz,mya)
		TSCstr = ""
		do i = 1, 3
			TSCstr = trim(adjustl(TSCstr))//ThisSubCase(i)
		end do
		vpi%ExpCase = TSCstr

	if( all(ThisSubCase=="L")) then
			ThisCase = 'A'
		elseif( (ThisSubCase(2)/="L") .and. (ThisSubCase(1)=="L") .and. (ThisSubCase(3)=="L") ) then
			ThisCase = 'B'
		elseif( (ThisSubCase(3)/="L") .and. all(ThisSubCase(1:2)=="L") ) then
			ThisCase = 'C'
		elseif( (ThisSubCase(1)/="L") .and. all(ThisSubCase(2:3)=="L") ) then
			ThisCase = 'D'
! Try this, in principle this case should only occur when setting atop, so should not
! end up in actual results (we don't expect Lambda-domination today):
                        ThisCase = 'A'
!                        vpi%uvw_linear = vpi%uvw
                        call CopyUVWtoLinear(vpi)
                        TSCstr = 'LLL'
                        ThisSubCase='L'
                        vpi%ExpCase = TSCstr
! End try
		elseif( (ThisSubCase(2)=="L") .and. (ThisSubCase(1)/="L") .and. (ThisSubCase(3)/="L")  ) then
			ThisCase = 'E'
		elseif( (ThisSubCase(3)=="L") .and. all(ThisSubCase(1:2)/="L") ) then
			ThisCase = 'F'
! Try this, in principle this case should only occur when setting atop, so should not
! end up in actual results (we don't expect Lambda-domination today):
                        ThisCase = 'A'
!                        vpi%uvw_linear = vpi%uvw
                        call CopyUVWtoLinear(vpi)
                        TSCstr = 'LLL'
                        ThisSubCase='L'
                        vpi%ExpCase = TSCstr
! End try
		elseif( (ThisSubCase(1)=="L") .and. all(ThisSubCase(2:3)/="L") ) then
			ThisCase = 'G'
		elseif( all(ThisSubCase/="L") ) then
			call lltb_error( 'everything zero? Impossible.. Case L in lltb_background_ht.' )
		end if

		!if(localfeedback>0)write(*,*)'pre uvw=',vpi%uvw
		!if(localfeedback>0)write(*,*)'pre uvw_linear=',vpi%uvw_linear

		call lltb_getallroots(vpi, ThisSubCase)
		call lltb_set_error(vpi)

		!if(localfeedback>0)write(*,*)'post uvw=',vpi%uvw
		!if(localfeedback>0)write(*,*)'post uvw_linear=',vpi%uvw_linear

		if(ThisCase=='A') then
			! later add a'(r,t) here.

			! This case works for A&B (also for k=0, \int dt is elliptic and encoded in void_uvw).
!			call void_uvw(vpi%workxyz,vpi%uvw, vpi%posrealroots, vpi%smallestrealroot)
!			call lltb_set_error(vpi)
			call tofa_mkl(vpi,mya, tc=mytc, t=myt)

		else if (ThisCase=='B') then

			call tofa_mkl(vpi,mya, tc=MYtc, t=myt)

		else if (ThisCase=='C') then
			! this case has one root: a= -x/y, which is positive for -y/x>0, i.e. omm>1.

			call tofa_km(vpi,mya,tc=mytc,t=myt)

		else if (ThisCase=='D') then

			call tofa_kl(vpi,mya,tc=mytc,t=myt)

		else if (ThisCase=='E') then
			! this case has no roots
			myt=mya/sqrt(abs(vpi%workxyz(2))) ! Since this should only occur if omk = 1 - omm - oml = 1
			mytc=dcmplx(myt,0._dl)
		else if (ThisCase=='F') then
			! this case has no big bang
			mytc=0.
			myt=0.
		else if (ThisCase=='G') then
			! this case has no roots
			myt=2._dl/3._dl*mya**(1.5_dl)/sqrt(vpi%workxyz(1))
			mytc=dcmplx(myt,0._dl)
		else
!			stop 'not coded yet'
			if(lfeedback > 2) write(*,*) 'no minkowsky!'
			mytc=dcmplx(-1.,0.)
			myt=-1._dl
		end if

		! roots are set now: reset error with these roots
!		call lltb_set_error(vpi)

		! C: tofa_km
		! D: tofa_kl
		! E: tofa_k
		! F: tofa_l
		! G: tofa_m
		! H: stop it, moron.


		! and now the expansions

		if (any(ThisSubCase=="S"))then
			if(TSCstr=="LSL")then ! B
				call tofa_pert_ml_k(vpi,mya, tc=mytc, t=myt)
				if ( isnan(myt) ) then ! Expansion blew up: real elliptic integral is more accurate.
!					call void_uvw(vpi%workxyz,vpi%uvw, vpi%posrealroots,vpi%smallestrealroot)
					!vpi%uvw_linear = vpi%uvw
                                        call CopyUVWtoLinear(vpi)
          vpi%ExpCase = "LLL"
          vpi%ExpCaseInv = "LLL"
					call tofa_mkl(vpi,mya, tc=mytc, t=myt)
!          write(0,*)'debugging right now: got nan for t at a. Calling H to see if we are beyond turning point.'
!          myt=HoverH0(vpi%xyz,mya)
!          stop 'debugging lltbbackgroundHtofa'
          if(isnan(myt)) then ! still..
            call lltb_error( 'isnan(myt) for case LSL in lltb_background_ht.f90.' )
          end if
				end if
			elseif(TSCstr=="LLS")then ! C
				if(ThisCase/="C")call lltb_error( 'Wrong cases for LLS' )
				call tofa_pert_mk_l(vpi,mya, tc=mytc, t=myt)
				if ( real(mytc) < -0.99_dl) then ! Expansion blew up: real elliptic integral is more accurate.
!					call void_uvw(vpi%workxyz,vpi%uvw, vpi%posrealroots,vpi%smallestrealroot)
					!vpi%uvw_linear = vpi%uvw
                                        call CopyUVWtoLinear(vpi)
          vpi%ExpCase = "LLL"
          vpi%ExpCaseInv = "LLL"
					call tofa_mkl(vpi,mya, tc=mytc, t=myt)
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
				call tofa_pert_m_l(vpi,mya, tc=mytc, t=myt)
				call tofa_pert_m_k(vpi,mya, tc=mytc, t=myt)

			! 6 SZ cases: LSZ, ZLS, SZL, LZS, SLZ, ZSL

			elseif(TSCstr=="LSZ")then ! G
				if(ThisCase/="G")call lltb_error( 'Wrong cases for LSZ' )
!				stop 'not coded yet: LSS'
!				myt=-1._dl
				call tofa_pert_m_k(vpi,mya, tc=mytc, t=myt)
			elseif(TSCstr=="LZS")then ! G
				if(ThisCase/="G")call lltb_error( 'Wrong cases for LZS' )
!				stop 'not coded yet: LSS'
!				myt=-1._dl
				call tofa_pert_m_l(vpi,mya, tc=mytc, t=myt)

! and the next four are only hypothetically: we never haver omm< 1.d-...
			elseif(TSCstr=="SLZ")then ! G
				if(ThisCase/="E")call lltb_error( 'Wrong cases for SLZ' )
write(*,*)'For a:',mya
				call lltb_error( 'not coded yet: SLZ' )
!				call tofa_pert_k_m(vpi,mya, tc=mytc, t=myt)
			elseif(TSCstr=="ZLS")then ! G
				if(ThisCase/="E")call lltb_error( 'Wrong cases for ZLS' )
				call lltb_error( 'not coded yet: ZLS' )
!				call tofa_pert_k_l(vpi,mya, tc=mytc, t=myt)

			elseif(TSCstr=="SZL")then ! G
				if(ThisCase/="F")call lltb_error( 'Wrong cases for SZL')
				call lltb_error( 'not coded yet: SLZ' )
!				call tofa_pert_l_m(vpi,mya, tc=mytc, t=myt)
			elseif(TSCstr=="ZSL")then ! G
				if(ThisCase/="F")call lltb_error( 'Wrong cases for ZSL')
				call lltb_error( 'not coded yet: ZLS' )
!				call tofa_pert_l_k(vpi,mya, tc=mytc, t=myt)

			else
				write(*,*)'This unknown case: ', TSCstr
				call lltb_error( 'Do not know what to do with this case.' )
			end if

		end if



		mytc=mytc*sqrt(vpi%workxyzfac)
if(isnan(myt))then
  write(0,*)'before normalization we have nan(t) in lltbbackgroundHtofa. a:',mya
  write(0,*)'expansion case:',TSCstr
  call lltb_error('isnan(t)')
end if
		myt=myt*sqrt(vpi%workxyzfac)
if(isnan(myt))then
  write(0,*)'after normalization we have nan(t) in lltbbackgroundHtofa. a, norm:',mya,sqrt(vpi%workxyzfac)
stop
end if

!			if(  (mya>vpi%smallestrealroot) )then
		if(  (mya>vpi%smallestlinroot) )then
			if(lfeedback>1)then
				write(*,*)'Scale factor larger than turning point: unphysical for our purposes.',maimag(vpi%uvw(1)),accuracy(dl) * acceptance_pure
				write(*,*)'uvw:',vpi%uvw
				writE(*,*)'realroot:',vpi%smallestrealroot
				writE(*,*)'linroot:',vpi%smallestlinroot
				write(*,*)'a',mya
                                write(*,*)'linroot-a',vpi%smallestlinroot-mya
			end if
			myt=-2._dl
		end if

		if(myt>0.)then

			if(abs(maimag(mytc))/abs(mreal(mytc))>accuracy(dl)*acceptance)then
!				if(abs(maimag(mytc))/abs(mreal(mytc)) < abs(get_nearest_cr(vpi%uvw,mya)/mya-1._dl) ) then
				if(  (abs(maimag(mytc))/abs(mreal(mytc)) < abs(vpi%smallestrealroot/mya-1._dl)) .or. & ! near pole here, next line on pole
				  (  (abs(vpi%smallestrealroot/mya-1._dl)<accuracy(dl)) .and. (abs(maimag(mytc))/abs(mreal(mytc)) < accuracy(dl)*acceptance_pole)  )&
				 ) then
					! then it's ok! We know it's difficult near the poles.
					if(lfeedback > 3) then
						write(*,*)'Accepting this imaginary time:'
						write(*,*)'Im(t)/Re(t):',abs(maimag(mytc))/abs(mreal(mytc)),' > ',accuracy(dl)*acceptance
						write(*,*)'Re(t):',abs(mreal(mytc))
						write(*,*)'xyz:',vpi%xyz
						write(*,*)'uvw:',vpi%uvw
						write(*,*)'a:',mya
						write(*,*)'with allowed imaginary part:',abs(get_nearest_cr(vpi%uvw,mya)/mya-1._dl)
					end if
					continue
				else if(abs(mreal(vpi%uvw(2)-vpi%uvw(3)))>accuracy(dl))then
					if(lfeedback > 2) then
						write(*,*)'Imaginary time with probably three positive roots: scale factor out of bounds.'
						write(*,*)'Im(t)/Re(t):',abs(maimag(mytc))/abs(mreal(mytc)),' > ',accuracy(dl)*acceptance
						write(*,*)'Re(t):',abs(mreal(mytc))
						write(*,*)'xyz:',vpi%xyz
						write(*,*)'uvw:',vpi%uvw
						write(*,*)'a:',mya
						write(*,*)'distance from pole:',vpi%smallestrealroot/mya-1._dl
						write(*,*)'pole:',vpi%smallestrealroot
					end if
					myt=-2._dl
				else if(mreal(vpi%uvw(1))>0 .and. maimag(vpi%uvw(1))<accuracy(dl)*acceptance)then
					if(lfeedback > 2) then
						write(*,*)'Imaginary time with positive real root: scale factor out of bounds.'
						write(*,*)'Im(t)/Re(t):',abs(maimag(mytc))/abs(mreal(mytc)),' > ',accuracy(dl)*acceptance
						write(*,*)'Re(t):',abs(mreal(mytc))
						write(*,*)'xyz:',vpi%xyz
						write(*,*)'uvw:',vpi%uvw
						write(*,*)'a:',mya
					end if
					myt=-2._dl
				else
					write(*,*)vpi%uvw
					call void_uvw(vpi%xyz,vpi%uvw)!, vpi%posrealroots, vpi%smallestrealroot)
					write(*,*)vpi%uvw
					write(*,*)vpi%xyz,a
					write(*,*)vpi%workxyz,mya
					call lltb_error( 'imaginary t' )
				end if
			end if
		end if

		if(present(t))t=myt


	end subroutine lltbbackgroundHtofa



	subroutine MyLittleDebugger(arr,arri,arrc,write)
		real(dl), optional, intent(in) :: arr(3)
		integer, optional :: arri(3)
		character(len=1), optional :: arrc(3)
		real(dl), save :: arrs(3)
		integer, save :: arris(3)
		character(len=1), save :: arrcs(3)
		logical, optional :: write

		if(.not.present(write)) then
			arrs=arr
			arris=arri
			arrcs=arrc
			return
		else
			write(*,*)arrs
			write(*,*)arrcs
			write(*,*)arris
			write(*,*)sqrt(arrs(2)*arrs(3))
		end if

	end subroutine MyLittleDebugger

end module lltb_background_ht
