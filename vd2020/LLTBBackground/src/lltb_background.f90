module lltb
use elliptics
use lltb_precision
use lltb_params
!use lltb_prefuncs
!use lltb_background_ht
use lltb_background_ht_ddok
use lltb_background_ht_ddok_funcs
!use lltb_background_ht_funcs
use lltb_background_aofht
!use lltb_background_aofht_funcs
use lltb_background_errors
use lltb_background_eqs
use lltb_background_limits
!use lltb_complex_tools

Interface lltb_precisiongoal
   module procedure lltb_precisiongoal_sl, lltb_precisiongoal_dl
end interface

private

public lltb_functions, lltb_functions_a, lltb_normalization, lltb_feedback, lltb_precisiongoal, lltb_attach_to_error, lltb_error, lltb_ricci

contains

	! lltb_background: returns all possible (optional) functions describing the
	! metric, as a function of: H0(r=infinity, t=t0), Omega_Lambda(r_inf, t0)
	! k(r) (subroutine k(r)), dkdr (subroutine dkdr(r),
	! tBB(r) (subroutine tBB(r)), dtBBdr (subroutine dtBBdr(r),
	! normalization Mtilde
	! coordinates r and t
	! Returning optionally:
	! R(r,t)        - Rltb
	! R'(r,t)       - Rp
	! \dot R'(r,t)  - Rpd
	! \ddot R'(r,t) - Rpdd
	! S(r,t)        - S
	! \dot S(r,t)   - Sd
	! \ddot S(r,t)  - Sdd
	! a(r,t)        - a
	! a'(r,t)       - ap
	! \dot a'(r,t)  - apd
	! \ddot a(r,t)  - add
	! \ddot a'(r,t) - apdd
	! H(r,t)        - H
        ! tturn         - age of the universe for which this radius experiences H=0
	subroutine lltb_functions_all(H0_inf, Lambda, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp, tturn)
		implicit none
		! ---------------------------------------------------------- !
		! INPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(in) :: H0_inf, Lambda
		interface
			function kofr(rr)
				use lltb_precision
				real(dl) :: kofr
				real(dl), intent(in) :: rr
			end function kofr
		end interface
		interface
			function dkdr(rr)
				use lltb_precision
				real(dl) :: dkdr
				real(dl), intent(in) :: rr
			end function dkdr
		end interface
		interface
			function tbbofr(rr)
				use lltb_precision
				real(dl) :: tbbofr
				real(dl), intent(in) :: rr
			end function tbbofr
		end interface
		interface
			function dtbbdr(rr)
				use lltb_precision
				real(dl) :: dtbbdr
				real(dl), intent(in) :: rr
			end function dtbbdr
		end interface
		real(dl), intent(in) :: r, t, Mtilde
		! ---------------------------------------------------------- !
		! OUTPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(out), optional :: Rltb, Rp, Rpd, S, Sd, a, H, ap, Hp, tturn, apd, Rpdd, Sdd, add, apdd
		! ---------------------------------------------------------- !
		! MEMORY:
		! ---------------------------------------------------------- !
! MOved to lltb_params, such that we can pass it to error report.
!		type LB_MYMEM
!			real(dl) :: H0_inf, Lambda, Mt, r, t, k
!			real(dl) :: t0, Mt2, Lambdaover3, H0_inf2
!			real(dl) :: Rltb, Rp, Rpd, S, Sd, a, H, ap, tturn, Hp, apd
!		end type LB_MYMEM
		type(LB_MYMEM), save :: lbmem
		! ---------------------------------------------------------- !
		! Local use:
		! ---------------------------------------------------------- !
		type(voidparams) :: vpi
		real(dl) :: thisk, thisdkdr, thisa, Mt2, H0_inf2, Lambdaover3
		real(dl) :: thisrltb, thisrp, thisrpd, thiss, thissd, thish, thishp
		real(dl) :: thisap, thisdtbbdr, thisdtdr, thisadot, thisdhda, thisdhpartdr
		real(dl) :: thistminustbb, thisdtdr_amax
                real(dl) :: thistturn
		logical :: model_cached, coords_cached
                real(dl) :: asigma, hsign,  thisadotterm, thisaprime_amax, thisadotterm_amax
                real(dl) :: thisadotdot, thisapdot, thishdot, thisapdot_overa, thishp_amax
                real(dl) :: thisapdot_overa_amax
                real(dl) :: thisrpdd, thisSdd, thisapdotdot
    real(dl) :: rescaleH, rescaleM, rst, rsr

!call lltb_feedback(3)

		model_cached=.true.
		coords_cached=.true.

                ! Need to do the dangerous check of kofr vs memory,
                ! because in PlotDensity we subsequently call this 
                ! routine for LTB and FLRW at exactly same coordinates,
                ! the only difference being VP%L => 0.
		thisk = kofr(r)


		! Cachable actions, anything r-independent:
		if(.not. all( (/lbmem%H0_inf, lbmem%Lambda, lbmem%mt/) == &
				 (/H0_inf, Lambda, Mtilde/) ) ) then

				lbmem%h0_inf = h0_inf
				lbmem%Lambda = Lambda
				lbmem%mt = Mtilde

				! do actions...
				! This one is a bit lame, but shows the idea:
				lbmem%mt2 = Mtilde**2
				lbmem%h0_inf2 = h0_inf**2
				lbmem%Lambdaover3 = Lambda/3._dl
                                
                                lbmem%t = t
                                lbmem%r = r
                                lbmem%k = thisk

				model_cached = .false.
				coords_cached = .false.
		else ! Else the model is the same, but maybe coordinates changed:

			if(.not. all((/r,t,thisk/)==(/lbmem%r,lbmem%t,lbmem%k/)))then
                           
                           lbmem%t = t
                           lbmem%r = r
                           lbmem%k = thisk

				coords_cached=.false.
			end if

		end if

    rescaleH = 1.d0!sqrt(1.d0/max(Lambda,lbmem%mt2))
    rescaleM = 1.d0!rescaleH


		if(.not.(model_Cached .and. coords_cached))then

                        call lltb_set_error(vpi,lbmem)

      rst = t / rescaleH
      rsr = r / rescaleM

			Mt2 = lbmem%mt2 * rescaleM**2
			H0_inf2 = lbmem%H0_inf2 * rescaleH**2
			Lambdaover3 = lbmem%Lambdaover3 * rescaleH**2

			! For error reporting purposes only:
			vpi%mtilde = lbmem%mt * rescaleM
			vpi%lambda = lbmem%lambda * rescaleH

			thisdkdr=dkdr(r) * rescaleM

			thistminustbb = (t - tbbofr(r)) / rescaleH
			if(thistminustbb.le.0.d0)then
                           write(*,*)'Input t:',t
                           call lltb_error('Beyond Big Bang in lltb_background (note that info below may be too old).')
                        end if

			! set omega's. Don't have to sum up to one.
			vpi%xyz(1) = 8._dl * pi * Mt2 / 3._dl
			vpi%xyz(2) = 2._dl * Mt2 * thisk
			vpi%xyz(3) = Lambdaover3


			! get a(r,t): always. This the weakest link for time consumption.
			call lltbbackgroundaofHt(vpi,thisa,thistminustbb,sigma=asigma)
                        if(isnan(thisa))call lltb_error('isnan(thisa)')
                        ! get tturn, which is set in aofHt:
                        thistturn=vpi%tturn
                        hsign=1._dl
if(lfeedback>3)then
  write(*,*)'After aofHt, ExpCase:'
  write(*,*)vpi%ExpCase
end if


			! get a'(r,t).
			thisdtbbdr = dtbbdr(r)
                        if(abs((thisa-vpi%smallestlinroot)/thisa) > accuracy(dl)*acceptance) then
                           call ddok_tofa(vpi,a=thisa, t=thisdtdr)
                           ! Now thisdtdr contains d t / d (Omega_k)
                           ! But, dt/dr = dt/dOm_k * dOm_k/dr   +                    
                           thisdtdr = thisdtdr * 2._dl * Mt2 * thisdkdr

!                        else ! When we approach the turning point too closely, take care of analytical limit (cancelling infinities):
!
!!!                           call ddok_tofamax(vpi, t=thisdtdr)
!                           thisaprime_amax = aprime_amax(vpi,a=vpi%smallestlinroot,yprime= 2._dl * Mt2 *thisdkdr,dtbbdr=0._dl, hsign=1._dl,adotterm=thisadotterm_amax)
!                           thisdtdr = -thisadotterm_amax
                        end if
                        ! I don't get it: the following solves a bug for small densityies at high z.
                        if(lfeedback<-10000)write(*,*)'thisdtdr:',thisdtdr
                        if(thistminustbb>vpi%tturn)then
                           hsign=-1._dl
                        end if

                        if((abs((thisa-vpi%smallestlinroot)/thisa) <= accuracy(dl)*acceptance) .or. &
                          (thistminustbb>vpi%tturn)   ) then
                           thisaprime_amax = aprime_amax(vpi,a=vpi%smallestlinroot,yprime= 2._dl * Mt2 *thisdkdr,dtbbdr=0._dl, hsign=1._dl,adotterm=thisadotterm_amax)

                           if((abs((thisa-vpi%smallestlinroot)/thisa) <= accuracy(dl)*acceptance))then
                              thisdtdr = -thisadotterm_amax
                           end if
                        end if

                        ! Are we beyond the turning point? Then an extra minus sign
                        ! (see ddok_tofa(..)) AND ddok_tmax(r)
                        if(thistminustbb>vpi%tturn)then
!                           thisaprime_amax = aprime_amax(vpi,a=vpi%smallestlinroot,yprime= 2._dl * Mt2 *thisdkdr,dtbbdr=0._dl, hsign=1._dl,adotterm=thisadotterm_amax)
                           !                           call ddok_tofamax(vpi, t=thisdtdr_amax)
                           !                           thisdtdr_amax = thisdtdr_amax * 2._dl * Mt2 * thisdkdr
                           
                           thisdtdr_amax = -thisadotterm_amax

                           thisdtdr = - thisdtdr  + 2._dl* thisdtdr_amax
                           if(isnan(thisdtdr))then
                              write(*,*)thisdtdr
                              writE(*,*)thisdtdr_amax
                              !call ddok_tofamax(vpi, t=thisdtdr_amax)
                              !writE(*,*)thisdtdr_amax
                              call lltb_error('isnan(dtdr_amax)')
                           end if
                        end if

			thish = hsign*HoverH0(vpi%xyz,thisa)
                        if(isnan(thish))call lltb_error('isnan(thish)') ! This should actually be impossible, secured in HoverH0 already.
			thisadot = thish * thisa

                        thishdot = thisa * 0.5_dl * dH2daoverH0(vpi%xyz,thisa)
                        !                        thisadotdot = -(3._dl * vpi%xyz(1) / thisa**3 + 2._dl * vpi%xyz(2)/thisa**2)
                       thisadotdot = adotdotoverH0(vpi%xyz,thisa)!-.5_dl * vpi%xyz(1) / thisa**2 + vpi%xyz(3)*thisa

                        if(abs((thisa-vpi%smallestlinroot)/thisa) <  accuracy(dl)*acceptance) then
                           !                           thisap = aprime_amax(vpi,a=thisa,yprime= 2._dl * Mt2 *thisdkdr,dtbbdr=thisdtBBdr, hsign=hsign,adotterm=thisadotterm)
                           thisap = aprime_amax(vpi,a=thisa,yprime= 2._dl * Mt2 *thisdkdr,dtbbdr=thisdtBBdr, hsign=1._dl,adotterm=thisadotterm)
                           
                           
                           thisapdot = 0.5_dl / thisadot *((-vpi%xyz(1)/thisa**2+2._dl*vpi%xyz(3)*thisa)*thisap + 2._dl* Mt2 *thisdkdr)
                           
                           
                           
                           thishp = Hprime_amax(vpi,a=thisa,yprime= 2_dl * Mt2 *thisdkdr, hsign=hsign,adotterm=thisadotterm,aprime=thisap)
                           !                           thishp = Hprime_amax(vpi,a=thisa,yprime= 2_dl * Mt2 *thisdkdr, hsign=1._dl,adotterm=thisadotterm,aprime=thisap)

                           if(thistminustbb>vpi%tturn)then

                              !                              thishp_amax = Hprime_amax(vpi,a=vpi%smallestlinroot,yprime= 2_dl * Mt2 *thisdkdr, hsign=1._dl,adotterm=thisadotterm,aprime=thisaprime_amax)
                              !                              thisapdot_overa_amax = (thishp_amax)
                              thisapdot_overa = (thishp + thisap / thisa * thish*hsign) ! everything for expansion 
                              thisap = 2._dl * thisaprime_amax - thisap ! this line works, although it has a small error.
                              thisapdot_overa = thisapdot_overa
  thishp = Hprime_amax(vpi,a=thisa,yprime= 2_dl * Mt2 *thisdkdr, hsign=hsign,adotterm=2._dl*thisadotterm_amax-thisadotterm,aprime=thisap)
                              
!                              thishp = thisapdot_overa - thisap / thisa * thish
                           end if
                           
                        else
                           thisap = - thisadot  * ( thisdtdr + thisdtbbdr)
                           
                           ! Using the correct aprime, also aprimedot and Hprime are continuous in 
                           ! a=u (adot=0). Checked a million times, turned out to have a bug that
                           ! caused a difference O(10^(-7) in 
                           ! ddok_tofamax, which led to a finite difference 
                           ! in Hp \propto ddok_tofamax/adot.. Corrected now.

                           !----!
                           ! Next apdot is correct, checked:
                           thisapdot = 0.5_dl / thisadot *((-vpi%xyz(1)/thisa**2+2._dl*vpi%xyz(3)*thisa)*thisap + 2._dl* Mt2 *thisdkdr)
                           !----!
                           
                           ! get H'(r,t)
                           
                           !----!
                           !			thishp = (0.5_dl* thisap * thisadotdot /thisadot + Mt2 * thisdkdr / thisadot / thisa)
                           ! The following is equal to the previous:
                           
                           thishp = thisapdot / thisa - thisap / thisa * thish
                           !----!
                           
                           
                        end if
                        

                        if(isnan(thisap))then
                           write(*,*) thisadot
                           write(*,*) thisdtdr
                           write(*,*) thisdtbbdr
                           call lltb_error('isnan(thisap)')
                        end if

                        if(isnan(thishp))then
                           write(*,*)'Thishp is NaN.'
!                           write(*,*)'thisdhda:',thisdhda
                           write(*,*)'thisap:',thisap
!                           write(*,*)'thisdhpartdr:',thisdhpartdr
                           write(*,*)'thisa:',thisa
                           write(*,*)'thishp:',thishp
                           write(*,*)'thisa:',thisa
                           !   call lltb_error('isnan(thishp)')
                        end if


			! Get all functions:
			thisrltb = rsr*thisa
                        if(isnan(thisrltb))call lltb_error('isnan(thisrltb)')
			thisrp = thisa + rsr*thisap
                        if(isnan(thisrp))call lltb_error('isnan(thisrp)')
                        ! thisrd = r*thisadot ! Not used.

!			thisrpd = thisadot + r*thish*thisap + r*thishp*thisa

!			thisrpd = thisadot - r*thish*thisap + r*thishp*thisa
			thisrpd =  thisadot + rsr*thish*thisap + rsr*thishp*thisa

                        if(isnan(thisrpd))call lltb_error('isnan(thisrpd)')
!			thisrpd = thisadot + r*thish*thisap - r*thishp*thisa

!                        thisrpd = thisrltb* thishp + thisrp*thish

			thiss = thisrp / (1._dl + 2._dl * rsr**2 * Mt2 * thisk)**(0.5_dl)
			thissd = thisrpd / (1._dl + 2._dl * rsr**2 * Mt2 * thisk)**(0.5_dl)

                        if(1._dl + 2._dl * rsr**2 * Mt2 * thisk .le. 0.d0)call lltb_error('Exceeding maximum radius of closed universe: this coordinate r does not exist.')
                        if(isnan(thiss))call lltb_error('isnan(thiss)')
                        if(isnan(thissd))call lltb_error('isnan(thissd)')

                        ! New: ddot's
                        thisapdotdot = (vpi%xyz(1) / thisa**3 + vpi%xyz(3))*thisap
                        thisRpdd = thisadotdot + rsr * thisapdotdot
                        thisSdd = thisrpdd / (1._dl + 2._dl * rsr**2 * Mt2 * thisk)**(0.5_dl)
if(lfeedback>4)then
  if(thiss<0)then
    call lltb_error('g_rr = S < 0')
    stop
  end if
end if

			lbmem%Rltb=thisrltb   
			lbmem%Rp=thisrp
			lbmem%Rpd=thisrpd
			lbmem%Rpdd=thisrpdd
			lbmem%S=thiss
			lbmem%Sd=thissd
			lbmem%Sdd=thissdd
			lbmem%a=thisa
			lbmem%ap=thisap
			lbmem%apd=thisapdot
			lbmem%add=thisadotdot
			lbmem%apdd=thisapdotdot
			lbmem%H=thish
			lbmem%Hp=thishp
			lbmem%tturn=thistturn
!else
!write(*,*)'Caching!' ! Tested: it works.
		end if


		! done:
		if(present(Rltb))Rltb=lbmem%rltb * rescaleM
		if(present(Rp))Rp=lbmem%rp
		if(present(Rpd))Rpd=lbmem%rpd / rescaleH
		if(present(Rpdd))Rpdd=lbmem%rpdd / rescaleH**2
		if(present(S))S=lbmem%s
		if(present(Sd))Sd=lbmem%sd / rescaleH
		if(present(Sdd))Sdd=lbmem%sdd / rescaleH**2
		if(present(a))a=lbmem%a
		if(present(ap))ap=lbmem%ap / rescaleM
		if(present(apd))apd=lbmem%apd / rescaleH / rescaleM
		if(present(add))add=lbmem%add / rescaleH**2
		if(present(apdd))apdd=lbmem%apdd / rescaleH**2 / rescaleM
		if(present(Hp))Hp=lbmem%hp / rescaleH / rescaleM
		if(present(H))H=lbmem%h / rescaleH
		if(present(tturn))tturn=lbmem%tturn * rescaleH

	end subroutine lltb_functions_all

        ! Same, using a as time parameter, t as output
	subroutine lltb_functions_a(H0_inf, Lambda, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, a, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, t, ap, apd, add, apdd, H, Hp, tturn)
		implicit none
		! ---------------------------------------------------------- !
		! INPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(in) :: H0_inf, Lambda
		interface
			function kofr(rr)
				use lltb_precision
				real(dl) :: kofr
				real(dl), intent(in) :: rr
			end function kofr
		end interface
		interface
			function dkdr(rr)
				use lltb_precision
				real(dl) :: dkdr
				real(dl), intent(in) :: rr
			end function dkdr
		end interface
		interface
			function tbbofr(rr)
				use lltb_precision
				real(dl) :: tbbofr
				real(dl), intent(in) :: rr
			end function tbbofr
		end interface
		interface
			function dtbbdr(rr)
				use lltb_precision
				real(dl) :: dtbbdr
				real(dl), intent(in) :: rr
			end function dtbbdr
		end interface
		real(dl), intent(in) :: r, a, Mtilde
		! ---------------------------------------------------------- !
		! OUTPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(out), optional :: Rltb, Rp, Rpd, S, Sd, t, H, ap, Hp, tturn, apd, Rpdd, Sdd, add, apdd
		! ---------------------------------------------------------- !
		! MEMORY:
		! ---------------------------------------------------------- !
! MOved to lltb_params, such that we can pass it to error report.
!		type LB_MYMEM
!			real(dl) :: H0_inf, Lambda, Mt, r, t, k
!			real(dl) :: t0, Mt2, Lambdaover3, H0_inf2
!			real(dl) :: Rltb, Rp, Rpd, S, Sd, a, H, ap, tturn, Hp, apd
!		end type LB_MYMEM
		type(LB_MYMEM), save :: lbmem
		! ---------------------------------------------------------- !
		! Local use:
		! ---------------------------------------------------------- !
		type(voidparams) :: vpi
		real(dl) :: thisk, thisdkdr, thisa, Mt2, H0_inf2, Lambdaover3
		real(dl) :: thisrltb, thisrp, thisrpd, thiss, thissd, thish, thishp
		real(dl) :: thisap, thisdtbbdr, thisdtdr, thisadot, thisdhda, thisdhpartdr
		real(dl) :: thistminustbb, thisdtdr_amax
                real(dl) :: thistturn
		logical :: model_cached, coords_cached
                real(dl) :: asigma, hsign,  thisadotterm, thisaprime_amax, thisadotterm_amax
                real(dl) :: thisadotdot, thisapdot, thishdot, thisapdot_overa, thishp_amax
                real(dl) :: thisapdot_overa_amax
                real(dl) :: thisrpdd, thisSdd, thisapdotdot

                real(dl) :: myt,mya

!call lltb_feedback(3)

                

                Mt2 = Mtilde**2
                
                Lambdaover3 = Lambda/3._dl
                
                vpi%xyz(1) = 8._dl * pi * Mt2 / 3._dl
                vpi%xyz(2) = 2._dl * Mt2 * kofr(r)
                vpi%xyz(3) = Lambdaover3
                
                call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz))
 
                ! SET t
                call lltbbackgroundHtofa(vpi,a=a,t=myt)

                if(present(t))t=myt
                
                ! Passing optional arguments to optional arguments, without explcicitly checking presence.
                ! Not sure this works with all compilers. Or is it Fortran standard? Works with ifort.
                call lltb_functions(H0_inf, Lambda, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, myt, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, mya, ap, apd, add, apdd, H, Hp, tturn)

                !write(*,*)'testing lltb_functions_a. a_in / a_out -1:',a/mya-1.d0

        end subroutine lltb_functions_a

        ! Wrapper: just necessary to force abs(r)
	subroutine lltb_functions(H0_inf, Lambda, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp, tturn)
		implicit none
		! ---------------------------------------------------------- !
		! INPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(in) :: H0_inf, Lambda
		interface
			function kofr(rr)
				use lltb_precision
				real(dl) :: kofr
				real(dl), intent(in) :: rr
			end function kofr
		end interface
		interface
			function dkdr(rr)
				use lltb_precision
				real(dl) :: dkdr
				real(dl), intent(in) :: rr
			end function dkdr
		end interface
		interface
			function tbbofr(rr)
				use lltb_precision
				real(dl) :: tbbofr
				real(dl), intent(in) :: rr
			end function tbbofr
		end interface
		interface
			function dtbbdr(rr)
				use lltb_precision
				real(dl) :: dtbbdr
				real(dl), intent(in) :: rr
			end function dtbbdr
		end interface
		real(dl), intent(in) :: r, t, Mtilde
		! ---------------------------------------------------------- !
		! OUTPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(out), optional :: Rltb, Rp, Rpd, S, Sd, a, H, ap, Hp, tturn, apd, Rpdd, Sdd, add, apdd


                call lltb_functions_all(H0_inf, Lambda, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						abs(r), t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp, tturn)


              end subroutine lltb_functions


        ! lltb_normalization returns Mtilde, Hlocal, delta0 and t0 
        ! as a function of H0(out),kb,kmax and Lambda.
        ! Normalization is such that a(L,t0)=1.
        !
        ! set kb= kb+kmax and kmax = 0 to normalize to a(0,t0)=1
        !
        subroutine lltb_normalization(H0,kb,kmax,Lambda, & 
             Mtilde,Hlocal,delta0,t0)
          implicit none
          real(dl), intent(in) :: H0, kb, kmax, Lambda ! input
          real(dl), intent(out) :: Mtilde, Hlocal, delta0, t0 ! output
          
          type(voidparams) :: vpi
          real(dl) :: Mt2
          real(dl) :: Lambdaover3
          real(dl) :: thisa
          real(dl) :: rho_out
          real(dl) :: rho_in
          real(dl) :: turnovera,turnovertime
          character(len=1) :: thissubcase(3)
          integer :: n_fullroots
          real(dl) :: aturn_fullroots

          Mtilde = sqrt(3._dl * (H0**2 - Lambda/3._dl)/8._dl/pi/(1._dl+0.75_dl * kb/pi))
          Mt2 = Mtilde**2

          Lambdaover3 = Lambda/3._dl

          vpi%xyz(1) = 8._dl * pi * Mt2 / 3._dl
          vpi%xyz(2) = 2._dl * Mt2 * kb
          vpi%xyz(3) = Lambdaover3

          call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz))
          if(vpi%smallestrealroot<1.d0)then ! then even the outside universe is collapsing. We don't want that.
             if(lfeedback>1)write(*,*)'Normalization not working: collapsing universe.'
             ! set all output to -1:
              Mtilde=-1.
              Hlocal=-1.
              delta0=-1.
              t0=-1.
              return
          end if
          !-------!
          ! SET T0
          call lltbbackgroundHtofa(vpi,a=1._dl,t=t0)

          ! Both at r=0 and r=inf, R'= a, a'=0
          rho_out = 1._dl ! times Mtilde^2 times Mp^2

          !-------!
          ! set a0_in
          ! vpi%xyz(1) = 8._dl * pi * Mt2 / 3._dl
          vpi%xyz(2) = 2._dl * Mt2 * ( kb + kmax )
          ! vpi%xyz(3) = Lambdaover3

          ! It is very very well possible that for this t0 and this kb+kmax
          ! the metric at this radius r=0 has already recollapsed to its big crunch.
          ! So before we ask for a(t), which may give an error (t>t_collapse)
          ! we must check what is the crunch time:

          thissubcase=GetOmegaCase(vpi%xyz)
          call lltb_getallroots(vpi,thissubcase)
          turnovera=min(vpi%smallestrealroot,vpi%smallestlinroot)
          aturn_fullroots = turnovera
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
                   call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz,turnovera))
                   turnovera = min(vpi%smallestrealroot,vpi%smallestlinroot)
                   ! End trouble.
                   ! HOWEVER, this solution should not hide the possible full roots:
                   if(n_fullroots>0 .and. vpi%posrealroots==0)then
                     vpi%posrealroots=n_fullroots
                     turnovera=aturn_fullroots
                   end if



          if(turnovera<1.e20_dl)then
             call lltbbackgroundHtofa(vpi,a=turnovera,t=turnovertime)
          else
             turnovertime=infinity !1.e30_dl
          end if
          
!          if(turnovertime>t0)then

          if(2._dl*turnovertime>t0)then

             call lltbbackgroundaofHt(vpi,thisa,t0)

          
             !-------!
             ! set H0_in
             Hlocal = HoverH0(vpi%xyz,thisa)
             if(t0>vpi%tturn)Hlocal=-Hlocal
             !-------!
             ! set d0
             rho_in = 1._dl / thisa**3 ! times Mtilde^2 times Mp^2
             ! a*3, because we demand k' = a' = 0 at r=0
             delta0 = rho_in / rho_out - 1._dl
          else


             !-------!
             ! set H0_in
             Hlocal = -infinity !-1.e30_dl

             
             !-------!
             ! set d0
             delta0 = infinity !1.e30_dl
          end if



        end subroutine lltb_normalization

        ! lltb_feedback: set verbosity of lltb module (0=default, higher is more)
        subroutine lltb_feedback(number)
          integer, intent(in) :: number
          
          lfeedback=number
        end subroutine lltb_feedback

        ! lltb_precisiongoal: can be called with real*4 and real*8
        ! sets the precisiongoal: default = 1.e-16, cannot set it smaller, but larger.
        ! Note that this is a goal, not an absolute limit.
        ! accuracy will normally vary between 1.e-11 and 1.e-16
        ! variation will rescale with the rescaled precisiongoal
        subroutine lltb_precisiongoal_sl(goal)
          real, intent(in) :: goal

          call lltb_precisiongoal_dl(dble(goal))

        end subroutine lltb_precisiongoal_sl

        subroutine lltb_precisiongoal_dl(goal)
          real(dl), intent(in) :: goal
          real(dl) :: orig_goal
          real(dl) :: power
          ! We initialize with precision goal = 1.d-16
          
          orig_goal = accuracy_p(dl) ! parameter accuracy
          
          if(goal < orig_goal)then
             write(*,'(A,ES14.5,A,ES14.5,A)')'You demand precisiongoal ',goal,' which is smaller than the hardcoded lowerlimit of ',orig_goal,'. That will not work.'
             stop 'Too high precisiongoal demanded. LLTBBackground.'
          end if

          power = log(goal)/log(orig_goal)

          ! Here we go. Be careful, this list must be complete. Compare to module lltb_precision.
          accuracy=accuracy_p**power
          lin_expansion=lin_expansion_p**power
          acceptance_pure=acceptance_pure_p**power
          acceptance=acceptance_p**power
          acceptance_pole=acceptance_pole**power
          accuracy_downscale=accuracy_downscale_p*power ! Note times, not to the power
          rel_acceptance = rel_acceptance_p**power

          call Elliptics_precisiongoal(goal)

        end subroutine lltb_precisiongoal_dl

        subroutine lltb_attach_to_error(msg)
          implicit none
          character(len=*), intent(in) :: msg

          call lltb_set_error(attachement=msg)

        end subroutine lltb_attach_to_error
        
        
#include "lltb_ricci.f90"
        
end module lltb
