module wlltb !
  use wlltb_constants
  use wlltb_types_funcs
  use wlltb_funcs
  use lltb, only: lltb_ricci
  use wlltb_errors, only: wlltb_error, wlltb_set_error  
  private

  public lltb_functions, lltb_functions_a, lltb_normalization, lltb_feedback, lltb_precisiongoal, lltb_attach_to_error
  public lltb_ricci, wlltb_error

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
	subroutine lltb_functions_all(H0_inf, Lambda, w, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp, tturn, For3,deltarhow,deltam)
		use lltb, only: oldlltb => lltb_functions
    use wlltb_funcs
    use wlltb_types
    use wlltb_integratormod
    implicit none
    
		! ---------------------------------------------------------- !
		! INPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(in) :: H0_inf, Lambda
                interface
			function myfunc(rr)
				use wlltb_constants
				real(dl) :: myfunc
				real(dl), intent(in) :: rr
			end function myfunc
		end interface
		procedure(myfunc)  :: kofr, dkdr, tbbofr, dtbbdr
		real(dl), intent(in) :: r, t, Mtilde
		real(dl), intent(in), optional :: w
		! ---------------------------------------------------------- !
		! OUTPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(out), optional :: Rltb, Rp, Rpd, S, Sd, a, H, ap, Hp, tturn, apd, Rpdd, Sdd, add, apdd, For3,deltarhow,deltam
        ! Same, using a as time parameter, t as output
		! ---------------------------------------------------------- !
		! LOCAL:
		! ---------------------------------------------------------- !
    real(dl) :: tbar, w_eos, t0, abar, Hbarout, rhoMbar, rhoMoutbar, y(index%max),dy(index%max), Eofr, kb
    real(dl) :: FMbaroverout, mrhowout, mrhomout, mrhomin, thist, thisa
!    real(dl) :: FMbaroutoverr3, Fwbaroutoverr3, FMbaroverout, Fwbaroverout, Hbarout
!    real(dl) :: FMbar, Fwbar
    type(funcscontainer) :: radfuncs
    type(LB_MYMEM) :: lbnew
    type(LB_MYMEM), save :: lbmem
    type(wlltbparams) :: vpi
    
    if(.not.present(w))then
      w_eos = -1.d0
    else
      w_eos = w
    end if

!if(present(deltam))lfeedback=4

    lbnew%H0_Inf = H0_Inf
    lbnew%Lambda = Lambda
    lbnew%w = w_eos
    lbnew%k = kofr(r)
    lbnew%dk = dkdr(r)
    lbnew%tbb = tbbofr(r)
    lbnew%dtbb = dtbbdr(r)
    lbnew%Mt = Mtilde
    lbnew%t = t
    lbnew%r = r
    
    vpi%dummy = 0.d0
    
    call wlltb_set_error(vpi,lbnew)

!if(.false.)then
    if(abs(w_eos + 1.d0)<1.d-30 .and. .not. present(deltam))then
      if(lfeedback>0)write(*,*)'Calling oldLTB'
      call oldlltb(H0_inf, Lambda, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, thisa, ap, apd, add, apdd, H, Hp, tturn)
      if(present(deltarhow))deltarhow = 0.d0
      if(present(For3))For3=4.d0*Mtilde**2*pi/3.d0 + Lambda / 6.d0 * thisa**3 ! 4 pi and 1/6, vs 8 pi and 1/3.
      if(present(a))a=thisa
      if(present(tturn).and.lfeedback>0)write(*,*)'tturn in wlltb:',tturn
      return
    end if
!end if

    radfuncs%kofr => kofr
    radfuncs%dkdr => dkdr
    radfuncs%tbbofr => tbbofr
    radfuncs%dtbbdr => dtbbdr
    
    


    ! Caching:
    if ( .not. lbnew == lbmem ) then ! see wlltb_types_funcs for comparison operator.
    
      if(abs(dkdr(0.d0))>1.d-30)write(0,*)'dkdr(0) /= 0: that is not consistent. The code will give wrong results for r-derivative functions.'
      
!      if(lfeedback>3)then
      if(.false.)then
        write(*,*)lbnew%H0_Inf,lbmem%H0_Inf
        write(*,*)lbnew%Lambda,lbmem%Lambda
        write(*,*)lbnew%w,lbmem%w
        write(*,*)lbnew%k,lbmem%k
        write(*,*)lbnew%dk,lbmem%dk
        write(*,*)lbnew%tbb,lbmem%tbb
        write(*,*)lbnew%dtbb,lbmem%dtbb
        write(*,*)lbnew%Mt,lbmem%Mt
        write(*,*)lbnew%t,lbmem%t
        write(*,*)lbnew%r,lbmem%r
      end if

      lbmem = lbnew
    
      call wlltb_initime(H0_inf, Lambda, w, radfuncs, Mtilde, r, &
            t0, tbar, abar, Hbarout, rhoMbar, rhoMoutbar, kb, FMbaroverout, y)

      lbmem%t0 = t0
      if(lfeedback>1)write(*,*)'wlltb tbar,t0,t:',tbar,t0,t
      if(lfeedback>1)write(*,*)'wlltb k(r),k''(r):',lbmem%k,lbmem%dk,lbmem==lbnew
      call wlltb_set_error(vpi,lbmem)

!if(t>t0*(1.d0+1.d-15))then
!  call wlltb_error('t>t0')
!else if (t>t0) then
!  thist = t0
!else
!  thist = t
!end if
! t can be > t0, when we try to find FLRW observer. Must be healthy though..
if(t>1.d2*t0)then
  call wlltb_error('t>t0')
end if

      if(any(abs(y)>1.d29))then
        write(0,*)'y:',y
        write(0,*)'r:',r
        write(0,*)H0_inf, Lambda, w, radfuncs%kofr(r), radfuncs%dkdr(r), Mtilde, r
        call wlltb_error('kanker')
        stop 'kanker'
      end if
      lbmem%t0 = t0
      
      call wlltb_derivs(index%max,tbar,y,dy)
      if(t>tbar)then
        if(lfeedback>1)write(*,*)'wlltb yin:',y
        call wlltb_integrator(H0_inf, Lambda, w, radfuncs, Mtilde, kb, &
						r, t, tbar, rhoMbar, y)
!						r, thist, tbar, rhoMbar, y)
        if(lfeedback>0)write(*,*)'wlltb yout:',y
      else
        if(isnan(t).or.(t<0.d0))then
          y=0.d0
          y(index%altb)=0.d0
          y(index%adot)=1.d31
          if(lfeedback>1)write(*,*)'wlltb t is sick:',t
        else
          if(lfeedback>1)write(*,*)'wlltb t<tbar:',t,'<',tbar
        end if
        ! else: continue with y(t=tbar)
      end if

      call wlltb_derivs(index%max,t,y,dy)

      lbmem%Rltb = r * y(index%altb)
      lbmem%a = y(index%altb)
!     lbmem%Rp = r * y(index%aprime) + y(index%altb)
      lbmem%Rp = y(index%raprime) + y(index%altb)
      lbmem%Rpd = dy(index%raprime) + dy(index%altb)
		  lbmem%Rpdd = radotdotprime(y) + dy(index%adot)
      lbmem%add = dy(index%adot)

      lbmem%F = r**3 * y(index%For3)
      lbmem%For3 = y(index%For3)
      lbmem%rhow = y(index%rhow)
      mrhowout = rhowout(y)
      
      lbmem%deltarhow = y(index%rhow)/mrhowout - 1.d0
      if(abs(w_eos + 1.d0)<1.d-30)lbmem%deltarhow=0.d0
      
      mrhomout = rhoMoutbar * (norms%abar/y(index%aout))**3
      mrhomin = rhoMbar * norms%abar**3 / lbmem%a**2 / lbmem%Rp
      lbmem%deltam = mrhomin/mrhomout - 1.d0
      ! if(present(deltam)) write(*,*)'deltam',mrhomout,mrhomin,rhoMbar,rhoMoutbar,y(index%aout),norms%a0, lbmem%a , lbmem%Rp, lbmem%deltam,y(index%raprime),y(index%altb)
      ! primes:
      if(abs(r)>0.d0)then
        lbmem%apd = dy(index%raprime)/r
        lbmem%ap = y(index%raprime)/r
        lbmem%apdd = radotdotprime(y) / r
        lbmem%hp = ( dy(index%raprime)/y(index%aLtb) - dy(index%altb)*y(index%raprime)/y(index%aLtb)**2  ) / r
      else
        lbmem%apd = 0.
        lbmem%ap = 0.
        lbmem%apdd = 0.
        lbmem%hp = 0.
      end if
        
      
      lbmem%H = dy(index%altb)/y(index%altb)

      Eofr = r**2 * kofr(r) * Mtilde**2

      lbmem%s = lbmem%Rp / (1.d0 + 2.d0 * Eofr)**0.5d0
      lbmem%sd = lbmem%Rpd / (1.d0 + 2.d0 * Eofr)**0.5d0
      lbmem%sdd = lbmem%Rpdd / (1.d0 + 2.d0 * Eofr)**0.5d0
      
      ! Last: tturn. Need to re-use integration variables, resetting everything.
      if(present(tturn) .and. kofr(r)<0.d0)then
      call wlltb_initime(H0_inf, Lambda, w, radfuncs, Mtilde, r, &
            t0, tbar, abar, Hbarout, rhoMbar, rhoMoutbar, kb, FMbaroverout, y)
      lbmem%t0 = t0
      
      call wlltb_derivs(index%max,tbar,y,dy)
        call wlltb_integratortturnonly(H0_inf, Lambda, w, radfuncs, Mtilde, kb, &
						r, t0, tbar, rhoMbar, y,lbmem%tturn)
      else
  		  lbmem%tturn = 1.d30
      end if

    else
!    write(*,*)'Caching!' ! Tested: it works.
      if(lfeedback>0)write(*,*)'wlltb same model: cached results.'
      continue    
    end if
    

		! done:
		if(present(Rltb))Rltb=lbmem%rltb
		if(present(Rp))Rp=lbmem%rp
		if(present(Rpd))Rpd=lbmem%rpd
		if(present(Rpdd))Rpdd=lbmem%rpdd
		if(present(S))S=lbmem%s
		if(present(Sd))Sd=lbmem%sd
		if(present(Sdd))Sdd=lbmem%sdd
		if(present(a))a=lbmem%a
		if(present(ap))ap=lbmem%ap
		if(present(apd))apd=lbmem%apd
		if(present(add))add=lbmem%add
		if(present(apdd))apdd=lbmem%apdd
		if(present(Hp))Hp=lbmem%hp
		if(present(H))H=lbmem%h
		if(present(tturn))tturn=lbmem%tturn
		if(present(For3))For3=lbmem%For3
		if(present(deltarhow))deltarhow=lbmem%deltarhow
		if(present(deltam))deltam=lbmem%deltam
    
    if(lfeedback>3)write(*,*)'returning from wlltb'

  

        
  end subroutine lltb_functions_all
  
  
	subroutine lltb_functions_a(H0_inf, Lambda, w, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, a, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, t, ap, apd, add, apdd, H, Hp, tturn, For3,deltarhow,deltam)
		implicit none
		! ---------------------------------------------------------- !
		! INPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(in) :: H0_inf, Lambda
		interface
			function kofr(rr)
				use wlltb_constants
				real(dl) :: kofr
				real(dl), intent(in) :: rr
			end function kofr
		end interface
		interface
			function dkdr(rr)
				use wlltb_constants
				real(dl) :: dkdr
				real(dl), intent(in) :: rr
			end function dkdr
		end interface
		interface
			function tbbofr(rr)
				use wlltb_constants
				real(dl) :: tbbofr
				real(dl), intent(in) :: rr
			end function tbbofr
		end interface
		interface
			function dtbbdr(rr)
				use wlltb_constants
				real(dl) :: dtbbdr
				real(dl), intent(in) :: rr
			end function dtbbdr
		end interface
		real(dl), intent(in) :: r, a, Mtilde
		real(dl), intent(in), optional :: w
		! ---------------------------------------------------------- !
		! OUTPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(out), optional :: Rltb, Rp, Rpd, S, Sd, t, H, ap, Hp, tturn, apd, Rpdd, Sdd, add, apdd, For3,deltarhow,deltam
		! ---------------------------------------------------------- !
		! MEMORY:
		! ---------------------------------------------------------- !
! MOved to lltb_params, such that we can pass it to error report.
!		type LB_MYMEM
!			real(dl) :: H0_inf, Lambda, Mt, r, t, k
!			real(dl) :: t0, Mt2, Lambdaover3, H0_inf2
!			real(dl) :: Rltb, Rp, Rpd, S, Sd, a, H, ap, tturn, Hp, apd
!		end type LB_MYMEM
!		type(LB_MYMEM), save :: lbmem
		! ---------------------------------------------------------- !
		! Local use:
		! ---------------------------------------------------------- !
!		type(voidparams) :: vpi
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

    

!    Mt2 = Mtilde**2
    
!    Lambdaover3 = Lambda/3._dl
    
!    vpi%xyz(1) = 8._dl * pi * Mt2 / 3._dl
!    vpi%xyz(2) = 2._dl * Mt2 * kofr(r)
!    vpi%xyz(3) = Lambdaover3
    
!    call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz))

    ! SET t
!    call lltbbackgroundHtofa(vpi,a=a,t=myt)

!    if(present(t))t=myt
    
    ! Passing optional arguments to optional arguments, without explcicitly checking presence.
    ! Not sure this works with all compilers. Or is it Fortran standard? Works with ifort.
    call lltb_functions(H0_inf, Lambda, w, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, myt, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, mya, ap, apd, add, apdd, H, Hp, tturn, For3,deltarhow,deltam)

                !write(*,*)'testing lltb_functions_a. a_in / a_out -1:',a/mya-1.d0

  end subroutine lltb_functions_a

        ! Wrapper: just necessary to force abs(r)
	subroutine lltb_functions(H0_inf, Lambda, w, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp, tturn, For3,deltarhow,deltam)
		implicit none
		! ---------------------------------------------------------- !
		! INPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(in) :: H0_inf, Lambda
		interface
			function kofr(rr)
				use wlltb_constants
				real(dl) :: kofr
				real(dl), intent(in) :: rr
			end function kofr
		end interface
		interface
			function dkdr(rr)
				use wlltb_constants
				real(dl) :: dkdr
				real(dl), intent(in) :: rr
			end function dkdr
		end interface
		interface
			function tbbofr(rr)
				use wlltb_constants
				real(dl) :: tbbofr
				real(dl), intent(in) :: rr
			end function tbbofr
		end interface
		interface
			function dtbbdr(rr)
				use wlltb_constants
				real(dl) :: dtbbdr
				real(dl), intent(in) :: rr
			end function dtbbdr
		end interface
		real(dl), intent(in) :: r, t, Mtilde
		real(dl), intent(in), optional :: w
		! ---------------------------------------------------------- !
		! OUTPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(out), optional :: Rltb, Rp, Rpd, S, Sd, a, H, ap, Hp, tturn, apd, Rpdd, Sdd, add, apdd, For3,deltarhow,deltam


    call lltb_functions_all(H0_inf, Lambda, w, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						abs(r), t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp, tturn, For3,deltarhow,deltam)


  end subroutine lltb_functions


  ! lltb_normalization returns Mtilde, Hlocal, delta0 and t0 
  ! as a function of H0(out),kb,kmax and Lambda.
  ! Normalization is such that a(L,t0)=1.
  !
  ! set kb= kb+kmax and kmax = 0 to normalize to a(0,t0)=1
  !
  subroutine lltb_normalization(H0,kb,kmax,Lambda,w, & 
       Mtilde,Hlocal,delta0,t0,tbar)
    use wlltb_types
    use wlltb_integratormod
		use lltb, only: oldlltb => lltb_normalization
    implicit none
    real(dl), intent(in) :: H0, kb, kmax, Lambda, w ! input
!		real(dl), intent(in), optional :: w          ! input
    real(dl), intent(out) :: Mtilde, Hlocal, delta0, t0, tbar ! output
!    real(dl), intent(out), optional :: tbar             ! output
!      real(dl) :: tbar
!    type(voidparams) :: vpi
    real(dl) :: Mt2
    real(dl) :: Lambdaover3
    real(dl) :: thisa
    real(dl) :: rho_out
    real(dl) :: rho_in, kbb, Foutor3, delta0ini
    real(dl) :: turnovera,turnovertime, FMbaroverout
    character(len=1) :: thissubcase(3)
    real(dl) ::  abar, Hbarout, rhoMbar, rhoMoutbar, y(index%max), dy(index%max), r, abar2Yprimebar, a2Yprime
    type(funcscontainer) :: radfuncs
    logical, save :: notwarned = .true.

    Mtilde = sqrt(3._dl * (H0**2 - Lambda/3._dl)/8._dl/pi/(1._dl+0.75_dl * kb/pi))

    Mt2 = Mtilde**2

        Hlocal=-1.


!      if(present(w))then

        radfuncs%kofr => dumfunc2
        radfuncs%dkdr => dumfunc3
        radfuncs%tbbofr => dumfunc
        radfuncs%dtbbdr => dumfunc

        r=0.d0
        !write(0,*)H0, Lambda, w, radfuncs%kofr(r), radfuncs%dkdr(r), Mtilde, r, kb
        call wlltb_initime(H0, Lambda, w, radfuncs, Mtilde, r, &
          t0, tbar, abar, Hbarout, rhoMbar, rhoMoutbar, kbb, FMbaroverout, y)
        thisa=t0


        if(abs(w + 1.d0)<1.d-30 )then
            call oldlltb(H0,kb,kmax,Lambda, & 
             Mtilde,Hlocal,delta0,t0)
             rhoMoutbar = 4*pi*Mtilde**2/3.d0
             rhoMbar = (delta0+1.d0)*rhoMoutbar
             delta0 = (rhoMbar + Lambda/6.d0)/(rhoMoutbar + Lambda/6.d0)-1.d0
            return
        end if
        

!        abar2Yprimebar = y(index%altb)**2 * (y(index%altb) + y(index%raprime))
        if(.not.any(y>1.d29) .and. allposdefsarepos(y) .and. t0>tbar .and. .not. isnan(tbar))then
          delta0ini = FMbaroverout - 1.d0


          call wlltb_integrator(H0, Lambda, w, radfuncs, Mtilde, kb, &
              0.d0, t0, tbar, rhoMbar, y)

          if((y(index%altb)>0.d0).and.(y(index%adot)<1.d30).and..not.any(isnan(y)))then
            call wlltb_derivs(index%max,t0,y,dy)            

            Hlocal = dy(index%altb)/y(index%altb)
          else

            Hlocal = 1.d30
          end if
  !        a2Yprime = y(index%altb)**2 * (y(index%altb) + y(index%raprime))

          ! outside: 

          thisa = (y(index%aout)/norms%a0) ! == 1

          Foutor3 = 0.5d0 * H0**2 * thisa**3 - Mt2 * kb * thisa 

          delta0 = y(index%For3) / Foutor3 * (y(index%aout) / y(index%altb))**3 - 1.d0
        else
          Mtilde = 1.d0
          Hlocal = 1.d30
          if(kmax>0.d0)then
            delta0 = -1.d0 
          else
            delta0 = 1.d30
          end if
          t0 = 0.d0
          tbar = 0.d0
        end if
        !writE(*,*)thisa, norms%a0, delta0, y(index%For3), Foutor3, kmax
!      else
!        stop 'write code for w=-1 case.. lltb_normalization.'
!      end if

!!      if(notwarned)then  
!        write(0,*)'Normalization is only a stub at this point.'
!        notwarned = .false.
!      end if
      
    return


!    call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz))
!    if(vpi%smallestrealroot<1.d0)then ! then even the outside universe is collapsing. We don't want that.
!       if(lfeedback>1)write(*,*)'Normalization not working: collapsing universe.'
!       ! set all output to -1:
!        Mtilde=-1.
!        Hlocal=-1.
!        delta0=-1.
!        t0=-1.
!        return
!    end if
    !-------!
    ! SET T0
!    call lltbbackgroundHtofa(vpi,a=1._dl,t=t0)

    ! Both at r=0 and r=inf, R'= a, a'=0
    rho_out = 1._dl ! times Mtilde^2 times Mp^2

    !-------!
    ! set a0_in
    ! vpi%xyz(1) = 8._dl * pi * Mt2 / 3._dl
!    vpi%xyz(2) = 2._dl * Mt2 * ( kb + kmax )
    ! vpi%xyz(3) = Lambdaover3

    ! It is very very well possible that for this t0 and this kb+kmax
    ! the metric at this radius r=0 has already recollapsed to its big crunch.
    ! So before we ask for a(t), which may give an error (t>t_collapse)
    ! we must check what is the crunch time:

!    thissubcase=GetOmegaCase(vpi%xyz)
!    call lltb_getallroots(vpi,thissubcase)
!    turnovera=min(vpi%smallestrealroot,vpi%smallestlinroot)

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
!             call lltb_getallroots(vpi, GetOmegaCase(vpi%xyz,turnovera))
!             turnovera = min(vpi%smallestrealroot,vpi%smallestlinroot)
             ! End trouble.


!    if(turnovera<1.e20_dl)then
!       call lltbbackgroundHtofa(vpi,a=turnovera,t=turnovertime)
!    else
!       turnovertime=infinity !1.e30_dl
!    end if
    
!          if(turnovertime>t0)then

!    if(2._dl*turnovertime>t0)then


!       call lltbbackgroundaofHt(vpi,thisa,t0)

    
       !-------!
       ! set H0_in
!       Hlocal = HoverH0(vpi%xyz,thisa)
!       if(t0>vpi%tturn)Hlocal=-Hlocal
       !-------!
       ! set d0
!       rho_in = 1._dl / thisa**3 ! times Mtilde^2 times Mp^2
!       delta0 = rho_in / rho_out - 1._dl
!    else


       !-------!
       ! set H0_in
!       Hlocal = -infinity !-1.e30_dl

       
       !-------!
       ! set d0
!       delta0 = infinity !1.e30_dl
!    end if

  contains
  
    function dumfunc(r)
      implicit none
      real(Dl), intent(in) :: r
      real(dl) :: dumfunc
      
      dumfunc = 0.d0
    end function dumfunc

    function dumfunc2(r)
      implicit none
      real(Dl), intent(in) :: r
      real(dl) :: dumfunc2
      
      if(r<0.9d0)then
        dumfunc2 = kmax + kb
      else
        dumfunc2 = kb
      end if        
    end function dumfunc2

    function dumfunc3(r)
      implicit none
      real(Dl), intent(in) :: r
      real(dl) :: dumfunc3
      
      ! IMPORTANT: these radii are fixed by hand, and
      ! must be such that initime finds the right kb and kmax.
      if(r<0.5d0)then
        dumfunc3 = 0.d0
      else if(r<0.9d0)then
        dumfunc3 = kmax
      else
        dumfunc3 = 0.d0
      end if        
    end function dumfunc3

  end subroutine lltb_normalization


  ! lltb_feedback: set verbosity of lltb module (0=default, higher is more)
  subroutine lltb_feedback(number)
    use lltb, only : oldlltb => lltb_feedback
    integer, intent(in) :: number
    
    lfeedback=number
write(*,*)number    
    call oldlltb(number)
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
    use elliptics , only : Elliptics_precisiongoal
    real(dl), intent(in) :: goal
    real(dl) :: orig_goal
    real(dl) :: power
    ! We initialize with precision goal = 1.d-16
    
!    orig_goal = accuracy_p(dl) ! parameter accuracy
    
    if(goal < orig_goal)then
       write(*,'(A,ES14.5,A,ES14.5,A)')'You demand precisiongoal ',goal,' which is smaller than the hardcoded lowerlimit of ',orig_goal,'. That will not work.'
       stop 'Too high precisiongoal demanded. LLTBBackground.'
    end if

    power = log(goal)/log(orig_goal)

    ! Here we go. Be careful, this list must be complete. Compare to module precision.
!    accuracy=accuracy_p**power
!    lin_expansion=lin_expansion_p**power
!    acceptance_pure=acceptance_pure_p**power
!    acceptance=acceptance_p**power
!    acceptance_pole=acceptance_pole**power
!    accuracy_downscale=accuracy_downscale_p*power ! Note times, not to the power
!    rel_acceptance = rel_acceptance_p**power

    call Elliptics_precisiongoal(goal)

  end subroutine lltb_precisiongoal_dl

  subroutine lltb_attach_to_error(msg)
    use lltb, only : oldlltb => lltb_attach_to_error
    implicit none
    character(len=*), intent(in) :: msg

    call wlltb_set_error(attachement=msg)
    call oldlltb(msg)

  end subroutine lltb_attach_to_error
        

end module wlltb
