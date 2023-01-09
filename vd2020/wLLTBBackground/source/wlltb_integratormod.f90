!*
! *  wlltb_integratormod.f90
! *  wLLTBBackground
! *
! *  Created by Wessel Valkenburg on 02/12/2011.
! *  Copyright 2011 Wessel Valkenburg. All rights reserved.
! *
! */

! 
module wlltb_integratormod
  use wlltb_constants
  use wlltb_types
  use wlltb_dvode
  use wlltb_errors
  implicit none

  private
  
  real(dl) :: Mt2, radius, w_eos, Eofroverr2, Eoverr2prime, rhoMbartimesa2Yp, adotsign, rhowout0, glob_t_asked, glob_t_bar
  real(dl) :: omegal, omegam, omegak, H0_inf_glob
  logical :: myfeedback=.false., stiff=.false.
  integer, parameter :: nevents=4
  real(dl), parameter :: adotepsilon = 1.d-8, myrtol=1.d-8, myatol=1.d-12
!  real(dl), parameter :: adotepsilon = 1.d-4, myrtol=1.d-4, myatol=1.d-8


  public :: wlltb_integrator, wlltb_integratortturnonly, wlltb_derivs, radotdotprime, rhowout

contains

  subroutine initintegrationconstants(H0_inf, Lambda, w, radfuncs, Mtilde, kb, &
						r, t, tbar, rhoMbarin, y)
    implicit none
    real(dl), intent(in) :: H0_inf, Lambda, w, Mtilde, r, t, tbar, rhoMbarin, kb
    type(funcscontainer) :: radfuncs
    real(dl), intent(in) :: y(:)

!    Mt2 = Mtilde**2
    Mt2 = sq(Mtilde)
    radius = r
    w_eos = w
    Eofroverr2 = Mt2 * radfuncs%kofr(r)
    Eoverr2prime = Mt2 * radfuncs%dkdr(r) 
!    rhoMbartimesa2Yp = rhoMbarin * norms%abar**3
    rhoMbartimesa2Yp = rhoMbarin * cub(norms%abar)

    adotsign = sign(1.d0,y(index%adot))
             
!    omegal =   Lambda /( 3._dl * H0_inf**2 )
    omegal =   Lambda /( 3._dl * sq(H0_inf) )
!    omegak = 2.d0 * Mtilde**2 / H0_inf**2 * kb
    omegak = 2.d0 * sq(Mtilde) / sq(H0_inf) * kb
    omegam = 1.d0 - omegak - omegal
    H0_inf_glob = H0_inf

    rhowout0 = Lambda / 8.d0 / pi

  end subroutine initintegrationconstants

  subroutine wlltb_integratortturnonly(H0_inf, Lambda, w, radfuncs, Mtilde, kb, &
						r, t, tbar, rhoMbarin, y,tturn)
    implicit none
    real(dl), intent(in) :: H0_inf, Lambda, w, Mtilde, r,t, tbar, rhoMbarin, kb
    real(dl), intent(inout) :: y(:), tturn
    type(funcscontainer) :: radfuncs
    ! DVODE:
    integer :: neq, itask, istats(31), iout, ierror, istate
    real(dl) :: atol(index%max), rstats(22), rtol, t_dvode, tout_dvode, dumdy(index%max), g0(nevents), dummy
    TYPE (VODE_OPTS) :: OPTIONS
    ! END DVODE
    logical :: repeat
    integer :: i
  
    if(size(y)/=index%max)call wlltb_error("vector y wrong size in wlltb_integratortturnonly.")
    if(t<tbar)call wlltb_error("asking t < tbar in wlltb_integratortturnonly.")
    
    call initintegrationconstants(H0_inf, Lambda, w, radfuncs, Mtilde, kb, &
						r, t, tbar, rhoMbarin, y)
    
    ! Calling DVODE:
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL, &
!        RELERR=RTOL,USER_SUPPLIED_JACOBIAN=.TRUE.)
    
    neq = index%max
    RTOL = myrtol
    ATOL = myatol
    ATOL(index%adot) = 1.d-1
    ITASK = 1
    ISTATE = 1
    t_dvode = tbar
    tout_dvode = t
    
    ! reset cache:
    if(derivs_cache(t,y,dumdy,.true.))then
      continue
    else
      continue
    end if

    OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, &
        RELERR=RTOL,NEVENTS=1) ! Nevents is size of gfun's result.
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., ABSERR=ATOL(1), &
!        RELERR=RTOL,NEVENTS=1)
    repeat = .true.
    i = 1
    
    stiff = .false.
    
    do while (repeat)
      CALL DVODE_F90(wlltb_derivs,NEQ,Y,t_dvode,tout_dvode,ITASK,ISTATE,OPTIONS,G_FCN=tturnGFUN)
      ! DVODE done.
      if(ISTATE==1)then ! nothing done: t=tout
        repeat=.false.
      else if(ISTATE==2)then ! success: tout reached.
        repeat=.false.
        tturn = 1.d30
      else if(ISTATE==3)then ! semi-success: root reached
        repeat=.false.
        tturn = t_dvode
        !write(0,*)'tturn found', tturn, t
      else if (ISTATE<1) then
        call wlltb_error("DVODE background Integration did not converge while searching tturn.")
      end if
      if(.not.repeat)exit
      i=i+1
      if(i>500)call wlltb_error('i>500 in wlltb_integratortturnonly')
    end do
    
    ! dvode cleanup:
    call RELEASE_ARRAYS
    if(associated(OPTIONS%rtol))deallocate(OPTIONS%rtol)
    if(associated(OPTIONS%atol))deallocate(OPTIONS%atol)

!    call wlltb_derivs(sizE(y),tout_dvode,y,dumdy)
!    y(index%raprime) = y(index%raprime) * dumdy(index%aLTB)

  contains
  
      ! Root finding:
    SUBROUTINE tturnGFUN(NEQ, T, Y, NG, GOUT)
      implicit none
      integer, intent(in) :: neq, NG
      real(dl), intent(in) :: t, y(NEQ)
      real(dl), intent(out) :: GOUT(NG)
      real(dl) :: dy(neq)
      real(dl) :: reladot, thissign
    
      if(ng/=1)call wlltb_error('NG /= 1 in tturngfun in wlltb_integratormod.f90')

      call wlltb_derivs(neq, t, y, dy)

      gout(1) =  dy(index%altb)
      !write(*,*)gout,t
    end subroutine tturnGFUN
  end subroutine wlltb_integratortturnonly

  subroutine wlltb_integrator(H0_inf, Lambda, w, radfuncs, Mtilde, kb, &
						r, t, tbar, rhoMbarin, y)
    implicit none
    real(dl), intent(in) :: H0_inf, Lambda, w, Mtilde, r, t, tbar, rhoMbarin, kb
    real(dl), intent(inout) :: y(:)
    type(funcscontainer) :: radfuncs
    ! DVODE:
    integer :: neq, itask, istats(31), iout, ierror, istate
    real(dl) :: atol(index%max), rstats(22), rtol, t_dvode, tout_dvode, dumdy(index%max), g0(nevents), dummy
    TYPE (VODE_OPTS) :: OPTIONS
    ! END DVODE
    logical :: repeat
    integer :: i
  
    if(size(y)/=index%max)call wlltb_error("vector y wrong size in wlltb_integrator.")
 if(t<tbar)then
    write(0,*)'t,tbar:',t,tbar
    call wlltb_error('t<tbar')
  end if
    if(t<tbar)call wlltb_error("asking t < tbar in wlltb_integrator.")
    
    call initintegrationconstants(H0_inf, Lambda, w, radfuncs, Mtilde, kb, &
						r, t, tbar, rhoMbarin, y)
    
    if(abs(radfuncs%dkdr(r))<1.d-14 .and. abs(radfuncs%kofr(r)-kb)<1.d-14)then
!      call integrate_FLRW(t,tbar,norms%abar,y)
!      return
    end if
 if(t<tbar)then
    write(0,*)'t,tbar:',t,tbar
    call wlltb_error('t<tbar')
  end if
  
  glob_t_asked=t
  glob_t_bar=tbar
    ! Calling DVODE:
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL, &
!        RELERR=RTOL,USER_SUPPLIED_JACOBIAN=.TRUE.)
    
    neq = index%max
    RTOL = myrtol
    ATOL = myatol
    ATOL(index%adot) = 1.d-1
    ITASK = 1
    ISTATE = 1
    t_dvode = tbar
    tout_dvode = t
    
    ! reset cache:
    if(derivs_cache(t,y,dumdy,.true.))then
      continue
    else
      continue
    end if

    OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, &
        RELERR=RTOL,NEVENTS=nevents) ! Nevents is size of gfun's result.
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., ABSERR=ATOL(1), &
!        RELERR=RTOL,NEVENTS=1)
    repeat = .true.
    i = 1
    
    stiff = .false.
    
    do while (repeat)
      CALL DVODE_F90(wlltb_derivs,NEQ,Y,t_dvode,tout_dvode,ITASK,ISTATE,OPTIONS,G_FCN=GFUN)
      ! DVODE done.
      if(ISTATE==1)then ! nothing done: t=tout
        repeat=.false.
      else if(ISTATE==2)then ! success: tout reached.
        repeat=.false.
      else if(ISTATE==3)then ! semi-success: root reached, continue integration.
        repeat=.true.

        call GFUN(NEQ, t_dvode, Y, Nevents, G0)
        call wlltb_derivs(sizE(y),t_dvode,y,dumdy)
        if(abs(g0(1))<myrtol)then  
          if(.not.stiff)then
            stiff=.true.
            ISTATE = 1
            rtol = myrtol * 1.d-4
            ATOL = myatol * 1.d-2
            ATOL(index%adot) = myatol * 1.d-4
!            OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, &
!              RELERR=RTOL,NEVENTS=nevents) ! Nevents is size of gfun's result.
            OPTIONS = SET_NORMAL_OPTS(SPARSE_J=.TRUE., ABSERR=ATOL(1), &
                RELERR=RTOL,NEVENTS=nevents)
!            y(index%adot)=abs(adotsquared(y))**0.5d0
            y(index%adot)=sqrt(abs(adotsquared(y)))
            dummy = radotprimelimit(dumdy(index%raprime))
!            write(*,*)'switching to second order',y(index%adot),adotsquared(y)
!            write(*,*)'g0:',g0
!            write(*,*)'y',y
          else 
            stiff=.false.
            adotsign = sign(1.d0,dumdy(index%altb))
            ISTATE = 1
            rtol = myrtol 
            ATOL = myatol
            ATOL(index%adot) = 1.d-1
            OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, &
                RELERR=RTOL,NEVENTS=nevents) ! Nevents is size of gfun's result.
!            write(*,*)'switching back to first order', adotsign,y(index%adot),dumdy(index%altb)
!            write(*,*)'g0:',g0
!            write(*,*)'y',y

!            if(adotsign>0)stop 'wrong adotsign'

          end if              
        end if
!        if(i>500)call wlltb_error('i>500 in wlltb_integrator')

        if(any(abs(g0(2:3))<myrtol) .or. i>500)then ! big crunch
          if(lfeedback>1)write(0,*)'Collapsed:',y
          if(lfeedback>1)write(0,*)'g0(2:3):',g0(2:3)
          if(lfeedback>1)write(0,*)'i:',i
          y(index%altb) = 0.d0
          y(index%adot) = 1.d31
          y(index%for3) = 1.d31
!          write(0,*)'Collapsed.'
          exit
        end if
if(lfeedback>3)write(*,'(20ES15.6)')y
if(lfeedback>3)write(*,'(20ES15.6)')dumdy
!if(lfeedback>3)write(*,'(20ES15.6)')g0
      else if (ISTATE<1) then
          if(lfeedback>2)then
            call GFUN(NEQ, t_dvode, Y, Nevents, G0)
            write(*,'(A,20ES14.5)')'roots:',g0
            write(*,'(A,20ES14.5)')'y:',y
            call wlltb_derivs(sizE(y),t_dvode,y,dumdy)
            write(*,'(A,20ES14.5)')'dy:',dumdy
          end if
          y(index%altb) = 0.d0
          y(index%adot) = 1.d31
          y(index%for3) = 1.d31
          write(0,*)"DVODE background Integration did not converge."
!          call wlltb_error("DVODE background Integration did not converge.")
          call wlltb_error_nocrash("DVODE background Integration did not converge.")
          exit
      end if
      if(.not.repeat)exit
      i=i+1
    end do

    ! dvode cleanup:
    call RELEASE_ARRAYS
    if(associated(OPTIONS%rtol))deallocate(OPTIONS%rtol)
    if(associated(OPTIONS%atol))deallocate(OPTIONS%atol)

!    call wlltb_derivs(sizE(y),tout_dvode,y,dumdy)
!    y(index%raprime) = y(index%raprime) * dumdy(index%aLTB)

  end subroutine wlltb_integrator
  
  function adotsquared(y)
    implicit none
    real(dl), intent(in) :: y(index%max)
    real(dl) :: adotsquared

    adotsquared = 2.d0 * y(index%for3) / y(index%aLTB) + 2.d0 * Eofroverr2

  end function adotsquared

  subroutine integrate_FLRW(t,tbar,abar,y)
    use wlltb_numtools
    implicit none
    real(dl), intent(in) :: t, tbar, abar
    real(dl), intent(inout) :: y(:)
    real(dl) :: amax, amin, da, nextt, lastt, mda, testa
    integer :: i, imax = 500
    real(dl), save :: lasta = norms%a0
    amax = norms%a0*1.1
    amin = abar
    testa = lasta ! educated guess: often this routine is called for subsequential times in a geodesic integration..
    ! Nice: it saves up to a factor of 10 for small t.
    
    nextt = wlltb_rombint(thisdtda,abar,testa,myatol)
    i=2
    da=0
    do while (abs(nextt/(t-tbar)-1.d0)>myrtol)
      i=i+1
      if (nextt>t-tbar) then
        amax = testa
      else
        amin = testa
      end if
      da = (amax + amin)*0.5d0 - testa
      if(i>1)mda = (t-tbar - nextt)/thisdtda(testa)
      if(mda+testa<amax .and. mda+testa>amin)da=mda
      testa = testa + da
      nextt = wlltb_rombint(thisdtda,abar,testa,myatol)
      if(i>imax)then
        call wlltb_error('rombint inversion is hanging in wlltb_integrator.')
        stop
      end if
!      write(*,*)testa,amax,amin,t-tbar,nextt,abs(nextt/(t-tbar)-1.d0),da
    end do
    
    lasta = testa
    
!write(*,*)testa,t, i
    y=0.d0
    y(index%aLTB) = testa
    y(index%adot) = 1/thisdtda(testa)
    y(index%raprime) = 0.d0
    y(index%aout) = testa
!    y(index%For3) = 0.5d0*testa*aoutdotsquared(y) - testa * H0_inf_glob**2 * omegak * 0.5d0
    y(index%For3) = 0.5d0*testa*aoutdotsquared(y) - testa * sq(H0_inf_glob) * omegak * 0.5d0
    y(index%rhow) = RHOWOUT(y)
    
    contains
    
      function thisdtda(a)
        implicit none
        real(dl) :: thisdtda
        real(dl), intent(in) :: a
        real(dl) :: y(index%max)
        
        y(index%aout) = a
        
        thisdtda = aoutdotsquared(y)**(-0.5d0)
        
      end function thisdtda
      
  end subroutine integrate_FLRW
  

  function aoutdotsquared(y)
    implicit none
    real(dl), intent(in) :: y(index%max)
    real(dl) :: aoutdotsquared
    real(dl) :: a, a0
    
    a0 = norms%a0
    a = y(index%aout)/a0

!    aoutdotsquared = H0_inf_glob**2 * (omegam / a + omegak + omegal / (a**((1+3*w_eos)) ))
    aoutdotsquared = sq(H0_inf_glob) * (omegam / a + omegak + omegal / (a**((1+3*w_eos)) ))
!    aoutdotsquared = sq(H0_inf_glob) * (omegam / a + omegak + omegal * sq(a) )

    if(aoutdotsquared<0.d0)call wlltb_error("Aoutdotsquared < 0: collapsing outside cosmology.")
  
  end function aoutdotsquared
  
  function eff_adot(y)
    implicit none
    real(dl), intent(in) :: y(index%max)
    real(dl) :: eff_adot
    real(dl) :: adot2
    
    adot2 = adotsquared(y) 

    if(adot2 <= 0.d0 .or. isnan(adot2) .and. .not.stiff)then 
!      write(*,*)'Turn point:',y(index%aLTB), adot2,radius,adotsign
!      eff_adot = -adotsign*abs(ADOT2)**0.5d0
      eff_adot = -adotsign*sqrt(abs(ADOT2))
    else
!      eff_adot = adotsign*((adot2))**0.5d0
      eff_adot = adotsign*(sqrt(adot2))
    end if

    if(stiff)then
      eff_adot = y(index%adot)
!      write(*,*)'stiff', adot
    end if
!      adot = y(index%adot)

  end function eff_adot
  
  
  function adotdot(y)
    implicit none
    real(dl), intent(in) :: y(index%max)
    real(dl) :: adotdot
    real(dl) :: Foverr3
    
    Foverr3 = y(index%for3)

    adotdot = -4.d0 * pi * y(index%aLTB) * wofrt(y) * y(index%rhow) - Foverr3 / y(index%aLTB)**2
!    adotdot = -4.d0 * pi * y(index%aLTB) * wofrt(y) * y(index%rhow) - Foverr3 / sq(y(index%aLTB))
!    adotdot = -4.d0 * pi * y(index%aLTB) * w_EOS * RHOWOUT - Foverr3 / y(index%aLTB)**2


  end function adotdot
  

  function Fprimeoverr2(y)
    implicit none
    real(dl), intent(in) :: y(index%max)
    real(dl) :: Fprimeoverr2
    real(dl) :: Yprime

    Yprime = (y(index%aLTB) + y(index%raprime))
    Fprimeoverr2 = 4.d0 * pi * (rhoMbartimesa2Yp + y(index%rhow) * sq(y(index%aLTB)) * Yprime)

  end function Fprimeoverr2

  function rhowout(y)
    implicit none
    real(dl), intent(in) :: y(index%max)
    real(dl) :: rhowout
    
    rhowout = rhowout0 * (norms%a0 / y(index%aout))**(3*(1+w_eos))
    
  end function rhowout

  function wofrt(y)
    implicit none
    real(dl), intent(in) :: y(index%max)
    real(dl) :: wofrt

    wofrt=w_eos * (rhowout(y)/y(index%rhow))

  end function wofrt
  

  subroutine wlltb_derivs(neq, t, y, dy)
    use wlltb_funcs, only: allposdefsarepos
    implicit none
    integer, intent(in) :: neq
    real(dl), intent(in) :: t, y(neq)
    real(dl), intent(out) :: dy(neq)
    ! local:
    real(dl) :: adot, adot2, Foverr3, Ha, Hr, radotprime, Fprime
    real(dl) :: r, rhom, rhow, Yprime, Ydotprime, mFprimeoverr2, wrt, mrhowout
    logical :: alert
    
!      write(*,*)y,dy,'...'
    if( derivs_cache(t,y,dy) ) then
      return
    end if
!    if(y(index%aLTB)<norms%abar)then
!      write(*,*)y,dy,'...'
!      write(*,*)'t,tbar:',glob_t_asked,glob_t_bar
!      stop 'a<abar.'
!    end if
!    if(abs(y(index%aLTB)/norms%abar-1.d0)<1.d-16)then
!      write(*,*)y,dy,'...'
!      write(*,*)'t,tbar:',glob_t_asked,glob_t_bar
!      stop 'a==abar.'
!    end if
    
    if(size(dy)/=size(y))call wlltb_error("dy and y have different sizes in wlltb_derivs")
    if(size(y)/=index%max)then
      write(0,*)"Size of y:",size(y)
      write(0,*)"NEQ:",neq
      call wlltb_error("vector y wrong size in wlltb_derivs")
    end if
    if(any(isnan(y)))then
!      write(0,*)"y:",y
!      call wlltb_error("isnan(y) in wlltb_derivs")
      alert=.true.
    end if

    r=radius
    mrhowout = rhowout(y)
    !    rhowout = omegal * (a0 / abar)**(3*(1+w)) * H0_inf**2 * 3.d0 / 8.d0 / pi
    wrt=wofrt(y)! w_eos * (mrhowout/y(index%rhow))
    alert = .false.

    if(.not.allposdefsarepos(y))alert=.true.

    adot = eff_adot(y)
!1
    dy(index%aLTB) = adot
                 
!if(.not.stiff.and. adotsign<0.)write(*,*)adot/y(index%adot)-1.d0,adot,y(index%adot)
                 
    Ha = adot / y(index%aLTB)
    

!2
    dy(index%adot) = adotdot(y)
    Foverr3 = y(index%for3)

    
    Yprime = (y(index%aLTB) + y(index%raprime))
  
!    mFprimeoverr2 = 4.d0 * pi * (rhoMbartimesa2Yp + y(index%rhow) * y(index%aLTB)**2 * Yprime)
    mFprimeoverr2 = 4.d0 * pi * (rhoMbartimesa2Yp + y(index%rhow) * sq(y(index%aLTB)) * Yprime)
!    radotprime = 1.d0 / adot * &
!        (   mFprimeoverr2 / y(index%aLTB) &
!          -  Foverr3 / y(index%aLTB)**2 * &
!              (3.d0 * y(index%aLTB) + y(index%raprime)) &
!          + r * Eoverr2prime &
!        )
    radotprime = 1.d0 / adot * &
        (   mFprimeoverr2 / y(index%aLTB) &
          -  Foverr3 / sq(y(index%aLTB)) * &
              (3.d0 * y(index%aLTB) + y(index%raprime)) &
          + r * Eoverr2prime &
        )
!write(*,*)radotprime
if(.false.)then
! Or equivalently:
    radotprime = 1.d0 / adot * &
        (   mFprimeoverr2 / y(index%aLTB) &
          -  3.d0 * Foverr3 / y(index%aLTB) &          
          + y(index%raprime) / y(index%aLTB) * Eofroverr2  &
          + r * Eoverr2prime &
        ) & 
        - y(index%raprime) * adot / 2.d0  / y(index%aLTB)

write(*,*)radotprime
! Or:
    radotprime = 1.d0 / adot * &
        (   mFprimeoverr2 / y(index%aLTB) &
          + y(index%raprime) * Eofroverr2 / y(index%aLTB) &
          + 3.d0 * Eofroverr2  &
          + r * Eoverr2prime &
        ) & 
        - adot * y(index%raprime) / 2.d0  / y(index%aLTB) - adot * 1.5d0
        
write(*,*)radotprime

! Or:
    radotprime =  1.d0 / adot * &
        (   (4.d0 * pi * rhoMbartimesa2Yp )/ y(index%aLTB) &
         + ( 4.d0 * pi * y(index%rhow) * y(index%aLTB)**2 * Yprime) / y(index%aLTB) &
          + y(index%raprime) * Eofroverr2 / y(index%aLTB) &
          + 3.d0 * Eofroverr2  &
          + r * Eoverr2prime &
        ) & 
        - adot * y(index%raprime) / 2.d0  / y(index%aLTB) - adot * 1.5d0
        
  write(17,*)t, &
           mFprimeoverr2 / y(index%aLTB) &
          ,-  3.d0 * Foverr3 / y(index%aLTB) &          
          , y(index%raprime) / y(index%aLTB) * Eofroverr2  &
          , r * Eoverr2prime &
        ,- y(index%raprime) * adot / 2.d0  / y(index%aLTB), adot, radotprime, y(index%adot)
  
  
  write(18,*)t,y(index%adot),(   mFprimeoverr2 / y(index%aLTB) &
          + y(index%raprime) * Eofroverr2 / y(index%aLTB) &
          + r * Eoverr2prime &
          + 3.d0 * Eofroverr2  &
        ), &
        mFprimeoverr2 / y(index%aLTB), &
          + y(index%raprime) * Eofroverr2 / y(index%aLTB), &
          + r * Eoverr2prime, &
          + 3.d0 * Eofroverr2 
end if

    if (isnan(radotprime) .or. (abs(radius)<1.d-30) .or. (abs(radotprime)>1.d30))radotprime = 0.d0
!3  
    if(stiff)then
      if(abs(radotprime/radotprimelimit()-1.d0)>1.d-2) radotprime =  radotprimelimit()
      if(abs(radotprime/radotprimelimit()-1.d0)<1.d-2) radotprime =  radotprimelimit(radotprime) ! then set limit again 
 !     write(*,*)radotprime,dy(index%raprime), radotprimelimit(), adot 
    end if
    dy(index%raprime) =  radotprime 
      
      
!write(21,*)t,adot,radotprime,adotsquared(y), y(index%for3), y(index%aLTB), Eofroverr2, dy(index%raprime)

!4
!    dy(index%f) = -4.d0 * pi * r**3 * y(index%aLTB)**2 * adot * wrt * y(index%rhow)
!    dy(index%for3) = -4.d0 * pi * y(index%aLTB)**2 * adot * w_eos*mrhowout 
    dy(index%for3) = -4.d0 * pi * sq(y(index%aLTB)) * adot * w_eos*mrhowout 
!    dy(index%for3) = -4.d0 * pi * y(index%aLTB)**2 * adot *  wrt * y(index%rhow)

    Ydotprime = adot + radotprime
  
    Ha = adot / y(index%aLTB)
    Hr = Ydotprime / Yprime

    if(isnan(Hr) .or. abs(Hr)>1.d30)then
      if(lfeedback>0)then
        write(0,*)'setting hr to zero:',Hr,Ha
      end if
      Hr = 0.d0
    end if
!5
    dy(index%rhow) = -(1.d0 + wrt)*y(index%rhow)*(Hr + 2.d0 * Ha)
    if(dy(index%rhow)>0.d0 .and. abs(w_eos+1.d0)>1.d-15)then
      if(lfeedback>0)then
        write(0,*)-(1.d0 + wrt),y(index%rhow),(Hr + 2.d0 * Ha)
        writE(0,*)Hr,Ha
        write(0,*)'E(r),E''(r),w,w(r,t):',Eofroverr2*r**2,Eoverr2prime *r**2, w_eos,wrt
        write(0,*)'dy(rhow) > 0. Wrong?',dy(index%rhow)
!!        stop
      end if
    end if
    if(abs(dy(index%rhow))>1.d29)then
      if(lfeedback>0)write(0,*)'wlltb capping it'
      dy(index%rhow)=sign(1.d29,dy(index%rhow))
    end if
!    if(abs(dy(index%rhow))>1.d29)then
!      write(0,*)'drhow trouble radotprimelimit:',radotprimelimit(),radotprime
!      write(0,*)'drhow trouble wrt:',wrt
!      write(0,*)'drhow trouble y(index%rhow):',y(index%rhow)
!      write(0,*)'drhow trouble Hr:',Hr
!      write(0,*)'drhow trouble Ha:',Ha
!      write(0,*)'drhow trouble Yprime:',Yprime
!      write(0,*)'drhow trouble Ydotprime:',Ydotprime
!      write(0,*)'drhow trouble Yprime/Ydotprime:',Yprime/Ydotprime
!    end if
!6
!    dy(index%aout) = aoutdotsquared(y)**0.5d0
    dy(index%aout) = sqrt(aoutdotsquared(y))

    if(alert)then
      dy(index%adot)=0.d0
      dy(index%raprime)=0.d0
      dy(index%for3)=0.d0
      dy(index%rhow)=0.d0
dy=0.d0
    end if

    if(any(isnan(dy)))then
      write(0,*)"dy:",dy
      write(0,*)"y:",y
      call wlltb_error("NaNs in dy(:) in wlltb_derivs.")
    end if
    if(any(abs(dy)>1.d30).and. (.not.abs(y(index%altb)/abs(dy(index%altb)) + myatol * sign(1.d0,dy(index%altb))) < myrtol ))then
      write(0,*)"dy:",dy
      write(0,*)"y:",y
      write(0,*)"gout(2):",y(index%altb)/abs(dy(index%altb)) + myatol * sign(1.d0,dy(index%altb))
      write(0,*)"stiff?", stiff
      call wlltb_error("Infs in dy(:) in wlltb_derivs.")
    end if
        
    if( derivs_cache(t,y,dy) ) then
      continue
    else
      continue
    end if
       
!write(*,*)y,dy,'.'       
!write(*,*)Hr,Ydotprime,'/',Yprime
!write(*,*)Ha, adot,'/', y(index%aLTB)

  end subroutine wlltb_derivs
  
  
  function Fdotprimeoverr2(y)
    implicit none
    real(dl), intent(in) :: y(:)
    real(dl) :: Fdotprimeoverr2
    real(dl) :: dy(index%max), Yprime, Ydotprime
    real(dl) :: dumt, adot, wrhow, a
    dumt = 1.d0
    
    call wlltb_derivs(index%max,dumt,y,dy)
    
    Yprime = (y(index%aLTB) + y(index%raprime))
    Ydotprime = (dy(index%aLTB) + dy(index%raprime))
    adot = dy(index%aLTB)
    a = y(index%aLTB)
    wrhow = rhowout(y) * w_eos
    
!    Fdotprimeoverr2 = -4.d0 * pi * wrhow * ( 2.d0 * a * adot * Yprime + a**2 * Ydotprime  )
    Fdotprimeoverr2 = -4.d0 * pi * wrhow * ( 2.d0 * a * adot * Yprime + sq(a) * Ydotprime  )
  
  end function Fdotprimeoverr2
  
  
 function radotdotprime(y)
    implicit none
    real(dl), intent(in) :: y(index%max)
    real(dl) :: radotdotprime
    real(dl) :: Foverr3, a
    
    a = y(index%aLTB)
    Foverr3 = y(index%for3)

!    radotdotprime = -4.d0 * pi * y(index%raprime) * w_EOS * RHOWOUT(y) &
!                    - Fprimeoverr2(y) / a**2 &
!                    + 3.d0 * Foverr3 / a**2 &
!                    + 2.d0 * Foverr3 / a**3 * y(index%raprime) 
    radotdotprime = -4.d0 * pi * y(index%raprime) * w_EOS * RHOWOUT(y) &
                    - Fprimeoverr2(y) / sq(a) &
                    + 3.d0 * Foverr3 / sq(a) &
                    + 2.d0 * Foverr3 / cub(a) * y(index%raprime) 
!    adotdot = -4.d0 * pi * y(index%aLTB) * w_EOS * RHOWOUT - Foverr3 / y(index%aLTB)**2


  end function radotdotprime

  
  ! Root finding:
  SUBROUTINE GFUN(NEQ, T, Y, NG, GOUT)
    implicit none
    integer, intent(in) :: neq, NG
    real(dl), intent(in) :: t, y(NEQ)
    real(dl), intent(out) :: GOUT(ng)
    real(dl) :: dy(neq)
    real(dl) :: reladot, thissign, adotepsiloncorrection
    
    if(ng/=nevents)call wlltb_error('NG /= Nevents in gfun in wlltb_integratormod.f90')

    call wlltb_derivs(neq, t, y, dy)
  
!     adotsquared = 2.d0 * y(index%for3) / y(index%aLTB) + 2.d0 * Eofroverr2

!    reladot = dy(index%altb)/min(abs(Eofroverr2),abs(y(index%for3) / y(index%aLTB)))
    reladot = eff_adot(y)/min(abs(Eofroverr2),abs(y(index%for3) / y(index%aLTB)))
    thissign = sign(1.d0,dy(index%altb))
    adotepsiloncorrection=1.d0
!    if(w_EOS==-1.d0)adotepsiloncorrection=1.d-4
    gout(1) =  reladot - thissign * adotepsilon * adotepsiloncorrection
    if(stiff)gout(1)=1.d0
!    gout(2) = reladot + adotepsilon
    
    gout(2) = y(index%altb)/abs(dy(index%altb)) + myatol * sign(1.d0,dy(index%altb))
    if(isnan(gout(2)).or.abs(gout(2))>1.d30)gout(2)=1.d30
!    write(*,*)gout(2),y(index%altb),adotsquared(y)

    gout(3) = y(index%rhow)
    if(y(index%rhow)<0.d0)then
!      write(*,*)'negative rhow!'
      gout(3)=0.d0
    end if
    if(w_eos==-1.d0)then
      gout(3)=1.d0 ! because then we couldn't care less about y(3)
    end if
    
    if(any(isnan(y)))then
      gout(2:3)=0.d0
    end if

    if(w_eos/=-1.d0)then
      gout(4) = dy(index%rhow)
    else
      gout(4)=1.d0
    end if
!  write(*,*)gout(4),'gfun'
  end subroutine GFUN
  
  function derivs_cache(t,y,dy,reset)
    implicit none
    real(dl), intent(in) :: t, y(:)
    real(dl), intent(inout) :: dy(:)
    logical :: derivs_cache
    logical, optional :: reset
    real(dl), save :: mt, my(index%max),mdy(index%max)
    
    if(present(reset))then
      if(reset)then
        mt=-1.d30
        my=-1.d30
        mdy=-1.d30
        derivs_cache = .false.
        return
      end if
    end if
    
    if((t==mt) .and. (all(y==my)))then
      dy=mdy
      derivs_cache = .true.
    else
      mdy=dy
      derivs_cache = .false.
    end if
    
            
  end function derivs_cache
  
  function radotprimelimit(radotprime)
    implicit none
    real(dl), intent(in), optional :: radotprime
    real(dl) :: radotprimelimit
    real(dl), save :: mem
    
    if(present(radotprime))then
      mem = radotprime
    end if
    
    radotprimelimit = mem
    
  end function radotprimelimit
  
  function sq(x)
    implicit none
    real(dl), intent(in) :: x
    real(dl) :: sq
    
    sq = x*x
  
  end function sq
    
  function cub(x)
    implicit none
    real(dl), intent(in) :: x
    real(dl) :: cub
    
    cub = x*x*x
  
  end function cub
    

end module wlltb_integratormod
