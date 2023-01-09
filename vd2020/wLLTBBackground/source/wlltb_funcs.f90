!*
! *  wlltb_funcs.f90
! *  wLLTBBackground
! *
! *  Created by Wessel Valkenburg on 29/11/2011.
! *  Copyright 2011 Wessel Valkenburg. All rights reserved.
! *
! */

! 
module wlltb_funcs
  use wlltb_constants
  use wlltb_errors
  implicit none
  
  private
  type outsidecosmology
    real(dl) :: omegal, omegak, omegam, H0, tbar, t0
  end type outsidecosmology
  
  type(outsidecosmology) :: prevCosmo
  logical :: firstcall = .true.
  
 interface operator (==)
    module procedure outsidecosmologycompare
  end interface 

  
  public :: wlltb_initime, allposdefsarepos

contains

!      call wlltb_initime(H0_inf, Lambda, w, radfuncs, Mtilde, &
!            t0, tbar, abar, Hbarout, y)
  subroutine wlltb_initime(H0_inf, Lambda, w, radfuncs, Mtilde, r, t0, tbar, abar, Hbarout, rhombar, rhoMout, kb, FMbaroverout, y)
    use wlltb_types
    use wlltb_numtools
    implicit none
    real(dl), intent(in) :: H0_inf, Lambda, w, Mtilde, r
    real(dl), intent(out) :: t0, tbar, abar, Hbarout, rhombar, rhoMout, y(:), kb, FMbaroverout
!    real(dl) :: wlltb_tbar
    type(funcscontainer) :: radfuncs
    
    integer :: i, imax
    real(dl) :: testr, omegal, omegam, omegak, rhowout
    real(dl), parameter :: aBB = norms%aBB, a0 = norms%a0
    real(dl) :: FMbaroutoverr3, Fwbaroutoverr3,Fwbaroverout,FMbar,Fwbar
    real(dl) :: FMbaroverr, Fwbaroverr,  FMbaroverr3, Fwbaroverr3
    real(dl) :: Eofr, Eprime, wrat, Fwbaroveroutprime, Fwoutoverr2, Fwprimeoverr2, rhow, kofr, dkdr
    real(dl) :: FMoutoverr2,FMbaroveroutprime, FMprimeoverr2, Eoverr2prime, hofa_wterm
    real(dl) :: Rdot, adot, adot2
    logical, parameter :: trouble = .false.
    logical :: inierror
    type(outsidecosmology) :: thisCosmo
    
!    logical :: trouble = .true.
    if(size(y)/=index%max)then
      call wlltb_error("Vector Y has wrong size in wlltb_initime.")
    end if
    
    abar = norms%abar
    if(lfeedback>3)write(*,*)'wlltb abar:',abar,norms%abar,1/norms%abar
    inierror = .false.
   
    imax = 10
    testr = 1.d0/H0_inf
    i = 0

      omegal = Lambda /( 3._dl * H0_inf**2 )
      omegam = 8.d0/3.d0*pi*Mtilde**2/H0_inf**2
      omegak = 1. - omegam - omegal

      if(wlltbCheckPositiveDefiniteness(norms%abar,norms%a0,MyHofa2,MydHofa2da))then
        ! reject model.
        inierror=.true.
        ! just some safe values for which myhofa won't crash:
        omegak=0.d0
        omegam=0.3d0
        omegal=0.7d0
        if(lfeedback>0)write(0,*)'Sick parameters: a=1 is in an unreachable branch of an already collapsed universe for this omegak.' 
      end if
      
! fuck me blind...
!    do while (abs(radfuncs%dkdr(testr))>1.d-30)
!      i=i+1

      ! write(*,*)testr, i, abs(radfuncs%dkdr(testr))
!      testr = testr * 10
!      if(i > imax)then
!        call wlltb_error("Cannot find zero derivative of dkdr. (wlltb_tbar)")
!      end if
!    end do

!    kb = radfuncs%kofr(testr)
    kb = omegak*0.5*H0_inf**2/Mtilde**2
    kofr = radfuncs%kofr(r)
    dkdr = radfuncs%dkdr(r)

!write(0,*)'testr, kb:',testr, kb
!write(0,*)'r, kofr:',r, kofr
!write(0,*)'dkdr:',dkdr
!write(0,*)'Mtilde:',Mtilde

!    omegal =   Lambda /( 3._dl * H0_inf**2 )
!    omegak = 2.d0 * Mtilde**2 / H0_inf**2 * kb
!    omegam = 1.d0 - omegak - omegal

    if(firstcall)thisCosmo=outsidecosmology(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      
    thisCosmo%omegak = omegak
    thisCosmo%omegam = omegam
    thisCosmo%omegal = omegal
    thisCosmo%H0 = H0_inf
    
    if(firstcall)prevCosmo = thisCosmo
    
    hofa_wterm = (3*(1+w))
    Hbarout = MyHofa(abar)

    if(.not. prevCosmo == thisCosmo .or. firstcall)then
    
      tbar = wlltb_rombint(Mydtda,aBB,abar,1.d-8)
      if(isnan(tbar))then 
        inierror=.true.
        write(*,*)aBB, abar, Mydtda(aBB), Mydtda(abar)
      end if
      t0 = wlltb_rombint(Mydtda,aBB,a0,1.d-8)
      if(isnan(t0))inierror=.true.

    else

      t0 = prevCosmo%t0
      tbar = prevCosmo%tbar
      
    end if

    thisCosmo%t0 = t0
    thisCosmo%tbar = tbar
    prevCosmo = thisCosmo
    if(firstcall)firstcall=.false.
    
!write(*,*)'t0,tbar:',t0,tbar
    if(trouble)write(*,*)'tbar,t0:',tbar,t0    
    rhoMout = omegam * (a0 / abar)**3 * H0_inf**2 * 3.d0 / 8.d0 / pi
    rhowout = omegal * (a0 / abar)**(3*(1+w)) * H0_inf**2 * 3.d0 / 8.d0 / pi
    if(trouble)write(*,*)'rhom,rhowout',rhoMout, rhowout    
    FMbaroutoverr3 = 4.d0/3.d0 * pi * abar**3 * rhoMout
    Fwbaroutoverr3 = 4.d0/3.d0 * pi * abar**3 * rhowout
    if(trouble)write(*,*)'FMbaroutoverr3,Fwbaroutoverr3',FMbaroutoverr3,Fwbaroutoverr3

    Eofr = Mtilde**2 * r**2 * kofr
    Eprime = 2 * Mtilde**2 * r * kofr + Mtilde**2 * r**2 * dkdr
    Eoverr2prime = Mtilde**2 * dkdr
    wrat = ( 1.d0 + w) / (1.d0 - 3.d0 * w)
      
!    FMbaroverout = -6.d0/5.d0 / abar**2 / Hbarout**2 / r**2  * Eofr + 1.d0
!    Fwbaroverout = (-6.d0/5.d0 / abar**2 / Hbarout**2 / r**2) * Eofr * wrat + 1.d0
    FMbaroverout = (-6.d0/5.d0 / abar**2 / Hbarout**2)  * Mtilde**2 * kofr + 1.d0
    Fwbaroverout = (-6.d0/5.d0 / abar**2 / Hbarout**2 ) * Mtilde**2 * kofr * wrat + 1.d0
    if(trouble)write(*,*)'FMbaroverout,Fwbaroverout',FMbaroverout,Fwbaroverout
    if(trouble)write(*,*)'abar,Hbarout:',abar,Hbarout
    if(trouble)write(*,*)'kofr:',kofr
     
    FMbar = FMbaroverout * FMbaroutoverr3 * r**3
    Fwbar = Fwbaroverout * Fwbaroutoverr3 * r**3
    FMbaroverr = FMbaroverout * FMbaroutoverr3 * r**2
    Fwbaroverr = Fwbaroverout * Fwbaroutoverr3 * r**2
    FMbaroverr3 = FMbaroverout * FMbaroutoverr3
    Fwbaroverr3 = Fwbaroverout * Fwbaroutoverr3

    FMbaroveroutprime = 6.d0 / (5.d0 * abar**2 * Hbarout**2  ) * ( - Eoverr2prime)
    FMoutoverr2 = 4.d0 * pi / 3.d0 * abar**3 * r * rhoMout
    FMprimeoverr2 = FMbaroveroutprime * FMoutoverr2 + 4.d0 * pi * abar**3 * rhoMout * FMbaroverout
    
    rhoMbar = FMprimeoverr2 / 4.d0 / pi / abar**3
    if(trouble)write(*,*)'rhombar:',rhoMbar

    Fwbaroveroutprime = wrat * Fmbaroveroutprime
    Fwoutoverr2 = 4.d0 * pi / 3.d0 * abar**3 * r * rhowout
    Fwprimeoverr2 = Fwbaroveroutprime * Fwoutoverr2 + 4.d0 * pi * abar**3 * rhowout * Fwbaroverout
    
    rhow = Fwprimeoverr2 / 4.d0 / pi / abar**3
    
    adot2 =  2.d0 * (FMbaroverr3 + Fwbaroverr3) / abar + 2.d0 * Mtilde**2 * kofr 
    if(adot2>0.d0)then
      adot = (adot2)**0.5d0
    else
      adot = 1.d30
    end if
    
    y(index%aLTB) = abar 
    y(index%aout) = abar 
    y(index%adot) = adot
    y(index%raprime) = 0.d0
    y(index%For3) = FMbaroverr3 + Fwbaroverr3
    y(index%rhow) = rhow
    
    if(inierror)then
      y(index%aLTB) = 0.d0
      y(index%aout) = 0.d0 
      y(index%adot) = 1.d30
      y(index%raprime) = 0.d0
      y(index%For3) = 1.d30
      y(index%rhow) = 1.d30
      t0 = 0.d0
      tbar = 0.d0
      abar = 0.d0
      Hbarout = 1.d30
      rhombar = 1.d30
      kb = 1.d30
      FMbaroverout = 1.d30
    end if
    
!    if(any(isnan(y)))trouble=.true.
    if(trouble)write(*,*)y
    if(trouble)write(*,*)adot, FMbar, Fwbar, rhow
    if(trouble)stop 'trouble'
    if(any(isnan(y)))then
      write(0,*)y
      call wlltb_error('isnan(y) in wlltb_initime')
      write(0,*)'isnan(y) in wlltb_initime'
    end if
    
  contains
  
    function Mydtda(a)
      implicit none
      real(dl) :: Mydtda
      real(dl), intent(in) :: a
    
      Mydtda = 1.d0 / (a * MyHofa(a))
    
    end function Mydtda

    function MyHofa(a)
      implicit none
      real(dl) :: MyHofa
      real(dl), intent(in) :: a
    
!      MyHofa = H0_inf * sqrt((omegam * (a0/a)**3 + omegak * (a0/a)**2 + omegal * (a0/a)**(3*(1+w))))
      MyHofa = sqrt(MyHofa2(a))
      if(isnan(MyHofa))then
        write(*,*)omegam * (a0/a)**3 + omegak * (a0/a)**2 + omegal * (a0/a)**(3*(1+w))
        write(*,*)omegam, omegak, omegal, a, a0, w, kb
        call wlltb_error('isnan MyHofa')
      end if
    end function MyHofa

    function MyHofa2(a)
      implicit none
      real(dl) :: MyHofa2
      real(dl), intent(in) :: a
    
!      MyHofa = H0_inf * sqrt((omegam * (a0/a)**3 + omegak * (a0/a)**2 + omegal * (a0/a)**(3*(1+w))))
      MyHofa2 = H0_inf**2 * ((omegam * (a0/a)**3 + omegak * (a0/a)**2 + omegal * (a0/a)**hofa_wterm))
!      if(MyHofa2<0.d0)then
!        write(*,*)omegam * (a0/a)**3 + omegak * (a0/a)**2 + omegal * (a0/a)**(3*(1+w))
!        write(*,*)omegam, omegak, omegal, a, a0, w, kb
!        call wlltb_error('MyHofa2 < 0.d0')
!      end if
    end function MyHofa2

    function MydHofa2da(a)
      implicit none
      real(dl) :: MydHofa2da
      real(dl), intent(in) :: a
    
!      MyHofa = H0_inf * sqrt((omegam * (a0/a)**3 + omegak * (a0/a)**2 + omegal * (a0/a)**(3*(1+w))))
      MydHofa2da = H0_inf**2 * ((-3/a0*omegam * (a0/a)**4 -2/a0* omegak * (a0/a)**3 - hofa_wterm/a0*omegal * (a0/a)**(hofa_wterm+1)))
      if(isnan(MydHofa2da))then
        write(*,*)omegam * (a0/a)**3 + omegak * (a0/a)**2 + omegal * (a0/a)**(3*(1+w))
        write(*,*)omegam, omegak, omegal, a, a0, w, kb
        call wlltb_error('isnan MydHofa2da')
      end if
    end function MydHofa2da
    
  end subroutine wlltb_initime
  
  
  function allposdefsarepos(y)
    implicit none
    real(dl) :: y (index%max)
    logical :: allposdefsarepos
    
    allposdefsarepos = (y(index%altb)>0.d0) .and. (y(index%for3)>0.d0) .and. (y(index%rhow)>0.d0)
    
  end function allposdefsarepos

  function outsidecosmologycompare(a,b)
    implicit none
    logical :: outsidecosmologycompare
    type(outsidecosmology), intent(in) :: a, b
    
    outsidecosmologycompare = myeq(a%omegam,b%omegam) .and. myeq(a%omegak,b%omegak) .and. myeq(a%omegal,b%omegal) .and. myeq(a%H0,b%H0)

  contains
  
    function myeq(x,y)
      implicit none
      real(dl), intent(in) :: x,y
      logical :: myeq
      
      if(x==0.d0 .or. y==0.d0)then
        myeq=abs(x-y)<1.d-14
      else      
        myeq = abs(x/y-1.d0)<1.d-14
      end if
!      if(y==0.d0 .and. abs(x)<1.d-30)myeq=.true.
      
    end function myeq
    
  end function outsidecosmologycompare

 ! wlltbCheckPositiveDefiniteness checks that func(x) is positive definite on range
 ! xmin < x < xmax, returns false (!) if it is positive definite. So, true if there is an error.
 ! Function dfunc(x) is f'(x), is optional, speeds up process.
! after a coarse probing, it finds the exact minimum, to be sure we didn't miss it. 
 function wlltbCheckPositiveDefiniteness(xmin,xmax,func,dfunc)
   interface
      function afunc(x)
        use wlltb_constants
        implicit none
        real(dl) :: afunc
        real(dl), intent(in) :: x
      end function afunc
   end interface
   procedure(afunc) :: func
   procedure(afunc), optional :: dfunc
   real(dl), intent(in) :: xmin,xmax ! range over which to probe func
   logical :: wlltbCheckPositiveDefiniteness
   real(dl) :: r, Eofr, minE, minr
   real(dl) :: rt, rb, dkt, dkb
   real(dl) :: fmem(3,2) ! steps times ( x, f(x) )
   integer :: i, imax

   imax = 100!0!1000
  
   do i = 1, imax
      r = xmin + dble(i-1)/dble(imax-1)*(xmax-xmin)
!      Eofr = 1._dl + 2._dl * r**2 * VoidKofr(r) * VP%Mtilde**2
      Eofr = func(r)
      wlltbCheckPositiveDefiniteness = 0._dl >= Eofr
      if(i==1)then
        minE = Eofr
        minr = r
      end if
      if(Eofr<minE)then
        minE=Eofr
        minr = r
      end if

!        write(123,'(50ES14.5)')r/VP%L,Eofr,myf()
      if(wlltbCheckPositiveDefiniteness)return
   enddo
   
   ! find minimum:
   if(minr==xmin)return
   if(minr==xmax)return
   
    rt = minr + 1.d0/dble(imax-1)*(xmax-xmin)
    rb = minr - 1.d0/dble(imax-1)*(xmax-xmin)
    if(rb<xmin)rb=xmin

    if(.not.present(dfunc))then
      fmem(1,:)=(/rb,func(rb)/)
      r=0.5*(rt+rb)
      fmem(2,:)=(/r,func(r)/)
      fmem(3,:)=(/rt,func(rt)/)
    else
      fmem=0.d0
    end if

    r=rt
    dkt=myf()
    r=rb
    dkb=myf()
    r=0.5*(rt+rb)
    
    
    do i = 1, imax
        if(myf()*sign(1.d0,dkt)>0)then
          rt=r
          fmem(3,:)=fmem(2,:)
        else
          rb=r
          fmem(1,:)=fmem(2,:)
        end if
        r=0.5d0*(rb+rt)
        Eofr = func(r)
        fmem(2,:)=(/r,Eofr/)
        wlltbCheckPositiveDefiniteness = 0._dl >= Eofr
        if(wlltbCheckPositiveDefiniteness)return
!        write(123,'(50ES14.5)')r/VP%L,Eofr,myf()
        if(abs(myf())<1.d-13)return
     enddo
!stop 'see 123'
    return
    
  contains
  
  function myf()
    implicit none
    real(dl) :: myf
    
    if(present(dfunc))then
      myf=dfunc(r)
    else
      if(abs((fmem(3,1)-fmem(1,1)))>1.d-16)then
        myf = (fmem(3,2)-fmem(1,2))/(fmem(3,1)-fmem(1,1))
      else
        if (abs(fmem(3,2)-fmem(1,2))<1.d-14)then
          myf = 0.d0
        else
          myf = 1.d30
        end if
      end if
    end if
    
  end function myf
 end function wlltbCheckPositiveDefiniteness

end module wlltb_funcs
