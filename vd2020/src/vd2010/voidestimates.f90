

  function AlessiosEstimate()
    implicit none
    real(dl) :: AlessiosEstimate, thisint, gamma, R2, alpha
!    integer :: i, nsteps
!    integer :: nsteps
!!$    nsteps = 5000
!!$
!!$    thisint=0.d0
!!$    dr = VP%L / dble(nsteps)
!!$    do i = 0, nsteps
!!$       thisr = dble(i)/dble(nsteps)*VP%L
!!$       thisint = thisint + dr * thisr * voidkofr(thisr)
!!$    end do

    gamma = (9.d0 * dsqrt(2.d0) / pi)**(1.d0/3.d0)
    R2 = 0.05d0
    alpha = pi * gamma**2 / 9.d0

    thisint =  4.d0 / 15.d0 * VP%kmax * VP%L**2 + VP%kb_om * VP%L**2 / 2.d0

    thisint = thisint * 6.d0 * R2 * gamma**2 * alpha**2 * VP%Mtilde**2
!    thisint = voidkofr(0.d0)*alpha**2 * gamma**2 * VP%L**2*VP%mtilde**2 / 12.5d0
    AlessiosEstimate = thisint
!!$    write(*,*)thisint, voidkofr(0.d0)*alpha**2 * gamma**2 * VP%L**2 *VP%Mtilde**2 / 12.5d0
!!$stop


  end function AlessiosEstimate

  subroutine TirthosTinyNumbers(kg2t02,krb2)
    implicit none
    real(dl) :: kg2t02, krb2
    real(dl) :: gamma, R2, alpha
    real(dl), parameter :: tau0 = (1.d0 / 6.d0 / pi)**(1.d0/6.d0)
    real(dl) :: tauF, Mb, g2, r
!    integer :: i, nsteps

    gamma = (9.d0 * dsqrt(2.d0) / pi)**(1.d0/3.d0)
    g2 = gamma**2
    R2 = 0.05d0
    alpha = pi * g2 / 9.d0
    Mb = VP%Mtilde
    r = VP%L

    tauF = tau0 - alpha * Mb * r
!!$    tau1 = -alpha * ( &
!!$                       R2 * g2 * ( &
!!$                                  tauF**2 * Mb * r * Voidkofr(r) &
!!$                                  + 2.d0 * alpha * tau0 * voidk1(VP%kmax,VP%L,r) * Mb**2 &
!!$                                  ) &
!!$                       - 1.2d0 * voidk2(VP%kmax,VP%L,r)* Mb**3 &
!!$           )
!!$   
    kg2t02 = Voidkofr(0.d0) * g2 * tau0**2
    krb2 = Voidkofr(0.d0) * (Mb * r)**2

  end subroutine TirthosTinyNumbers

  function void_B21()
    implicit none
    real(dl) :: void_B21, gamma, R2, alpha
    real(dl), parameter :: tau0 = (1.d0 / 6.d0 / pi)**(1.d0/6.d0)
    real(dl) :: tauF, Mb, tau1, g2, r, LogFrac, thisL, thisrfork
!    integer :: i, nsteps
!!$    nsteps = 5000
!!$
!!$    thisint=0.d0
!!$    dr = VP%L / dble(nsteps)
!!$    do i = 0, nsteps
!!$       thisr = dble(i)/dble(nsteps)*VP%L
!!$       thisint = thisint + dr * thisr * voidkofr(thisr)
!!$    end do

!write(*,*)tau0,(VP%Mtilde_paper*VP%t0bg)**(1.d0/3.d0),(VP%Mtilde*VP%t0)**(1.d0/3.d0)
! Yes, they are all equal.

    gamma = (9.d0 * dsqrt(2.d0) / pi)**(1.d0/3.d0)
    g2 = gamma**2
    R2 = 0.05d0
    alpha = pi * g2 / 9.d0
    Mb = VP%Mtilde!_paper
    r = VP%L!_paper
    thisL = r
    thisrfork = r * Mb / VP%Mtilde


    tauF = (tau0 - alpha * Mb * r)
    tau1 = -alpha * ( &
                       R2 * g2 * ( &
                                  tauF**2 * Mb * r * Voidkofr(thisrfork) &
                                  + 2.d0 * alpha * tau0 * voidk1(VP%kmax,thisL,r) * Mb**2 &
                                  ) &
                       - 1.2d0 * voidk2(VP%kmax,thisL,r)* Mb**3 &
           )

    LogFrac = -2.d0 * tau1 / tauF + 2.d0 * alpha * R2 * g2 * tauF * Mb*r * Voidkofr(thisrfork) + 2.d0 * alpha**2 * R2 * g2 * voidk1(VP%kmax,thisL,r) * Mb**2


    void_B21 = dexp(LogFrac)
    ! = (1 + z_LTB) / (1 + z_FLRW) = T_out / T_in = 1 / (Tv / Tf)

    return

  contains
    
    function Voidk1(kmax,L,r)
      implicit none
      real(dl) :: Voidk1, kmax,L,r 
      
      Voidk1 = kmax * (r**2 / 2.d0 + L**2 * ( -(r/L)**6/3.d0 +(r/L)**10/10.d0 ) ) !+    VP%kb_om *  r**2 / 2.d0 
    end function Voidk1

    function Voidk2(kmax,L,r)
      implicit none
      real(dl) :: Voidk2, kmax,L,r 
      
      Voidk2 = kmax * (r**3 / 3.d0 + L**3 * ( -(r/L)**7*2.d0/7.d0 +(r/L)**11/11.d0 ) )  !+  VP%kb_om * r**3 / 3.d0
    end function Voidk2


  end function Void_B21

  function simpler(z)
    implicit none
    real(dl) :: simpler, thisz, thisint, z, dz
    integer :: i, nsteps
    nsteps = 10000

    thisint=0.d0
    dz = z / dble(nsteps)
    do i = 0, nsteps
       thisz = dble(i)/dble(nsteps)*z
       thisint = thisint + dz  / (VP%omk * (1.d0 + thisz)**2+(1.d0 - VP%omk) * (1.d0 + thisz)**3)**(0.5d0)
    end do
    
!    simpler = 2.d0/3.d0 * thisint / VP%t0! * VP%Mtilde
    simpler =thisint / VP%H0! * VP%Mtilde

  end function simpler

  
  function void_B12()
    implicit none
    real(dl) :: void_B12, gamma, alpha
    real(dl), parameter :: tau0 = (1.d0 / 6.d0 / pi)**(1.d0/6.d0)
!    real(dl) :: tau0 
    real(dl) :: Mb, g2, r!, tauF, tau1, LogFrac!, thisL!, taur!, myk
    integer, parameter :: R2nmaxX = 22
    integer :: R2nmax 
!    real(dl), parameter :: R2n(R2nmax)=(/ 1.d0/20.d0, -3.d0/2800.d0, 23.d0/504000.d00, -947.d0 / 3880800000.d0, 3293.d0 / 224224.d5  /)
    real(dl), parameter :: R2n(R2nmaxX)=(/ 0.05,-0.0010714285714285715,0.00004563492063492064,-2.440218511647083e-6,1.4686206650492366e-7,-9.50993765619616e-9,6.473542653639692e-10,-4.569666554548725e-11, 3.3160977010058107e-12,-2.4592808133854967e-13, 1.856126993926667e-14,-1.4213047418173125e-15,1.1016217975387359e-16,-8.62694185023713e-18,6.816125240403192e-19,-5.427168189992692e-20,4.350653109976447e-21,-3.5086670294556637e-22,2.844806405203960395510995068614405309e-23,-2.317641992709655202149241651520015778820919e-24, 1.896350199276954934359155834823826573233873e-25,-1.557734537805811766816831444807557477260300e-26/)
!    integer :: i, nsteps
!    integer :: sptaun ! semipublic n for tau 
!    real(dl) :: sptau ! semipublic tau
!    real(dl), external :: voidrombint
	integer :: integrand
	
    real(dl) :: int1, int2, otherterm

R2nmax = R2nmaxX 

    gamma = (9.d0 * dsqrt(2.d0) / pi)**(1.d0/3.d0)
    g2 = gamma**2

    alpha = pi * g2 / 9.d0
    Mb = VP%Mtilde
    r = VP%L

!    int1 = voidrombint(integrand1,0.d0, r, 1.d-12) * Mb**2 * 2.d0 * alpha**2 / tau0**2
!    int2 = voidrombint(integrand2,0.d0, r, 1.d-12) * Mb**2 * 2.d0 * alpha**2 / tau0**2

	integrand = 1
    int1 = voidrombint(0.d0, r, 1.d-12) * Mb**2 * 2.d0 * alpha**2 / tau0**2
	integrand = 2
    int2 = voidrombint(0.d0, r, 1.d-12) * Mb**2 * 2.d0 * alpha**2 / tau0**2

!    write(*,'(4ES23.14)')int1,int2, VP%kmax, VP%L*VP%Mtilde
!    write(*,*)fr(r),fr1(r),fr2(r)
!stop
    otherterm = 0.d0

    void_B12 = dexp(int1 + int2 + otherterm)

!    void_B12 = 8.d0/5.d0 * VP%kmax * VP%Mtilde**2*VP%L**2 *g2 * alpha**2 * R2n(1) - 256.d0 / 63.d0 * VP%kmax**2 * VP%Mtilde**2*VP%L**2 * alpha**2 * g2**2 * R2n(2) * tau0**2
!    void_B12 = dexp(void_B12)
  contains

	function f(r)
		real(dl) :: f, r
		
		if(integrand==1)then
			f=integrand1(r)
		else if (integrand == 2) then
			f=integrand2(r)
		else
			stop 'Integrand not set yet in void_B12'
		end if
	end function f
	

    function integrand1(ri)
      implicit none
      real(dl) :: ri, integrand1

      integrand1 = ri * fr(ri)

    end function integrand1
    
    function integrand2(ri)
      implicit none
      real(dl) :: ri, integrand2
      real(dl) :: term1, term2, x

      x = tau0**2 * g2 * voidkofr(ri)

      term1 = 3.d0 * x * fr1(ri) + 2.d0 * x**2 * fr2(ri) - fr(ri)

      term2 = 1.d0 + fr(ri) + ri * voiddkdr(ri) * g2 * tau0**2 * fr1(ri)

      integrand2 = ri * term1 * term2

    end function integrand2
    
    function fr(rf1)!, intau)
      real(dl) :: fr, rf, rf1 , myf
      real(dl) :: x, thisk, thistau
      integer :: fi
!      real(dl), optional :: intau      
!      if(present(intau))then
!         thistau = intau
!      else
         thistau = tau0
!      end if
         rf = rf1

      thisk = voidkofr(rf)

      x = thistau**2 * g2 * thisk


      myf = 0.d0      
      do fi = 1, R2nmax

         myf = myf + R2n(fi) * x**fi
      end do


      fr = myf
    end function fr

    function fr1(rf)!, intau)
      real(dl) :: fr1, rf, myf
      real(dl) :: x, thisk, thistau
      integer :: fi
!!$      real(dl), optional :: intau      
!!$      if(present(intau))then
!!$         thistau = intau
!!$      else
      thistau = tau0
!!$      end if
      
      thisk = voidkofr(rf)
      x = thistau**2 * g2 * thisk

      myf = 0.d0      
      do fi = 1, R2nmax
         myf = myf + dble(fi) * R2n(fi) * x**(fi-1)
      end do
      fr1 = myf
    end function fr1
   
   function fr2(rf)!, intau)
      real(dl) :: fr2, rf, myf
      real(dl) :: x, thisk, thistau
      integer :: fi
!!$      real(dl), optional :: intau      
!!$      if(present(intau))then
!!$         thistau = intau
!!$      else
      thistau = tau0
!!$      end if
      
      thisk = voidkofr(rf)
      x = thistau**2 * g2 * thisk

      myf = 0.d0      
      do fi = 2, R2nmax
         myf = myf + dble(fi * (fi-1)) * R2n(fi) * x**(fi-2)
      end do
      fr2 = myf
    end function fr2

	include 'vd2010/voidrombint.f90'

 

 end function void_B12
