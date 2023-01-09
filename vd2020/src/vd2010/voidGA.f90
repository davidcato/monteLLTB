
  ! COnversions:
  ! h == H0 / 100, H0 == H_{0, out) -> defines age.
  ! rho_crit = Mtilde^2 / 3  = 1.9 h^2 10^(-29) g cm^(-3)
  ! 1 Mpc = 3 10^{24) cm
  ! M_o = solar mass = 2 10^(33) g
  !
  ! -> rho_crit = 8.6 h^2 10^(10) M_o / Mpc^3 
  !             = Mtilde^2 
  !
  ! Gauge choice (Mtilde == Mtilde_paper):
  ! L(r) = 'gravitating mass' == 4 pi Mtilde^2 r^3  / 3
  ! rho(r) == Mtilde^2 r^2  /  (R'(r,t) R^2(r,t))
  ! M'(r) == 4 pi Mtilde^2 r^2 / sqrt( 1 + 2 r^2 Mtilde^2 k(r))
  ! M(r) = 4 pi \int_0^r dr S(r,t) R^2(r,t) rho(r,t)
  !
  ! Dim[M(r)] (Mtilde /= Mpaper):
  ! M(r) = xx * Mtilde^2 * Dim[r]^3
  ! r(mpc) = r(code) * VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde * S(r>L,t)/S(robs,t)
  ! however, rho(M0/mpc^3) is in terms of H0out, -> Mpc are in terms of outside observer.
  ! rho(mpc) = rho(code) * ( VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde )^3
  ! M(r(Mpc)) = M(r) * ( VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde )^3
  ! M(M_o) = M(r(Mpc)) / (2 Mtilde^2) * 8.6 h^2 10^10
  !        = M(r)  * ( VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde )^3 
  !               / (2 Mtilde_paper^2) * 8.6 h^2 10^10
  !        = M(r) * Vpc^3 * h^2 / VP%Mtilde_paper**5 * VP%Mtilde**3 * 4.3 * 10^1
  function VoidMassInsideRadius(rmin,rmax)
    implicit none
    real(dl) :: VoidMassInsideRadius
    real(dl), intent(in) :: rmin,rmax
!    real(dl), external :: rombint
    REAL(DL) :: testm, myrob, barem
    real(dl) :: h!, Sobs,Sout, Srat3


    h = VP%H0 / 100.d0

!!$    Sobs = avoidS(VP%r0,VP%t0)
!!$    Sout = avoidS(VP%L*1.1_dl,VP%t0)
!!$    Srat3 = (Sout / Sobs)**3

!    testm=rombint(MyIntegrand,rmin,rmax,VoidDistPrec)
    barem=voidrombint(rmin,rmax,VoidDistPrec) ! in units of Mtilde.
!    testm=testm * Vpc**3 * h**2 / VP%Mtilde_paper**5 * VP%Mtilde**3 * 43._dl
    
    ! M_sun = 1.116 10^57 GeV (Kolb and Turner)
    !       = 1.69548 10^81 sec^-1
    !       = 1.774 10^-39 Mpc^-1
    !       = 1.989 10^33 g
    !
    ! G_N = 6.672 10^-8 cm^3 g^-1 s^-2
    !
    ! barem is in units of (Mtilde and H0 [km/s/Mpc]) * r^3
    ! barem * c[km/s]^3 * (Mpc/cm)^3 is in units (km/s/Mpc)^2 * cm^3
    ! barem * c[km/s]^3 * (Mpc/cm)^3 * (km/Mpc)^2 is in units (s)^-2 * cm^3
    ! barem * c[km/s]^3 * (Mpc/cm)^3 * (km/Mpc)^2 / G_N is in units g
    ! = barem * 1.24627 10^58
    !
    ! barem * 6.26579 10^24 = barem [M_sun]
    ! 
    !      

    testm = barem * 6.226579d24 * (VP%Mtilde_paper/VP%Mtilde)**2

    VoidMassInsideRadius = testm

  contains 

    function f(r)
      real(dl), intent(in) :: r
      real(dl) :: f
      
      f = 4._dl * pi * VP%Mtilde**2  * r**2 / (1._dl + 2._dl * r**2 * VoidKofr(r) * VP%Mtilde**2)**(0.5_dl)
      
    end function f
	
    include 'vd2010/voidrombint.f90'
	
  end function VoidMassInsideRadius

  function VoidGravMassInsideRadius(rmin,rmax)
    real(dl) :: VoidGravMassInsideRadius
    real(dl), intent(in) :: rmin,rmax
    REAL(DL) :: testm, myrob
    real(dl) :: h


    h = VP%H0_paper / 100.d0

    testm=8._dl * pi * VP%Mtilde**2 * (rmax**3 - rmin**3) / 3._dl
    testm=testm * Vpc**3 * h**2 / VP%Mtilde_paper**5 * VP%Mtilde**3 * 43._dl

    VoidGravMassInsideRadius = testm

  end function VoidGravMassInsideRadius

  ! VoidGARadius returns the coordinate radius at which
  ! the density transits from over- to under-dense or
  ! vice versa, at time tin.
  function VoidGARadius(tin)
    real(dl) :: VoidGARadius
    real(dl), intent(in), optional :: tin
    integer :: i, imax
    real(dl) :: r, t, rmin,rmax, dir, a_FLRW3, bestresult(2) ! r,myfunc
    integer :: inisign, thissign, counter
    logical :: yahoo
real(dl) :: dummy, tturn
    ! reset myfunc:
    t = 0.9*VP%t0
    dummy = MyFunc(0.d0)

    if(present(tin))then
       t=tin
    else
       t=VP%t0
    end if

    a_FLRW3 = voida(VP%L*1.1d0,t)**3

    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=0._dl,t=1.d-10*VP%t0bg,tturn=tturn)
    if(t>2*tturn)then
      VoidGARadius = 1.d30
      VOidError = .true.
      return
    end if

    ! over- or under-dense? inisign.
    inisign = nint(sign(1._dl,MyFunc(0.d0)))

! Provoke error:
!if(VoidTestingIntegrationOutput)write(*,*)dummy

    imax = 100
    rmin=0._dl
    yahoo = .false.
    do i = 1, imax
       r = dble(i-1)/dble(imax-1) * VP%L

       rmax = r
       !       write(*,*) MyFunc(r)
       if(MyFunc(r)==0._dl)then
          yahoo = .true.
       end if
       thissign = nint(sign(1._dl,MyFunc(r)))

!       write(*,'(20ES24.15)')r/VP%L,MyFunc(r),  aVoidRp(r,t) , voida(r,t)**2, voidkofr(r)
!       write(123,'(20ES24.15)')r/VP%L,MyFunc(r),  aVoidRp(r,t) , voida(r,t)**2, voidkofr(r)


       if(thissign * inisign == -1)then ! Sign change w.r.t. inisign: crossed the transition.
!          write(*,*)'Found flip at r/L =',r/VP%L
          exit
       end if
       rmin = r
    end do

    if(.not.yahoo)then
       dir = -dble(inisign)
       counter = 0
       bestresult = (/r,abs(MyFunc(r))/)
       ! Added minimal value 1.d-6, otherwise we sometimes
       ! demand precision higher then double from MyFunc.
!       do while (abs(MyFunc(r)/max(r,1.d-6)) > VoidDistPrec)
       do while (abs(MyFunc(r)) > VoidDistPrec)
       
       if(abs(rmin/rmax-1)<1.d-14) exit

          if(dir*MyFunc(r) > 0.d0)then
if(abs(rmax/r-1.d0)<1.d-15)r=r*1.1
             rmax = r
          else
if(abs(rmin/r-1.d0)<1.d-15)r=r*0.9
             rmin=r
          end if
          r = 0.5_dl*(rmax+rmin)
          
          if(abs(MyFunc(r))<bestresult(2)) bestresult = (/r,abs(MyFunc(r))/)

            

          counter = counter + 1
          if(counter > 1000)then
            ! This is only a parameter inversion: accuracy need not be more than 1:10000.
            if (abs(MyFunc(r)) < sqrt( VoidDistPrec ) ) exit
            if(bestresult(2)<1.d0)then
              r=bestresult(1)
              exit
            end if
             call VoidErrorReport()
             stop 'VoidGARadius hanging.'
          end if
!       write(*,'(20ES24.15)')r/VP%L,MyFunc(r),  aVoidRp(r,t) , voida(r,t)**2, voidkofr(r),VP%L, abs(MyFunc(r)/max(r,1.d-6)) !, VP%L2_1, abs(MyFunc(r))
!       write(*,'(20ES24.15)')r,rmax,rmin,  aVoidRp(r,t) , voida(r,t)**2, voidkofr(r),VP%L, (MyFunc(r)/max(r,1.d-6)) , MyFunc(r)
!       write(123,'(20ES24.15)')r/VP%L,MyFunc(r)/r,  aVoidRp(r,t) , voida(r,t)**2, voidkofr(r),VP%L, abs(MyFunc(r)/max(r,1.d-6)) 
       end do
    end if
	if(VoidTesting)then
		write(*,*)'Precision:',MyFunc(r)
		write(*,*)'r/L at which delta density changes sign:',r/VP%L
	end	if
    VoidGARadius = r

      if(VoidGARadius<0.d0)then
        writE(0,*)'VoidGARadius<0'
        call VoidErrorReport()
        stop
      end if

!    stop 'GAradius'
!    write(*,*)'Exit GARadius:',r, VP%L

  contains
    
    function MyFunc(r)
      real(dl) :: MyFunc, r
      real(dL), save :: MyFuncMEM=0.d0,rMEM=0.d0, tMEM=0.d0
      
      if(t/=tMEM .or. r/=rMEM)then
        MyFunc = (aVoidRp(r,t) * voida(r,t)**2 - a_FLRW3)/(aVoidRp(r,t) * voida(r,t)**2 + a_FLRW3)
        MyFuncMEM = MyFunc
        rMEM = r
        tMEM = t
      else
        MyFunc = MyFuncMEM
      end if
    end function MyFunc

  end function VoidGARadius


  function VoidLogGAMassM0(rflip)
    real(dl) :: VoidLogGAMassM0
    real(dl), optional :: rflip
    real(dl) :: thisrflip

    if(.not.VP%is_pure_flrw)then
       if(present(rflip))then
          thisrflip = rflip
       else
          thisrflip = VoidGARadius()
       end if
       
       VoidLogGAMassM0 = log10( VoidMassInsideRadius(0.d0,thisrflip) )
    else
       VoidLogGAMassM0 = 0.d0
    end if

  end function VoidLogGAMassM0

  function VoidComovGARadiusMPC(rflip)
    real(dl) :: VoidComovGARadiusMPC
    real(dl), optional :: rflip
    real(dl) :: thisrflip

    if(present(rflip))then
       thisrflip = rflip
    else
       thisrflip = VoidGARadius()
    end if

    VoidComovGARadiusMPC = VoidComovDistMPC(0._dl,thisrflip,VP%t0,VP%r0)

  end function VoidComovGARadiusMPC

  subroutine VoidTestGA()
    real(dl) :: thismass
    real(dl) :: kmem, lmem, lmpc, rflip
    character(len=64) :: myformat

    myformat = '(A32,ES24.15)'
    Write(*,*)'units are M_sun and Msun/Mpc^3.'

    lmem = VP%L
    VP%L=0.d0

    lmpc = VoidComovDistMPC(0._dl,lmem,VP%t0,VP%L*1.1_dl)
    thismass = VoidMassInsideRadius(0._dl,lmem)
    write(*,myformat)'Mass (flrw):',thismass
    write(*,myformat)'Rho derived (flrw):',thismass / lmpc**3 / pi * 0.75_dl
    write(*,myformat)'If omega_k=0 Should be:', 8.6d10 * (VP%H0_paper / 100.d0)**2
    thismass = VoidGravMassInsideRadius(0._dl,lmem)
    write(*,myformat)'Gravitating mass (flrw):',thismass
    VP%L = lmem
    thismass = VoidMassInsideRadius(0._dl,lmem)
    write(*,myformat)'Mass (ltb):',thismass
    thismass = VoidGravMassInsideRadius(0._dl,lmem)
    write(*,myformat)'Gravitating mass (ltb):',thismass
    

    rflip = VoidGARadius()
    thismass = VoidMassInsideRadius(0._dl,rflip)
    write(*,myformat)'Mass (GA,ltb):',thismass
    thismass = VoidGravMassInsideRadius(0._dl,rflip)
    write(*,myformat)'Gravitating mass (GA,ltb):',thismass

    VP%L = 0.d0
    thismass = VoidMassInsideRadius(0._dl,rflip)
    write(*,myformat)'Mass (GA,flrw):',thismass
    thismass = VoidGravMassInsideRadius(0._dl,rflip)
    write(*,myformat)'Gravitating mass (GA,flrw):',thismass
    VP%L = lmem
    write(*,myformat)'Robs / r_GA:',VP%r0 / rflip

    lmpc =  VoidComovDistMPC(0._dl,rflip,VP%t0,VP%r0)
    write(*,myformat)'Size of GA (comov MPC for obs):',lmpc

    write(*,myformat)'Function VoidLogGAMassM0()',VoidLogGAMassM0()
    write(*,myformat)'Function VoidComovGARadiusMPC()',VoidComovGARadiusMPC()
    stop
    
  end subroutine VoidTestGA

  
