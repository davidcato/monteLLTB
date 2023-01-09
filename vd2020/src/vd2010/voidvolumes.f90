subroutine VoidMakeCritDens()
  implicit none
  real(dl) :: thisval

  call VoidTestVolumes()
  call sett0bg_Mt()
  
  VP%L=1._dl
  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=0._dl,t=1.d-10*VP%t0bg,tturn=vp%tturn)
  if(vp%t0bg.ge.2._dl*vp%tturn)voiderror=.true.
  VP%L=0._dl
  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=0._dl,t=1.d-10*VP%t0bg,tturn=thisval)
  if(vp%t0bg.ge.2._dl*thisval)voiderror=.true.
  if(voiderror)stop 'Times are fucked.'
  call init_FLRW()
  call VoidSEtL()
  
end subroutine VoidMakeCritDens


function DecelerationQ(r,t)
  implicit none
  real(dl) :: DecelerationQ
  real(dl) :: add, r, t, H, a

    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,add=add, a=a, H=H)
    
    DecelerationQ = - add * a / H**2 
    
end function DecelerationQ

subroutine VoidTestVolumes()
  implicit none
  real(dl) :: lmem, kmaxmem
  integer :: i, imax
  lmem=VP%L
  kmaxmem=VP%kmax

  write(*,*)'Void initialized. Starting volumes.'
  write(*,*)'VoidError:',VoidError
!  if(voiderror)then
!     stop 'Cannot perform volume routines in case of VoidError. Makes no sense.'
!  end if
!write(*,*)'L=VP%L/10'
!VP%L = 1.d-1*lmem
!VP%kmax=kmaxmem*4.d0
!  call VoidVolMakeFlat()
!write(*,*)'L=VP%L',VP%L
!VP%L=lmem
!VP%kmax=kmaxmem
  
  VoidError=.false.


!   imax = 1000
!  do i = 1, imax
!     VP%L = lmem * 10**(dble(i-50)/10.d0)
     call VoidVolMakeCurved(VP%kb_om)
!  end do

!  call VoidVolMakeFlat()
!  call VoidMAssComps
!  call VoidTimeComps
!  call VoidRadiusComps()
!  call VoidVolCurvTest
!  stop 'Finished VoidTestVolumes'
     write(*,*)'Finished VoidTestVolumes'

end subroutine VoidTestVolumes


subroutine VoidVolMakeFlat()
  implicit none
  real(dl) :: kb_ommem, ktop, kbot, dk, thism, thisl
  integer :: i
  
  kb_ommem = VP%kB_om
  ktop=2._dl*1._dl / (2._Dl * VP%Mtilde**2  / voidH(1.1*VP%L,VP%t0)**2 / voida(1.1*VP%L,VP%t0)**2 )
  kbot=-ktop
  VP%kb_om = 0.d0
  do i = 1, 1000
     if (.not.VoidCheckKofr()) then
        thism = VoidUnnormMassInsideRadius(0.d0,VP%L)
        thisl = VP%L**3 / 3._dl
     else
        thism = 1.d0
        thisl = 0.d0
     end if
     if(thism<thisl)then
        ktop = VP%kb_om
     else
        kbot = VP%kb_om
     end if
     VP%kb_om = 0.5_dl*(ktop + kbot)        
     if(dabs(thism/thisl-1.d0)<1.d-14)then
        go to 10
     end if
!write(*,*)thism,thisl,VP%kb_om
  end do
  stop 'VoidVolMakeFlat is hanging'

10 Continue

  write(*,*)'New kb:',VP%kb_om
  write(*,*)'m(r):',thism
  write(*,*)'l(r):',thisl
  write(*,*)'Omega_k:',2._Dl * VP%Mtilde**2 * VP%kb_om / voidH(1.1*VP%L,VP%t0)**2 / voida(1.1*VP%L,VP%t0)**2 
  write(*,*)'t0(Gyr):',VP%t0_gyr
end subroutine VoidVolMakeFlat

subroutine VoidVolMakeCurved(kb)
  implicit none
  real(dl), intent(in) :: kb
  real(dl) :: kb_ommem, ktop, kbot, dk, thism, thisl,kmaxmem, thisq
  integer :: i, FBmem
  
  kb_ommem = VP%kB_om
  kmaxmem=VP%kmax
  VP%kmax=0.d0
  ktop=1.d2
  kbot=-ktop
  VP%kb_om = kb
  do i = 1, 1000
     if (.not.VoidCheckKofr()) then
        thism = VoidUnnormMassInsideRadius(0.d0,VP%L)
        thisl = VP%L**3 / 3._dl
     else
        thism = 1.d0
        thisl = 0.d0
     end if
     if(thism<thisl)then
        ktop = VP%kmax
     else
        kbot = VP%kmax
     end if
     VP%kmax = 0.5_dl*(ktop + kbot)        
     if(dabs(thism/thisl-1.d0)<1.d-14)then
        go to 10
     end if
!write(*,*)thism,thisl,VP%kmax
  end do
  stop 'VoidVolMakeCurved is hanging'

10 Continue

  call sett0bg_Mt()
  thisq = DecelerationQ(r=VP%L*1.1d0,t=VP%t0)
  write(*,*)'Old kb:',kb_ommem
  write(*,*)'New kb:',VP%kb_om
  write(*,*)'Old kmax:',kmaxmem
  write(*,*)'New kmax:',VP%kmax
  write(*,*)'Ratio:',VP%kmax/kmaxmem
!write(*,*)VP%L,VP%kmax,VP%kb_om
  write(*,*)'m(r):',thism
  write(*,*)'l(r):',thisl
  write(*,*)'Omega_k:',2._Dl * VP%Mtilde**2 * VP%kb_om / voidH(1.1*VP%L,VP%t0)**2 / voida(1.1*VP%L,VP%t0)**2 
  write(*,*)'DecelerationQ:',thisq
  kb_ommem=VP%kb_om
  VP%kb_om=0.d0
  call sett0bg_Mt()
  thisq = DecelerationQ(r=VP%L*1.1d0,t=VP%t0)
  write(*,*)'Flat DecelerationQ:',thisq
  VP%kb_om=kb_ommem
  call sett0bg_Mt()
!!$  kb_ommem=VP%kb_om
!!$  VP%kb_om=0.d0
!!$  FBmem=VoidFeedback
!!$  VoidFeedback=0
!!$  call sett0bg_Mt()
  write(*,*)'t0(Gyr):',VP%t0_gyr
!!$  VP%kb_om=kb_ommem
!!$  call sett0bg_Mt()
!!$  VoidFeedback=FBmem
  write(*,*)'delta0:',VP%delta0
end subroutine VoidVolMakeCurved

subroutine VoidMassComps()
  implicit none
  real(dl) :: lmem
  real(dl) :: MLTB, MFLRW

  lmem = VP%L

  MLTB = VoidNumMassInsideRadius(0.d0,lmem)

  VP%L = 0.d0

  MFLRW = VoidNumMassInsideRadius(0.d0,lmem)

  VP%L = lmem

  write(*,*)'Mass comparison:'
  write(*,*)MLTB, MFLRW, MLTB/MFLRW - 1.d0


  
end subroutine VoidMassComps

subroutine VoidVolCurvTest()
  implicit none
  integer :: i, imax
  real(dl) :: r, int1, int2

  imax = 1000000
  int1=0
  int2=0
  do i = 1, imax
     r=dble(i-1)/dble(imax-1)*VP%L*1.1
     int1=int1+1._Dl/dble(imax-1)*VP%L*1.1*r**2
     int2=int2+1._Dl/dble(imax-1)*VP%L*1.1*r**2/(1._dl + 2._dl * r**2 * VoidKofr(r) * VP%Mtilde**2)**(0.5_dl)
     
!     write(*,'(4ES15.6)')r,int1,int2,int1/int2-1.d0
  end do
     write(*,'(4ES15.6)')r,int1,int2,int1/int2-1.d0

end subroutine VoidVolCurvTest

subroutine VoidTimeComps()
  implicit none
  integer :: i, imax
  real(dl) :: t

  write(*,*)'LTB, FLRW, rel diff:'

  imax = 10000
!  do i = imax, 1, -1
  do i = 1, imax
     t = (dble(i)/dble(imax)) * VP%t0 * 2.
!     t = nt(z=1089.d0*dble(i)/dble(imax))
     call VoidCompVols(t)
  end do
  
end subroutine VoidTimeComps

subroutine VoidCompVols(t)
  implicit none
  real(dl), intent(in) :: t
  real(dl) :: lmem, kmem, kbmem, L2_1mem, d0_2mem
  real(dl) :: volLTB, volFLRW!, surfLTB, surfFLRW

  lmem = VP%L
  kmem = VP%kmax
  L2_1mem = VP%L2_1
  d0_2mem = VP%d0_2
  kbmem = VP%kb_om

  volLTB = VoidVol(t=t,r=lmem)
!  surfLTB = 4._dl * pi * lmem**2 * voida(lmem,t)

  VP%L = 0.d0

VP%kb_om=0.d0  
  volFLRW = VoidVol(t=t,r=lmem)
!  surfFLRW = 4._dl * pi * lmem**2 * voida(lmem,t)
VP%kb_om=kbmem

  VP%L = lmem

  write(*,'(5ES15.4)')t/VP%t0, volLTB, volFLRW, volLTB/volFLRW - 1.d0, voidH(0.d0,t)/abs(voidH(0.d0,t))

end subroutine VoidCompVols


function VoidVol(t,r)
  implicit none
  real(dl), intent(in) :: t
  real(dl), intent(in), optional :: r
  real(dl) :: VoidVol
  real(dl) :: rr, dum
!  REAL(DL), external :: rombint

  if(present(r))then
     rr=r
  else
     rr=VP%L
  end if

!  VoidVol = rombint(VoidVolIntegrand, 0.d0, rr, VoidDistPrec)
  VoidVol = voidrombint( 0.d0, rr, VoidDistPrec)

contains

  function f(r)
    IMPLICIT none
    real(dl), intent(in) :: r
    real(dl) :: f
    
    f = avoidS(r,t) * voida(r,t)**2 * r**2

  end function f
  
  include 'vd2010/voidrombint.f90'
  
  
end function VoidVol



  function VoidNumMassInsideRadius(rmin,rmax)
    real(dl) :: VoidNumMassInsideRadius
    real(dl), intent(in) :: rmin,rmax
    real(dl), external :: rombint
    REAL(DL) :: testm, myrob
    real(dl) :: h!, Sobs,Sout, Srat3


    h = VP%H0 / 100.d0

!!$    Sobs = avoidS(VP%r0,VP%t0)
!!$    Sout = avoidS(VP%L*1.1_dl,VP%t0)
!!$    Srat3 = (Sout / Sobs)**3

!    testm=rombint(MyIntegrand,rmin,rmax,VoidDistPrec)
    testm=voidrombint(rmin,rmax,VoidDistPrec)
!    testm=testm * Vpc**3 * h**2 / VP%Mtilde_paper**5 * VP%Mtilde**3 * 43._dl

    VoidNumMassInsideRadius = testm

  contains 

    function f(r)
      real(dl), intent(in) :: r
      real(dl) :: f
      
      f = 8._dl * pi * VP%Mtilde**2  * r**2 / (1._dl + 2._dl * r**2 * VoidKofr(r) * VP%Mtilde**2)**(0.5_dl)
      
    end function f
	
	include 'vd2010/voidrombint.f90'
	
  end function VoidNumMassInsideRadius

  function VoidUnnormMassInsideRadius(rmin,rmax)
    real(dl) :: VoidUnnormMassInsideRadius
    real(dl), intent(in) :: rmin,rmax
    real(dl), external :: rombint
    REAL(DL) :: testm, myrob
    real(dl) :: h!, Sobs,Sout, Srat3
    real(dl) :: myprec

!    MyPrec = VoidDistPrec
    MyPrec = 1.d-12
!    testm=rombint(MyIntegrand,rmin,rmax,MyPrec)
    testm=rombint(rmin,rmax,MyPrec)

    VoidUnnormMassInsideRadius = testm

  contains 

    function f(r)
      real(dl), intent(in) :: r
      real(dl) :: f
      
      f =  r**2 / (1._dl + 2._dl * r**2 * VoidKofr(r) * VP%Mtilde**2)**(0.5_dl)
      
    end function f
	
	include 'vd2010/voidrombint.f90'

  end function VoidUnnormMassInsideRadius
  
  function VoidVolPhysRadius(t)
    real(dl) ::  VoidVolPhysRadius
    real(dl), intent(in) :: t
    real(dl), external :: rombint
    REAL(DL) :: testm, myrob
    real(dl) :: h!, Sobs,Sout, Srat3
    real(dl) :: myprec, rmin, rmax

    rmin=0.d0
    rmax=VP%L

    MyPrec = VoidDistPrec
!    MyPrec = 1.d-10
!    MyPrec=1.d-4
!    testm=rombint(MyIntegrand,rmin,rmax,MyPrec)
    testm=rombint(rmin,rmax,MyPrec)

     VoidVolPhysRadius = testm

  contains 

    function f(r)
      real(dl), intent(in) :: r
      real(dl) :: f
      
      f = avoidS(r,t)
      
    end function f

	include 'vd2010/voidrombint.f90'

  end function VoidVolPhysRadius
  
  subroutine VoidRadiusComps()
    implicit none
    real(dl) :: t, LTBr, FLRWr, kmaxmem, kbmem, L2_1, d0_2, kmax_2
    integer :: i, imax
    real(dl) :: zmax, thisa
    
    kmaxmem=VP%kmax
    kbmem=VP%kb_om
    d0_2 = VP%d0_2
    kmax_2 = VP%kmax_2
    L2_1 = VP%L2_1

    void_zmax=void_zmax_fixed
    call init_FLRW()
    zmax=maxval(VPNR_FLRW%ztr(1,:))

    open(1234,file="RadiiRatio.txt",Status="REPLACE")

    imax = 1000
    do i = 1, imax

!       t=dble(i)/dble(imax)*VP%t0

       thisa=1.d0/(1.d0+zmax*dble(i)/dble(imax))
       t=voidtofa(VP%L,thisa)
!       write(*,*)'a_LTB:',voida(VP%L,t)
       
       LTBr = VoidVolPhysRadius(t)

       VP%kb_om=0.d0
       VP%kmax=0.d0
       VP%kmax_2 = 0.d0
       VP%d0_2 = 0.d0
       VP%L2_1 = 0.d0
       
       t=voidtofa(VP%L,thisa)
!       write(*,*)'a_FLRW:',voida(VP%L,t)
       FLRWr = VoidVolPhysRadius(t)

       VP%kb_om=kbmem
       VP%kmax=kmaxmem
       VP%kmax_2 = kmax_2
       VP%d0_2 = d0_2
       VP%L2_1 = L2_1
       
       if(modulo(i,10)==1)write(*,'(A)', advance="no")'.'

       write(1234,'(5ES24.15)')t/VP%t0,voida(VP%L*1.1d0,VP%t0)/voida(VP%L*1.1d0,t)-1.d0,LTBr,FLRWr,LTBr/FLRWr-1.d0
    end do
    
    close(1234)
    write(*,*)'RadiiRatios.txt written.'

  end subroutine VoidRadiusComps
