subroutine VoidGetColdSpot()
  use nrodeint
!  use camb
  implicit none
  real(dl) :: t0mem
  integer :: ii

  real(dl) :: DeltaTGeodesic
  real(dl) :: DeltaPhi, NewtGauge, rhox, rhobar, ax, abar, Rpx,Rpbar, deltarho, deltarho0
  real(dl) :: thisda, thisz, Lmpc, Lmem

!  real(dl) :: aV_em, aV_L
!  real(dl) :: aF_em, aF_L, tL, rL

  character(len=1024) :: outputline 
  t0mem = VP%t0

  ! A little bit extra accuracy in the geodesics (slower):
  setlboost = 1.d-3

!  write(*,'(A)')"There must be something in my definition of the times at which to comare the geodesics leaving the cold spot, that makes it identical to adiabatic photon-dust perturbations. I'm thinking of just putting the Bardeen value in the paper and forgetting about the numerics. But if we can understand why the numerics give Phi/3 ........ would be nice."

  vp%t1100=nt(z=1089.d0)
!  call TestColdSpot()

!  open(185,file="coldspot_of_z.txt")
  open(185,file="coldspot_of_L.txt")
  Lmem=VP%L
  ! loop over patch size:
  do ii = 1,1!  50 ! gives only z=1089 if i = 1,1
!     thisz = 1100.d0*(1._dl-dble(ii)/100._dl)
     thisz=1089.d0
     VP%L = dble(51-ii)/50.d0*Lmem
!     vp%t1100=nt_FLRW(z=1100.d0*(1._dl-dble(ii)/100._dl))
     vp%t1100=nt(z=thisz)
     
     thisda = DA(z=thisz)

!     DeltaPhi = VoidColdPoisson(VP%t1100)
     NewtGauge = VoidColdNewtGauge(VP%t1100,r=0.d0)
!     DeltaTGeodesic = VoidColdGeodesic(VP%t1100,r=0.d0)
     

     ! Now DeltaTGeodesic is the ordinary + integrated Sachs-Wolfe effect,
     ! assuming velocities are zero, i.e. 2 Phi + int dot 2 Phi 

!     VP%coldspot = DeltaTGeodesic - 2._DL * DeltaPhi + DeltaPhi/3._dl
!     =
!     VP%coldspot =  1._DL/3._dl * DeltaPhi
     VP%coldspot =  1._DL/3._dl * NewtGauge
     VP%coldspotmk = VP%coldspot * 2.726d6

     call VoidDistFuncs(r=0.d0,t=VP%t1100,a=ax,Rp=Rpx)
     call VoidDistFuncs(r=1.1d0*VP%L,t=VP%t1100,a=abar,Rp=Rpbar)
     rhox = 1._dl / Rpx / ax**2
     rhobar = 1._dl / Rpbar / abar**2     
     deltarho = (rhox-rhobar)/rhobar

     call VoidDistFuncs(r=0.d0,t=VP%t0,a=ax,Rp=Rpx)
     call VoidDistFuncs(r=1.1d0*VP%L,t=VP%t0,a=abar,Rp=Rpbar)
     rhox = 1._dl / Rpx / ax**2
     rhobar = 1._dl / Rpbar / abar**2     
     deltarho0 = (rhox-rhobar)/rhobar

     if(ii==1)then
        write(*,'(A,ES14.5)')'#        Coldspot (dT/T):',VP%coldspot
        write(*,'(A,F14.5)' )'# Coldspot (dT in \mu K):',VP%coldspotmk
        write(*,'(A,ES14.5)' )'#    Underdensity in LSS:',deltarho
!        Lmpc = VP%L * (VPc/1000.d0) / VP%Mtilde_paper * VP%Mtilde / voida(VP%L*1.1,vp%t0) * voida(0.d0,vp%t0)
        Lmpc = VoidComovDistMPC(VP%r0,VP%L,VP%t0)
        write(*,'(A,F14.5)' )'#Size of LTB-patch [Mpc]:',Lmpc
        write(*,'(A,F14.5)' )'#    Comov DA(LSS) [Mpc]:',thisda*(1+thisz)
        write(*,'(A,F14.5)' )'# Angular size of object (degrees):',2*Lmpc / (thisda*(1+thisz)) / pi * 180
!        write(outputline,'(3(A,ES14.2),2(A,F14.2),A)')' $',VP%delta0,'}$ & $',deltarho0,'}$ & $',deltarho,'}$ & ',2*Lmpc / (thisda*(1+thisz)) / pi * 180,' & ',VP%coldspot*2.726*1.d6,' \\'
        write(outputline,'(2(A,ES14.2),2(A,F14.2),A)')' $',VP%delta0,'}$ & $',deltarho,'}$ & ',VP%coldspot*2.726*1.d6,' & ',2*Lmpc / (thisda*(1+thisz)) / pi * 180,' \\'

        call MyReplaceChars(outputline,'E','\times 10^{')
        call MyReplaceChars(outputline,'{-0','{-')
        call MyReplaceChars(outputline,'.00\times','\times')
        write(*,'(A)')trim(outputline)
     end if



!     write(outputline,'(15ES24.15)')1100.d0*(1._dl-dble(ii)/100._dl),VP%coldspot!,DeltaPhi/3._Dl, DeltaTGeodesic, NewtGauge/3._dl!, VP%Mtilde, VP%Mtilde**2, voida(r=VP%L*1.1,t=VP%t1100)
     Lmpc = VoidComovDistMPC(0.d0,VP%L,VP%t0,VP%L)
     write(outputline,'(15ES24.15)')thisz,Lmpc,VP%coldspot!,DeltaPhi/3._Dl, DeltaTGeodesic, NewtGauge/3._dl!, VP%Mtilde, VP%Mtilde**2, voida(r=VP%L*1.1,t=VP%t1100)

     write(*,'(A)')trim(outputline)
     write(185,'(A)')trim(outputline)

     VP%t0 = t0mem
     VP%tdir=1._dl
     call Init_FLRW

  end do
  close(185)
!  stop 'Done?'
 

end subroutine VoidGetColdSpot

function VoidColdPoisson(t1100,r)
  real(dl), intent(in) :: t1100
  real(dl), intent(in), optional :: r
  real(dl) :: thispot, VoidColdPoisson, FullPot
  real(dl), save :: tsave, FullPotSave

  if(present(r))then
     thispot = VoidColdPoisson_r(t1100,r)
  else
     thispot = 0.d0
  end if
  
  if(t1100==tsave)then
     FullPot = FullPotSave
  else
     FullPot = VoidColdPoisson_r(t1100,VP%L)
  end if
  
  VoidColdPoisson = FullPot -thispot 

  FullPotSave = FullPot
  tsave = t1100

end function VoidColdPoisson

function VoidColdPoisson_r(t1100,r)
  use nrodeint
  use ode_path
  implicit none
  real(dl), intent(in) :: t1100
  real(dl), intent(in), optional :: r
  real(dl) :: VoidColdPoisson_r
  real(dl) :: vars(2),dvars(2),  intvarstart, intvarexactlist(1), stopval, nreps, dummy
  integer :: stopvar
  real(dl) :: bump,rflrw
  ! Integrate from r to L, knowing that Phi(L)=0.

  bump = 1.d-3

  vars=0._dl
!  intvarstart=0._dl
  stopvar=0
  stopval=VP%L
  if(present(r))stopval=r
  ! Translate this r to rflrw:
  dummy = VoidGrad2Phi(r=stopval,t=t1100,rflrw=rflrw)
!  stopval = rflrw
  intvarexactlist=stopval
  nreps = 1.d-4!VoidDistPrec
  nrsuccess=.true.

  ! Bump first step, because if y=dy=0 nrodeint can't start:
 
  intvarstart = bump*VP%L
  call VoidPhiDerivs(intvarstart,vars,dvars)
  vars = vars + bump*VP%L*dvars

  call odeint(vars,intvarstart, intvarexactlist,stopvar,stopval,nreps,nreps,0.d0,VoidPhiDerivs,rkqs)
  if(.not.nrsuccess)write(*,*)'Poisson is wrong.'
  VoidColdPoisson_R = - vars(1) / VP%L
end function VoidColdPoisson_r


function VoidGrad2Phi(r,t,rflrw)
  implicit none
  real(dl) :: VoidGrad2Phi
  real(dl), intent(in) :: r
  real(dl), intent(in), optional :: t
  real(dl), intent(out), optional :: rflrw
  real(dl), save :: thist, Hbar, abar, RLTBbars,Rpbars,  rhobar
  real(dl) :: RLTBx,Rpx, mydelta, ax, drflrwdr, dr, rhox, RLTBbar,Rpbar
  real(dl) :: thisr, thisrtop, thisrbot
  integer :: counter

  if(present(t))then
     thist=t
     VoidGrad2Phi=0._dl
     call VoidDistFuncs(obs='FLRW',r=1._dl,t=thist,H=Hbar,RLTB=RLTBbars,Rp=Rpbars,a=abar)
     ! r**2 = 1**2 = 1
     rhobar = 1._dl / Rpbars / RLTBbars**2
     return
  end if

  thisr = max(r,1.d-10)

  call VoidDistFuncs(r=max(r,1.d-10),t=thist,RLTB=RLTBx,Rp=Rpx, a=ax)
  
  RLTBbar = r*RLTBbars
  Rpbar = r*(Rpbars - abar)+abar
  if((abs(RLTBbar/RLTBx-1._dl)>1.d-10).and.(r>2.d-10))then
!  if(.false.)then
     thisr = r
     thisrtop = VP%L
     thisrbot = 1.d-10
     counter = 0
     do while (abs(RLTBbar/RLTBx-1._dl)>1.d-10)
        drflrwdr = Rpx/abar 
        if((counter > 5).and. (counter<800))then
           dr = (r - RLTBx/abar)/drflrwdr
        else
            dr=-1._dl
           dr=sqrt(dr) ! force nan
        end if

        if((RLTBx/RLTBbar-1._dl)>0.)then
           thisrtop = thisr
           if(isnan(dr))dr=0.5_dl*(thisrbot-thisr)
        else
           thisrbot = thisr
           if(isnan(dr))dr=0.5_dl*(thisrtop-thisr)
        end if
 
        
        if(thisr+dr>thisrtop)dr = 0.5_dl * (thisrtop - thisr)
        if(thisr+dr<thisrbot)dr = 0.5_dl * (thisrbot - thisr)

        thisr = thisr + dr

        call VoidDistFuncs(r=max(thisr,1.d-10),t=thist,RLTB=RLTBx,Rp=Rpx, a=ax)
        counter = counter + 1
        if(counter > 1000)exit
!!$        if(counter > 1000 .and. abs(RLTBx/abar/r-1._dl)<1.d-2 )exit
!!$        if(counter > 1000)then
!!$           write(*,*)thisr,r,dr
!!$           write(*,*)RLTBx/RLTBbar-1._dl,RLTBbar,RLTBx
!!$           stop 'counter > 1000'
!!$        end if
     end do

  end if

!  rhobar = 1._dl / Rpbar / RLTBbar**2
  rhox = thisr**2 / Rpx / RLTBx**2

  mydelta = (rhox-rhobar)*VP%Mtilde**2 * 8._dl * pi


  VoidGrad2Phi = 0.5_dl * abar**2 * mydelta

  if(present(rflrw))rflrw=thisr

end function VoidGrad2Phi

subroutine VoidPhiDerivs(r,y,dy)
  implicit none
  real(dl), intent(in) :: r
  real(dl), intent(in) :: y(:)
  real(dl), intent(out) :: dy(:)

  if(size(y)/=2)stop 'Wrong dimension in VoidPhiDerivs.'

! y(1) = r * Phi: (this one is much faster)
  dy(1) = y(2)
  dy(2) = r * VoidGrad2Phi(r)

! y(1) = Phi:
!  dy(1) = y(2) / max(r,1.d-10)**2
!  dy(2) = max(r,1.d-10)**2 * VoidGrad2Phi(r)

end subroutine VoidPhiDerivs

subroutine VoidChiDerivs(r,y,dy)
  implicit none
  real(dl), intent(in) :: r
  real(dl), intent(in) :: y(:)
  real(dl), intent(out) :: dy(:)
  real(dl) :: thisr

!  if(size(y)/=3)stop 'Wrong dimension in VoidEDerivs.'
  if(size(y)/=7)stop 'Wrong dimension in VoidEDerivs.'

  ! y(1) = Chi' / r
  ! y(2) = Chi_tau
  ! y(3) = Chi_tau' / r
  ! y(4) = Chi_tautau
  ! y(5) = Chi_tautau' / r
  ! dy = y'

  thisr = max(r,1.d-10)

  dy(1) = VoidE(thisr) / thisr

  dy(2) = thisr * y(3)
  dy(3) = VoidE_tau(thisr) / thisr
!write(*,'(9ES15.4)')thisr/VP%l,y,dy

  dy(4) = thisr * y(5)
  dy(5) = VoidE_tautau(thisr) / thisr
!write(*,'(11ES14.5)')r/VP%L,y,dy

  dy(6) = thisr * voidkofr(thisr)
  dy(7) = voiddkdr(thisr)

end subroutine VoidChiDerivs

function VoidE(r,t)
  use wlltb
  implicit none
  real(dl) :: VoidE
  real(dl), intent(in) :: r
  real(dl), intent(in), optional :: t
  real(dl), save :: abar, thist
  real(dl) :: Sx, Rx, thisr, ax, apx
  real(dl) :: thiscurvterm

  if(present(t))then
     thist=t
     VoidE=0._dl
     call VoidDistFuncs(obs='FLRW',r=1._dl,t=thist,a=abar)
     return
  end if

  thisr = max(r,1.d-10)

  call VoidDistFuncs(r=thisr,t=thist,RLTB=Rx,S=Sx, a = ax)
  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=thisr,t=thist,ap=apx)
  
  ! New curvature:
!  Sx = Sx * (1.d0 + 2.d0 * r**2 * VP%kb_om * VP%Mtilde**2)**0.5_dl
  thiscurvterm=(1.d0 + 2.d0 * r**2 * VP%kb_om * VP%Mtilde**2)

  VoidE = (Sx/abar)**2*thiscurvterm - (ax/abar)**2
!!$  write(*,*)'E:'
!!$  write(*,*)'1:',VoidE
!!$  voidE = (ax/abar)**2*(-2._dl * voidkofr(thisr) * thisr**2 * VP%Mtilde**2 + 2._dl * thisr * apx/ax)
!!$ write(*,*)'2:',VoidE,2._dl * voidkofr(thisr) * thisr**2 * VP%Mtilde**2 ,2._dl * thisr * apx/ax

end function VoidE

function VoidE_tau(r,t)
  use wlltb
  implicit none
  real(dl) :: VoidE_tau, VoidE_t
  real(dl), intent(in) :: r
  real(dl), intent(in), optional :: t
  real(dl), save :: abar, thist, Hbar
  real(dl) :: Sx, Rx, Sdx, Rdx, thisr, ax, Hx, adx, apdx, apx
  real(dl) :: thiscurvterm

  if(present(t))then
     thist=t
     VoidE_tau=0._dl
     call VoidDistFuncs(obs='FLRW',r=1._dl,t=thist,a=abar, H=Hbar)
     return
  end if

  thisr = max(r,1.d-10)

  call VoidDistFuncs(r=thisr,t=thist,RLTB=Rx,S=Sx,H=Hx,a=ax,Sd=Sdx)
  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=thisr,t=thist,ap=apx, apd=apdx)

  ! New curvature:
  ! Sx = Sx * (1.d0 + 2.d0 * r**2 * VP%kb_om * VP%Mtilde**2)**0.5_dl
  ! Sdx = Sdx * (1.d0 + 2.d0 * r**2 * VP%kb_om * VP%Mtilde**2)**0.5_dl

  thiscurvterm=(1.d0 + 2.d0 * r**2 * VP%kb_om * VP%Mtilde**2)
  
  Rdx = r*ax*Hx
  adx = ax*Hx
!  VoidE_tau = abar* ( 2._dl*(1._dl/abar)**2 * Sx*Sdx - 2._dl*(1._dl/r/abar)**2 * Rx * Rdx)
!  VoidE_tau = abar* ( 2._dl*(1._dl/abar)**2 * Sx*Sdx - 2._dl*(1._dl/abar)**2 * ax * adx)
  VoidE_tau =  2._dl*(1._dl/abar) * Sx*Sdx*thiscurvterm - 2._dl*(1._dl/abar) * ax * adx
  VoidE_tau = VoidE_tau +  ( -2._dl*(Sx)**2*Hbar/abar*thiscurvterm  + 2._dl*(ax)**2 * Hbar/abar)
!!$  write(*,*)'E_tau:'
!!$  write(*,*)'1:',voidE_tau
!!$
!!$  voidE_t = (-2._dl * voidkofr(thisr) * thisr**2 * VP%Mtilde**2 + 2._dl * thisr * apx/ax) * &
!!$       2._dl * ( adx * ax / abar**2 - ax**2/abar**2 * Hbar) &
!!$       + ax**2 / abar**2 * 2._dl * thisr * ( apdx/ax - apx * adx / ax**2)
!!$  voidE_tau = voidE_t * abar
!!$
!!$  write(*,*)'2:',voidE_tau

end function VoidE_tau

function VoidE_tautau(r,t) 
  use wlltb
  implicit none
  real(dl) :: VoidE_tautau, VoidE_taut, VoidE_tt
  real(dl), intent(in) :: r
  real(dl), intent(in), optional :: t
  real(dl), save :: abar, thist, Hbar, addbar, adbar
  real(dl) :: Sx, Rx, Sdx, Rdx, thisr, ax, Hx, adx, apdx, apx, thisvoidE_tau, Sddx
  real(dl) :: addx, apddx, kofr, Mt2
  real(dl) :: thiscurvterm


  if(present(t))then
     thist=t
     VoidE_tautau=0._dl
     call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=VP%L*2.,t=thist,a=abar, H=Hbar, add=addbar)
     adbar = abar*Hbar
     return
  end if

  thisr = max(r,1.d-10)
  kofr = voidkofr(thisr)
  MT2 = VP%Mtilde**2

  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=thisr,t=thist,a=ax, H=Hx, add=addx, apd=apdx,apdd=apddx, ap=apx)

  adx = ax*Hx

!!$  VoidE_tautau = (-2*thisr*(2*(2*adbar**2 - abar*addbar)*ax**2*kofr*Mt2*thisr + &
!!$           abar*(3*adbar*adx*apx + &
!!$              abar*(-2*adx*apdx - addx*apx + 2*adx**2*kofr*Mt2*thisr)) + &
!!$           ax*(-4*adbar**2*apx + abar**2*(-apddx + 2*addx*kofr*Mt2*thisr) + &
!!$              abar*(2*addbar*apx + adbar*(3*apdx - 6*adx*kofr*Mt2*thisr)))))/abar**2
!!$  
!!$ 
!!$  ! Double checked:
!!$  write(*,*)'E_tautau:'
!!$  write(*,*)'1:',VoidE_tautau
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=thisr,t=thist,a=ax, H=Hx, add=addx, apd=apdx,apdd=apddx, ap=apx, Sdd=Sddx, Sd=Sdx, S=Sx)

  ! New curvature:
!  thiscurvterm=(1.d0 + 2.d0 * r**2 * VP%kb_om * VP%Mtilde**2)**0.5_dl
!  Sx = Sx * thiscurvterm
!  Sdx = Sdx * thiscurvterm
!  Sddx = Sddx * thiscurvterm

    thiscurvterm=(1.d0 + 2.d0 * r**2 * VP%kb_om * VP%Mtilde**2)



    VoidE_tautau = (-2*(abar*(abar*addx - 3*adbar*adx)*ax + (2*adbar**2 - abar*addbar)*ax**2 + &
          ( abar**2*(adx**2 - Sdx**2) + abar*(-(abar*Sddx) + 3*adbar*Sdx)*Sx + &
           (-2*adbar**2 + abar*addbar)*Sx**2) )*thiscurvterm    )/abar**2 
!!$  write(*,*)'2:',VoidE_tautau

end function VoidE_tautau

function VoidColdNewtGauge(t1100,r)
  use wlltb
  use nrodeint
  use ode_path
  implicit none
  real(dl), intent(in) :: t1100
  real(dl), optional, intent(in) :: r
  real(dl) :: VoidColdNewtGauge
!  real(dl) :: vars(3), dvars(3)
  real(dl) :: vars(7), dvars(7)
  real(dl) :: intvarstart, intvarexactlist(1), stopval, nreps, dummy
  integer :: stopvar
  
  real(dl) :: lphi, chippchip,chitauatau, Rx, Sx, abar, abardot, Hbar, ax, offset, apx
  real(dl) :: bump, Psi, Phi, smalluphi
  real(dl), parameter :: R2 = 0.05d0, g2=(9.d0*dsqrt(2.d0)/pi)**(2.d0/3.d0)

  ! set times in VoidE and VoidE_tau = \partial_tau E:
  dummy = VoidE(r=0._dl,t=t1100)
  dummy = VoidE_tau(r=0._dl,t=t1100)
  dummy = VoidE_tautau(r=0._dl,t=t1100)

  offset = 0.d0
  bump = 1.d-4 

  ! vars(1) = Chi' / r
  ! vars(2) = Chi_tau
  ! vars(3) = Chi_tau' / r
  ! vars(3)=0._dl
  ! vars(1:2)=offset
  vars=0._dl

  ! We integrate from L to 0 (or r) because
  ! we know the boundary conditions at r>L:
  ! Chi = Chi' = Chi_tau = 0.
  nrsuccess=.true.
  intvarstart=VP%L*(1._dl-bump)
  stopvar=0
  stopval=0.d0
  if(present(r))stopval = r !max(r,stopval)
  intvarexactlist=stopval
  nreps = 1.d-5!VoidDistPrec!1.d-4!
!  nreps = VoidDistPrec!1.d-4!

  ! Bump first step, because if y=dy=0 nrodeint can't start:
  call VoidChiDerivs(intvarstart,vars,dvars)
  vars = vars - bump*VP%L*dvars

  call odeint(vars,intvarstart, intvarexactlist,stopvar,stopval,nreps,nreps,0.d0,VoidChiDerivs,rkqs)
  vars(1:2) = vars(1:2) - offset

!!$  if(r==0.d0)then
!!$     write(*,*)'at r=0:'
!!$     write(*,*)'chi''/r:',vars(1)
!!$     write(*,*)'chi_tau:',vars(2)
!!$  end if


  call VoidDistFuncs(obs='FLRW',r=1._dl,t=t1100,a=abar, H=Hbar)
  call VoidDistFuncs(r=stopval,t=t1100,RLTB=Rx,S=Sx, a = ax)
  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=stopval,t=t1100,ap=apx)
 
  ! VoidE = (Sx/abar)**2 - (Rx/r/abar)**2
  ! Rx / r = ax

  ! New curvature:
  Sx = Sx * (1.d0 + 2.d0 * r**2 * VP%kb_om * VP%Mtilde**2)**0.5_dl

  lphi = 0.5_dl - (Sx/abar)**2 / 6._dl - 1._dl/3._dl*(ax/abar)**2
  ! lphi = -lphi
  ! write(*,*)'1:',lphi
!  lphi = 0.5_dl - (ax/abar)**2*0.5_dl *( 1._dl - 2._dl/3._dl * voidkofr(stopval)*stopval**2*VP%Mtilde**2 + 2._dl/3._dl * stopval * apx/ax)
  ! write(*,*)'2:',lphi
  ! Checked: def 1 == def 2.

  chippchip = (VoidE(r=stopval) + 3._dl * vars(1)) / 6._dl

  abardot = abar*Hbar
  chitauatau = 0.5_dl * vars(2) * abardot


  VoidColdNewtGauge = lphi + chippchip + chitauatau

  Phi = VoidColdNewtGauge
  ! VoidColdNewtGauge =  -1._dl *  VoidColdNewtGauge 
  if(.not.nrsuccess)write(*,*)'NewtGauge is wrong.'
!  write(*,'(3ES15.4,A,3ES15.4)')lphi,chippchip,chitauatau,' = ', VoidColdNewtGauge
  !stop

  ! testing: vars(4) = Chi_tautau
  Psi = 0.5_dl * (vars(4) + vars(2) * abardot)
!  write(*,'(2ES15.4,A,3ES15.4)')0.5*vars(4), 0.5*vars(2)*abardot,' = ', Psi
!  write(*,*)'Phi, Psi:',VoidColdNewtGauge, Psi
!  stop
  
  smalluphi=vars(6)*6.d0*R2*g2*VP%Mtilde**2*(pi/9.d0*g2)**2
!  write(*,'(5ES24.15)')abs(-Phi/Psi-1._Dl),abs(abs(smalluphi)/abs(Phi)-1._Dl),Phi,smalluphi
  ! Tested: smalluphi = Phi = Psi for small delta0

  if(nrsuccess .and. (abs(-Phi/Psi-1._dl)>1.d-2) .and. (abs(Psi)>1.d-9))write(*,'(A,3ES14.5)')'Phi and Psi differ more than 1% in VoidColdNewtGauge. Is this model in the perturbative regime?',Phi,Psi, smalluphi
  
end function VoidColdNewtGauge

subroutine testE(t1100)
  use wlltb
  implicit none
  integer :: i, imax
  real(dl) :: t1100
  real(dl) :: lphi, chippchip,chitauatau, abar, abardot, Hbar, offset
  real(dl) :: Sx, Rx, Sdx, Rdx, thisr, ax, Hx, adx,Rpx, E, Etau

  imax = 110
  open (184,file="testtmp.txt")

  do i = 1, imax
     thisr = dble(i-1)/dble(imax-1)*VP%L
  
  
     call VoidDistFuncs(obs='FLRW',r=1._dl,t=t1100,a=abar, H=Hbar)
     call VoidDistFuncs(r=thisr,t=t1100,RLTB=Rx,S=Sx,H=Hx,a=ax,Sd=Sdx,Rp=Rpx)
     E= voidE(r=thisr)
     Etau = voidE_tau(r=thisr)
     write(184,'(I5,20E24.15)')i,thisr/VP%L,thisr,abar,Hbar,Rx,Sx,Hx,ax,Sdx, Rpx, E, Etau
  end do

  close(184)
  stop
 

end subroutine testE

function VoidColdGeodesic(t1100,r,norestore)
  use wlltb
  implicit none
  real(dl) :: VoidColdGeodesic
  real(dl), intent(in) :: t1100
  real(dl), intent(in), optional :: r
  logical, intent(in), optional :: norestore
  real(dl) :: r0mem, ax_em,Rpx_em, ax_L,Rpx_L, rhox_em,rhox_L, geodt, t0mem
  real(dl) :: aF_em, aF_L,tL,rL, aV_L, aV_em, abar, tF_em, myfac

integer :: i
real(dl) :: thisr

  myfac=1._dl

  r0mem = VP%r0
  t0mem = VP%t0

  if(present(r))then
    VP%r0 = r
  else
    VP%r0 = 0.d0
  end if
  
  VP%t0 = t1100
  VP%tdir = -1._dl
  VP%LTB_FLRW_r => VP%L
  VP%deadon = VP%LTB_FLRW_r
  
  voiderror = .false.

  call init_flrw
  
  
  call dothemagic()
  
  tL = NT(r=VP%L*myfac)
  rL = NR_FLRW(t=tL)
!!$write(*,*)'tl:',tl
!!$write(*,*)'rl:',rl
!!$write(*,*)'vp%l:',vp%l
!!$write(*,*)maxval(VPNR%ztr(1,:)),minval(VPNR%ztr(1,:))
!!$write(*,*)maxval(VPNR%ztr(2,:)),minval(VPNR%ztr(2,:))
!!$write(*,*)maxval(VPNR%ztr(3,:)),minval(VPNR%ztr(3,:))
!!$
!!$do i = -200, 100
!!$   rl=dble(i-1)/dble(99.)*VP%L
!!$   tL = NT(r=rl)
!!$   write(*,'(10ES24.15)')rl/VP%L,NZ(r=rl),NZ_FLRW(r=rl),NZ(r=rl)/NZ_FLRW(r=rl)-1._dl,NZ_FLRW(r=rl),NZ(r=rl)/NZ_FLRW(t=tl)-1._dl
!!$end do
!!$
!!$stop

  ! the stretching of the photons, a = scale factor for photon:
!  aV_em = 1_dl/(1._dl+nz(r=VP%r0))
  aV_em = 1_dl/(1._dl+nz(t=t1100))  
  aV_L = 1_dl/(1._dl+nz(r=VP%L*myfac))

  ! redshifting of flrw photons:
!  aF_em = 1._dl/(1._dl+nz_FLRW(r=VP%r0))
  aF_em = 1._dl/(1._dl+nz_FLRW(t=t1100))
  aF_L = 1._dl/(1._dl+nz_FLRW(t=tL))
  
  geodt = aV_em / aV_L / (aF_em / aF_L) - 1._dl

  
  VoidColdGeodesic = geodt
  if(.not.present(norestore))then
     ! restore:
     VP%r0 = r0mem
     VP%t0 = t0mem
     VP%deadon = VP%LTB_FLRW_r
     VP%tdir=1._dl
     call Init_FLRW
     call dothemagic
  end if
  
end function VoidColdGeodesic
 
subroutine TestColdSpot()
  implicit none

  integer :: i, imax, j, jmax
  real(dl) :: r, newtg, geodt, Rx,Rpx, rhox, Mr, Sx, poisson, ax
  real(dl) :: aF_em, aF_L,tL,rL, aV_L, aV_em, abar, rFLRW, rhobar
  real(Dl) :: thisz, zmax

  open (184,file="coldspot.txt")

!  jmax = 100
  jmax=1
  zmax=1089.d0
!  zmax=50_dl
  do j = jmax, 1, -1
!  do j = 1, jmax

     write(*,'(A)', advance="no")'Test ColdSpot(r)'
     
     thisz=zmax
     if(jmax/=1)then
        thisz=thisz*dble(j-1)/dble(jmax-1)
     end if
     vp%t1100=nt(z=thisz)
     if(vp%t1100<0.d0)then
        write(*,*)'t1100 < 0: ',VP%t1100
        write(*,*)'minval(z):',minval(VPNR%ztr(1,:))
        write(*,*)'maxval(z):',maxval(VPNR%ztr(1,:))
        stop
     end if

     imax = 50
     mr=0  
     call VoidDistFuncs(obs='FLRW',r=1.e-1_dl,t=VP%t1100,a=abar)
     call VoidDistFuncs(r=VP%L*1.1,t=VP%t1100,a=abar,Rp=Rpx, S=Sx)
     
     rhobar = 1._dl / Rpx / abar**2
     
!     do i = 1, imax
     do i = -imax+2, imax
        r = dble(i-1)/dble(imax-1)*VP%L*1.2
        newtg = VoidColdNewtGauge(vp%t1100,r=abs(r))
!        poisson = VoidColdPoisson(vp%t1100,r=r)
        geodt = VoidColdGeodesic(vp%t1100,r=r,norestore=.true.)
!        geodt = VoidColdGeodesic(vp%t1100,r=r)
        !if(r/vp%l>0.64 .and. r/vp%l<0.69)call lltb_feedback(-1)
        call VoidDistFuncs(r=r,t=VP%t1100,RLTB=Rx,Rp=Rpx, S=Sx, a=ax)
        !call lltb_feedback(0)
        if(i>1)Mr = Mr + VP%L/dble(imax-1)*rhox * Sx*Rx**2
        rFLRW = Rx/abar
        !     rhox = r**2 / Rpx / Rx**2 * VP%Mtilde**2 * 8 * pi
        rhox = 1._dl / Rpx / ax**2
        
        ! Poisson routine is wrong for r/=0.
        write(184,'(10ES24.15,F5.1,5ES24.15)')r/VP%L,r,rFLRW,newtg/3.*2.726d6,(rhox-rhobar)/rhobar, geodt*2.726d6, thisz
        write(*,'(A)', advance="no")'.'
     end do
     writE(*,*)''
     write(184,*)''
  end do

  close(184)
  stop

end subroutine TestColdSpot
