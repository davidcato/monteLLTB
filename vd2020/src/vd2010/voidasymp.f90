

subroutine voidasymp()
  integer :: i, imax
  real(dl) :: z, zmin, zmax
  real(dl) :: t, r, hz, zz, mte, S, Sd, SF, SFd, zf, tf, rf, zf2, rf2, zf3, zf_eff
  real(dl) :: TL, TF_obs, RF_obs, ZF_obs, RLTB, RFLRW
  real(dl), target :: lr

  imax = 1000
  zmin = 0.d0 
  zmax = 3000.d0!void_zmax_fixed

  mte = VoidTVOverTF()

  
  ! Prepare data for new zmax
  if(zmax>void_zmax)void_zmax = zmax
!  SetLBoost = 1.d-2
  call init_flrw
  call dothemagic

  call VoidSetTF(RF_obs)
  TF_obs = nt_FLRW(r=RF_obs)
  ZF_obs = z_FLRW(RF_obs)



  open(126, file="test_LTB_asymp.txt")
  do i = 1, imax
     z = zmin + (zmax - zmin) / dble(imax-1) * dble(i-1)
     r = nr(z)
!  zmax = 3.d0*z
     VP%deadon = r
!  SetLBoost = 1.d-2
!  call init_flrw
!  SetLBoost = 1.d0
!  call dothemagic
!  z = nz(r)
     t = nt(r=r)
     rf = r + RF_obs
     zf = (1+z_FLRW(rf))

     tf = nt_FLRW(r=rf)
     rf2 = nr_flrw(t=tf)
!     rf3 = nr_flrw(z=zf)

     zf2 = nz_flrwt(t)
     zf3 = nz_flrwt(tf)
     zf_eff = (1.d0+zf3)/(1.d0+ZF_obs) - 1.d0

     hz = voidH(r,t) * VP%Mtilde_paper / VP%Mtilde
     zz = ( VP%H0 / hz )**(2.d0/3.d0) * (1.d0+z) - 1.d0


     call VoidDistFuncs(z=z,S=S, Sd=Sd, RLTB = RLTB)
     call VoidDistFuncs(obs='FLRW',r=r,t=t,S=SF, Sd=SFd, Rltb=RFLRW)
     call VoidDistFuncs(obs='FLRW',z=zf3,S=SF, Sd=SFd, Rltb=RFLRW)

     lr = r
     VP%LTB_FLRW_r => lr
     mte = VoidTVOverTF()    

!     write(126,'(20ES19.7)')z, 1.d0/(1+z)*(1+zf), 1.d0/(1+z)*(1+zf2), 1.d0/(1+z)*(1+zf3), t*(Sd/S - SFd/SF), mte, log(dabs(mte-1)) - log(dabs( t*(Sd/S - SFd/SF))), RLTB, RFLRW, zf, r
     write(126,'(20ES19.7)')t, tf, r, VP%L, rf, rf2, nr_FLRW(t=TL) , VP%L, (1+z)/(1+zf_eff), z, t*(Sd/S - SFd/SF), mte

!     if( dabs(t*(Sd/S - SFd/SF)) .lt. VoidDistPrec)write(*,*)'Whooptidoo!!!',mte
  end do
!  VP%LTB_FLRW_r => VP%L

  mte = VoidTVOverTF()
  write(*,*) mte
  write(*,*)'VP%omk:',VP%omk
  close(126)
  stop 'VoidAsymp is done.'

end subroutine voidasymp

subroutine VoidSetTF(RF2_obs)
  implicit none
  integer :: counter
  real(dl), optional :: RF2_obs
  real(dL) :: tf, rf_obs, tf_obs, zf_obs, t0mem, vplmem, tl!, sf
  real(dl), parameter :: tf_precision = voidminprec, epsilon = tf_precision
  !
  real(dl) :: ttop, tbot, f, df, dt, prevf, thistr(2)
  integer, parameter :: newtraphstart = 3
  real(dl), parameter :: increase=1.2d0
  integer :: i, mystat
! trying:
  real(dl) :: rcomov, fac1, thisaoft
  logical :: skipnewtraph
  real(dl) :: myLTB_FLRW_r, r0mem, DAgoal
  integer, parameter :: TFFeedback = 0
  type(vpnr_type) :: myltbVPNR



!!$  if(VP%r0/=0.d0)then
!!$     call VoidBigAlert('RObs/=0 : must skip VoidSetTF, because this needs to be rewritten for this case. Cannot do CMB!!')
!!$     return
!!$  end if
!!$

! Why did we go to 1100 even for compensated voids?
! WV new 10/02/2011
  IF(VP%LTB_FLRW_r == VP%L) then
     void_zmax=1.5_dl*VP%zb
  end IF

  ! WV new 12/02/2011: improve compatibility with small zb:
  if(abs(vp%zb)<0.1)then
     setlboost = max(abs(VP%zb),1.d-5)
  end if

  if(present(RF2_obs))RF2_obs = VP%r0 ! 0.d0

  ! Time of Last Scattering (Synchronous / proper time gauge)

  if(TFFeedback>1)write(*,*)'Getting TL.'
  TL = nt(r=VP%LTB_FLRW_r)

  if(TFFeedback>1)write(*,*)'Getting debugging DA.'
  thistr=(/TL,VP%LTB_FLRW_r/)
  DAgoal = DA_unnorm(0.d0,thistr)

  ! The trick: FLRW observer is at r=0. So for him, DA(z) = Rltb(z).
  ! Therefore, we need to put him on the end point of a geodesic
  ! that passes through r = Rltb(z)/a(z) = DA_LTB(r=L,t=tl) / a_FLRW(t=tl).
  if(VP%r0/=0.d0)then
     if(TFFeedback>1)write(*,*)'Getting DA for off-center obs.'
     thistr=(/TL,VP%LTB_FLRW_r/)
     myLTB_FLRW_r = DA_unnorm(0.d0,thistr) / voida(t=TL,r=VP%LTB_FLRW_r)
  else
     if(TFFeedback>1)write(*,*)'Getting DA for center obs.'
     myLTB_FLRW_r = VP%LTB_FLRW_r
  end if

  r0mem=VP%r0
  t0mem = VP%t0
  VP%r0 = 0.d0

! Stop the stupid numerical inversion of the integration:
if(.false.)then
  ! Difference in coordinate distance to that time, for Void- and FLRW-obs.
  ! Nope: difference in angular diameter distance to L, in case of r0/=0.
  RF_obs = nr_FLRW(t=TL) - myLTB_FLRW_r

  ! Time of FLRW obs at which he is at coord-dist r from Last Scattering.
  TF_obs = nt_FLRW(r=RF_obs)

  ! Redshift of FLRW obs to LSS.
  ZF_obs = z_FLRW(RF_obs)

  ! But this is not the answer. delta r / delta t is not time independent.
  ! Changing r today by a translation delta, does not
  ! change r_LSS by the same delta. Because r /= physical distance.


  VP%deadon = myLTB_FLRW_r


  ! Make sure that RF_obs > r=0, such that R_LSS 
  ! is in integration results.
  ! This physics is FLRW, so r0 is irrelevant.
  t0mem = VP%t0
!  if(RF_obs>=VP%r0)then
  if(RF_obs>=0.d0)then

     tbot = nt_FLRW(r=myLTB_FLRW_r)
     TF_obs = nt_FLRW(r=RF_obs)
     ZF_obs = z_FLRW(RF_obs)
     VP%t0 = TF_obs
     call init_flrw
     if(TFFeedback>1)write(*,*)'RF_obs is inside initial results, moving on.'

  else
     if(TFFeedback>1)write(*,*)'Increasing t0 to get rf inside results.'

     do i = 1, 10
        VP%t0 = t0mem
        voiderror = .false. ! Otherwise we wouldn't have been here in the first place.
        RF_obs = -1.d0
        counter = 0
!        do while (RF_obs<VP%r0)
        do while (RF_obs<0.d0)

           call VoidDistFuncs(obs='FLRW',r=0._dl,t=vp%t0,a=thisaoft)
           
           ! Move FLRW observer to the future: his distance to 
           ! R_LSS was too small.
           VP%t0 = VP%t0 * increase**(1.d0/dble(i))
           call init_flrw
           RF_obs = nr_FLRW(t=TL) - myLTB_FLRW_r
           if(vp%zb<0.1)then
              if(VP%t0/t0mem>20)then
                 setlboost=setlboost*0.5_dl
                 voiderror=.true.
!                 writE(*,*)'This safety measure in voidasymp has never been tested.'
!                 call voiderrorreport()
! seems to work just fine
                 exit
              end if
           end if
           if(thisaoft>1.d10)voiderror=.true. ! then one of next steps lltb would crash. Sick model.
           if(VP%t0>1.d20)voiderror=.true.
           if(voiderror)exit
        end do
        if(minval(VPNR_FLRW%ztr(2,:))>TL)then
           void_zmax = void_zmax * 10.d0
           if(i.gt.3)then
              voiderror = .true.
              exit
           end if
           cycle
        end if
        if(.not.voiderror)exit

     end do
     if(TFFeedback>1)write(*,*)'Finished loop.'
     if(voiderror)then
        if(voidfeedback>0 .or. TFFEedback > 1)write(*,*)'Cannot find tF: error before we found RF_obs>0.'
        VP%r0 = r0mem
        return
     end if
     TF_obs = nt_FLRW(r=RF_obs)
     ZF_obs = z_FLRW(RF_obs)
     tbot = t0mem

  end if


  ! now start Newton-Raphson loop:
  tf = TF_obs
  f = RF_obs - epsilon
  dt = 1.d0
  counter = 0
  ttop = TF_obs
!  tbot = set already

!  We want RF(TF) to go to zero.
!  RF is the radius at which an FLRW observer sits at time TF,
!  such that his observed angular diameter distance
!  to the last scattering surface == that of the void observer.
!
!  RF is can be made zero by choosing an observation time TF, 
!  for which the photon (passing throuhg r=0,t=TF) passes
!  the LSS at VP%LTB_FLRW_r at time TL.
!  So the quest is for TF.
!
!  What we are doign here is essentially finding the solution of the geodesic-deviation-equation
! with initial condition at t=TL: r_1 = VP%LTB_FLRW_r + RF, r_2 = VP%LTB_FLRW_r
! down to t=TF: r_1 = nr_FLRW(t=TF), r_2 = r0 = 0.
! If we knew the separation of the geodesics as a function of time,
! We could give TF as a function of RF: TF is the time for 
! which r_1(TF) - deviation(TF) = 0, with deviation(TL)=RF.
! Unfortunately, I don't know deviation(TF) for kb /= 0.
! For kb=0, deviation(TF) = deviation(TL) * a(TF)/a(TL), 
! with deviation expressed as a physical distance (not comoving).

  do while(dabs(f/myLTB_FLRW_r).gt.tf_precision)
     skipnewtraph = .false.
     prevf = f
     ! FLRW = r0 independent. Even worse, FLRW obs == at r=0.
     f = RF_obs - epsilon !- VP%r0
     df = (f-prevf)/dt
!!$     if((dabs(df)>1.d30).or.(counter<newtraphstart))skipnewtraph = .true.
!skipnewtraph = .true.
!!$     if(f.gt.0.d0)then
!!$        ttop = tf
!!$        dt = -dabs(f/df)
!!$!        if((tf+dt<tbot).or.(counter.lt.newtraphstart))then
!!$        if((tf+dt<tbot).or.skipnewtraph)then
!!$           dt = -0.5d0*(tf - tbot)
!!$        end if
!!$     else
!!$        tbot = tf
!!$        dt = dabs(f/df)
!!$!        if((tf+dt>ttop).or.(counter.lt.newtraphstart))then
!!$        if((tf+dt>ttop).or.skipnewtraph)then
!!$           dt = 0.5d0*(ttop - tf)
!!$        end if
!!$     end if
     ! Retry that code:
!     dt=f/df
!     if((tf+dt>ttop).or.(tf+dt<tbot).or.skipnewtraph)then
        if(f.gt.0.d0)then
           ttop = tf
           dt = -abs(0.5d0*(tf - tbot))
        else
           tbot = tf
           dt = abs(0.5d0*(ttop - tf))
        end if
!     end if
!!$        

     if((dt<0.d0).and.(abs(tf/tbot-1.d0)<abs(f/df/tf*1.d-1)).and.(f/myLTB_FLRW_r>1.d-6))then
        tbot = tbot*0.9d0
     end if
     if((dt>0.d0).and.(abs(tf/ttop-1.d0)<abs(f/df/tf*1.d-1)).and.(f/myLTB_FLRW_r<-1.d-6))then
        ttop = ttop*1.1d0
     end if


     if(ttop==tbot .or. counter > maxsteps .or. ((counter .gt. newtraphstart) .and.(isnan(dt) .or. dabs(dt) < 1.d-30 .or. dabs(dt) > 1.d30 .or. dabs(df) < 1.d-30)))then
        if(dabs(f/myLTB_FLRW_r).lt.voiddistprec**0.5d0)then
           exit
        else
           if(voidfeedback>2)then
              write(*,*) 'Voiddistances got stalled while searching t_{F, obs}.'
              write(*,*)'tbot,t,ttop,error,goal'
              write(*,'(10ES25.16)')tbot,tf,ttop,dabs(f/myLTB_FLRW_r),voiddistprec*(10**(3*dble(counter)/dble(maxsteps)))
           end if
           VoidError=.true.
           exit
        end if
     end if

     tf = tf + dt
     VP%t0 = tf
     voiderror = .false.
     call init_flrw
     RF_obs = nr_FLRW(t=TL) - myLTB_FLRW_r
     counter = counter + 1

!     call VoidDistFuncs(t=TL,r=nr_FLRW(t=TL),Rltb=zf_obs,a=tf_obs,obs="FLRW")
!     write(*,'(6ES24.15)')DAgoal,zf_obs,nr_FLRW(t=TL),myLTB_FLRW_r, VP%t0,  RF_obs/myLTB_FLRW_r
!     write(*,*)f,vp%r0,rf_obs

  end do
  if(TFFeedback>1)write(*,*)'Original tf:',tf,'old rf_obs:',RF_obs
end if 
! END Stop the stupid numerical inversion of the integration:
! Now, just integrate from r=RL,t=TL, to r=0,t=??.
  vplmem = VP%L
  VP%L = 0.d0
  call voidasymp_integratebackwards(TL,myLTB_FLRW_r, & ! input
                                              VP%tf) ! output

  ! Before we do anything:
  ! backup VPNR to myltbVPNR,
  ! to avoid doing double integration.
  !myltbVPNR = VPNR
!  myltbVPNR%ztr => VPNR%ztr
!  myltbVPNR%ddztr => VPNR%ddztr
!  myltbVPNR%ztrerr => VPNR%ztrerr
!  myltbVPNR%rtz => VPNR%rtz
!  myltbVPNR%ddrtz => VPNR%ddrtz
!  myltbVPNR%kmaxz => VPNR%kmaxz
!  myltbVPNR%kminz => VPNR%kminz
!  myltbVPNR = VPNR ! to copy other fields than pointers
  call associateVPNRToThisVPNR(VPNR,myltbVPNR)
  !write(*,*)VPNR%ztr(1,1),myltbVPNR%ztr(1,1)
!  nullify(VPNR%ztr,VPNR%ddztr,VPNR%ztrerr,VPNR%rtz,VPNR%ddrtz,VPNR%kmaxz, VPNR%kminz)
  call nullifyThisVPNR(VPNR)
  allocate(VPNR%ztr(1,1),VPNR%ddztr(1,1),VPNR%ztrerr(1,1),VPNR%rtz(1,1),VPNR%ddrtz(1,1),VPNR%kmaxz(1), VPNR%kminz(1))
  !write(*,*)VPNR%ztr(1,1),myltbVPNR%ztr(1,1)


  ! Repair zmax:
  void_zmax = void_zmax_fixed

  VP%t0 = VP%tf
  if(TFFeedback>1)write(*,*)'Doing FLRW integration for new initial time.'
  call init_flrw()
  
  RF_obs = nr_FLRW(t=TL) - myLTB_FLRW_r
!  write(*,*)'New tf:',VP%tf,'new rf_obs:',RF_obs
!  stop
  VP%L = vplmem

  RF_obs = nr_FLRW(t=TL) - myLTB_FLRW_r

  TF_obs = nt_FLRW(r=RF_obs)

  ZF_obs = z_FLRW(RF_obs)

!  VP%tf = VP%t0
  VP%t0 = t0mem
  
  VP%r0 = r0mem

  if(TFFeedback>1)write(*,*)'restoring ltb integration.'
! don't do this:
!  call dothemagic
! but this:
  !write(*,*)VPNR%ztr(1,2),myltbVPNR%ztr(1,2)
!  nullify(VPNR%ztr,VPNR%ddztr,VPNR%ztrerr,VPNR%rtz,VPNR%ddrtz,VPNR%kmaxz, VPNR%kminz)
    if(associated(VPNR%ztr))deallocate(VPNR%ztr,stat=mystat) ! mystat, such that we can just call this when %ztr had never been allocated.
    if(associated(VPNR%ddztr))deallocate(VPNR%ddztr,stat=mystat)
    if(associated(VPNR%ztrerr))deallocate(VPNR%ztrerr,stat=mystat)
    if(associated(VPNR%rtz))deallocate(VPNR%rtz,stat=mystat) ! mystat, such that we can just call this when %ztr had never been allocated.
    if(associated(VPNR%ddrtz))deallocate(VPNR%ddrtz,stat=mystat)
    if(associated(VPNR%kmaxz))deallocate(VPNR%kmaxz,stat=mystat) ! mystat, such that we can just call this when %ztr had never been allocated.
    if(associated(VPNR%kminz))deallocate(VPNR%kminz,stat=mystat)

!  deallocate(VPNR%ztr,VPNR%ddztr,VPNR%ztrerr,VPNR%rtz,VPNR%ddrtz,VPNR%kmaxz, VPNR%kminz)
  nullify(VPNR%ztr,VPNR%ddztr,VPNR%ztrerr,VPNR%rtz,VPNR%ddrtz,VPNR%kmaxz, VPNR%kminz)
  VPNR%ztr => myltbVPNR%ztr
  VPNR%ddztr => myltbVPNR%ddztr
  VPNR%ztrerr => myltbVPNR%ztrerr
  VPNR%rtz => myltbVPNR%rtz
  VPNR%ddrtz => myltbVPNR%ddrtz
  VPNR%kmaxz => myltbVPNR%kmaxz
  VPNR%kminz => myltbVPNR%kminz
  VPNR = myltbVPNR ! to copy other fields than pointers

  !write(*,*)VPNR%ztr(1,2),myltbVPNR%ztr(1,2)
  nullify(myltbVPNR%ztr,myltbVPNR%ddztr,myltbVPNR%ztrerr,myltbVPNR%rtz,myltbVPNR%ddrtz,myltbVPNR%kmaxz, myltbVPNR%kminz)
  !write(*,*)VPNR%ztr(1,2)


  if(present(RF2_obs))RF2_obs = RF_obs

  if(TFFeedback>1)write(*,*)'t0 vs tf:',VP%t0, VP%tf

  if(TFFeedback>1)write(*,*)'Voidsettf done.'

end subroutine VoidSetTF

subroutine VoidSetAsymp
  integer :: i, imax
  real(dl) :: z, zmin, zmax
!  real(dl) :: t, hz, zz, mte, S, Sd, SF, SFd, zf, tf, rf, zf2
  real(dl) :: t, S, Sd, SF, SFd, zf, rf, zf2
  real(dl), target :: r
!  real(dl), target :: lr
  real(dl), parameter :: increase=2d0
  real(dl), parameter :: realzmax = 1100.d0
    integer :: thetime(3)
    if(voidfeedback>0)then
       call system_clock(count=thetime(1),count_rate=thetime(3))
    end if

       !    VP%LTB_maxr => VP%L
       IF(VP%LTB_maxr == VP%L) then
          VP%LTB_FLRW_r => VP%L
          VP%deadon = VP%LTB_FLRW_r

          call VoidSetTF()

          return
       end IF

!!$VoidLTB_FLRW_r = nr(z=1100.d0)
!!$VP%LTB_FLRW_r => VoidLTB_FLRW_r
!!$call VoidSetTF()
!!$return

!!$       else
!!$          if(.not.voiderror)call VoidSetAsymp()
!!$       end IF

  imax = 1000
  zmin = VP%zB
  zmax = void_zmax_fixed/2.d0
  void_zmax = 2 * zmax
  call init_flrw
  call dothemagic
!  write(*,*)'Test the duration of this routine as a function of ''increase''.'
!  write(*,*)'And adjust voidtvovertf according to new insights.'
!!$  ! Prepare data for new zmax
!!$  void_zmax = zmax
!!$!  SetLBoost = 1.d-2
!!$  call init_flrw
!!$  call dothemagic
!write(*,*)'zmax:',void_zmax,maxval(VPNR_FLRW%ztr(1,:))

  z = zmin
i=1
r = nr(zmin)
  do while (z<realzmax)
     z = z * increase**i
     if(z.gt.realzmax)z=realzmax
     if(z.gt.zmax)then
        zmax = zmax * 5
        void_zmax = 2 * zmax
        call init_flrw
!write(*,*)'zmax:',void_zmax,maxval(VPNR_FLRW%ztr(1,:))
        call dothemagic
     end if
     if(voiderror)then
        if(voidfeedback>0)write(*,*)'Cannot find asymptotic FLRW. Error before reaching zmax. (1)'
        return
     end if
     r = nr(z)

     ! new 22.04.2010
     VP%LTB_FLRW_r => r
     !write(*,'(20ES20.11)')r,zmax,void_zmax, VP%L, z, VP%zB
!     write(*,'(A,20ES20.11)')'Going.', z, maxval(VPNR%ztr(1,:)), maxval(VPNR_FLRW%ztr(1,:)),r,maxval(VPNR%ztr(3,:)),maxval(VPNR_FLRW%ztr(3,:))
     call VoidSettF()
!     write(*,*)'Returned.', size(VPNR_FLRW%ztr(1,:))
     
     if(voiderror)then
        if(voidfeedback>0)write(*,*)'Cannot find asymptotic FLRW. Error before reaching zmax. (2)'
        return
     end if

!  zmax = 3.d0*z
!  VP%deadon = r
!  SetLBoost = 1.d-2
!  call init_flrw
!  SetLBoost = 1.d0
!  call dothemagic
!  z = nz(r)
     t = nt(r=r)
     rf = r
     zf = z_FLRW(rf)
!     tf = nt_FLRW(zf)
     zf2 = nz_flrwt(t)
!write(*,*)zf,zf2
!     hz = voidH(r,t) * VP%Mtilde_paper / VP%Mtilde
!     zz = ( VP%H0 / hz )**(2.d0/3.d0) * (1.d0+z) - 1.d0

     if(voiderror)then
        if(voidfeedback>0)write(*,*)'Cannot find asymptotic FLRW. Error before reaching zmax. (3)'
        return
     end if

     call VoidDistFuncs(z=z,S=S, Sd=Sd)
!     call VoidDistFuncs(obs='FLRW',r=r,t=t,S=SF, Sd=SFd)
     call VoidDistFuncs(obs='FLRW',z=zf,S=SF, Sd=SFd)

!     lr = r
!     VP%LTB_FLRW_r => lr
!     mte = VoidTVOverTF()    

!     write(126,'(20ES19.7)')z, 1.d0/(1+z)*(1+zf), 1.d0/(1+z)*(1+zf2), t*(Sd/S - SFd/SF), mte, r, t, rf, tf, zf, zf2

!     write(*,*)dabs(t*(Sd/S - SFd/SF)),  1.d0/(1+z)*(1+zf2)

     if( dabs(t*(Sd/S - SFd/SF)) .lt. VoidDistPrec)then
        VoidLTB_FLRW_r = r
        VP%LTB_FLRW_r => VoidLTB_FLRW_r 
        if(voidfeedback>0)then
           write(*,*)'LTB_FLRW_r / VP%L:', VP%LTB_FLRW_r / VP%L
           write(*,*)'epsilon:',  dabs(t*(Sd/S - SFd/SF)) 
           write(*,*)'and z(r):',z
        end if
        return
     end if

     if(voiderror)then
        if(voidfeedback>0)write(*,*)'Cannot find asymptotic FLRW. Error before reaching zmax. (4)'
    if(feedback>2)then
       call system_clock(count=thetime(2))
       write(*,*)'Done AsympVoid in',dble(thetime(2)-thetime(1))/dble(thetime(3)),'seconds.'
    end if
        return
     end if
     i=i+1
!write(*,'(A,I3,3ES15.4)')' debug info:',i, z, dabs(t*(Sd/S - SFd/SF))
  end do
  !  mte = VoidTVOverTF()
  !  write(*,*) mte
  !  write(*,*)'VP%omk:',VP%omk
  !  close(126)
  if(voidfeedback>1)write(*,*)'It looks like this void has no asymptotic FLRW. Could that be?.'
       if(voidfeedback>1)then
          write(*,*)'LTB_FLRW_r / VP%L:', r
          write(*,*)'LTB_FLRW_r / VP%L:', VP%L
          write(*,*)'and z(r):',z
       end if
  if(voidfeedback>0)write(*,*)'Taking simplistic z = 1100 for LSS, match to FLRW there.'
  VoidLTB_FLRW_r = r
  VP%LTB_FLRW_r => VoidLTB_FLRW_r 
  
  VP%deadon = VP%LTB_FLRW_r
  
  call dothemagic
    if(feedback>2)then
       call system_clock(count=thetime(2))
       write(*,*)'Done AsympVoid in',dble(thetime(2)-thetime(1))/dble(thetime(3)),'seconds.'
    end if


end subroutine VoidSetAsymp


! voidasympdipole finds the r at which an observer observes
! exactly the dipole in the CMB that is observed by WMAP.
! subroutine VoidAsympDipole()
!   integer :: i, imin,imax, counter, j
!   real(dl) :: rol, deltat, rolbot, roltop, deltatbot, deltattop, rolbtmem(2,2),SetLBoostmem
!   real(dL), parameter :: WMAPDipole = 3.346d-3, COBETemp = 2.726, Dipole_over_T = WMAPDipole / COBETemp
!   logical :: allok, success
!   real(dL) :: f, df, drol, tfmem, t0mem, r0mem
  


! !  VoidFeedback = 0
!   if (voiderror) return
!   IF(VoidFeedback>0)write(*,*)'Getting Dipole estimate.'
!   tfmem = VP%tf
!   t0mem = VP%t0
!   r0mem = VP%r0
!   if(VP%r0/=0.d0)then
!      write(*,*)'Not going to search for r_max, here is dipole at r0 in stead.'
!      rol=VP%r0/VP%L
!      call VoidAsympGetDipole(rol,deltat)
!      call VoidDipoleOutput(rol,deltat*COBETemp)


!      ! Fix damage:
! !     call VoidSetTF()
!      VP%t0 = tfmem
!      call init_flrw()
!      VP%tf = tfmem
!      VP%t0 = t0mem
!      call DoTheMagic()
!      return
!   end if
       
!   ! New, add loop with increasing accuracy for 
!   ! rare cases
! !  rolbtmem(1,:)=(/rolbot,deltatbot/)
! !  rolbtmem(2,:)=(/roltop,deltattop/)
!   success = .true.
!   SetLBoostmem = SetLBoost
  
!   do j = 1, 4
       
!     success = .true.
!   !!$  imin = 1
!   !!$  imax = 10
!   !!$  do i = imin, imax
!   !!$     rol = 0.2d0*(dble(i)-dble(imin))/(dble(imax)-dble(imin))
!   !!$     call VoidAsympGetDipole(rol,deltat)
!   !!$     write(*,*)rol,deltat
!   !!$  end do

!     rolbot = 0.d0
!     deltatbot = 0.d0

!     roltop = 0.05d0
!     if(VoidFeedback>1)write(*,'(A)', advance="no")'*'
!     call VoidAsympGetDipole(roltop,deltat)
!     deltattop = deltat

!     allok = .false.
!     do while (deltattop < Dipole_over_T)
!        if(VoidFeedback>1)write(*,'(A)', advance="no")'-'
!   !     write(*,*)'I''m here, rol:',roltop
!        rolbot = roltop
!        deltatbot = deltattop
!        roltop = roltop + 0.1d0
!        call VoidAsympGetDipole(roltop,deltattop)
!   !     write(*,'(40ES14.5)')roltop,deltattop,Dipole_over_T
!        if(deltattop > Dipole_over_T)exit
!        if(roltop >= 1.d0)then
!           allok = .true.
!           exit
!        end if
!     end do
!     if(deltattop>1.d29)voiderror=.true.
    
!     if(voiderror)then
!        counter = 0
!        do while (voiderror)
!           if(VoidFeedback>1)write(*,'(A)', advance="no")'+'
!           voiderror = .false.
!           roltop = rolbot + (roltop - rolbot)*0.9
!           call VoidAsympGetDipole(roltop,deltattop)
!       if(deltattop>1.d29)voiderror=.true.
!           if(.not.voiderror)exit
!           counter = counter + 1
!           if(counter > 1000)then
!             if(deltatbot<1.d29)then
!               rol=rolbot
!               roltop=rolbot
!               deltattop=deltatbot
!               voiderror = .false.
!               exit
!             end if
!             call VoidErrorReport()
!              stop 'VoidAsympDipole did not make it'
!           end if
!        end do
!     end if

    

!     if(allok)then
!       rol = 1.d0
!     else if (roltop==rolbot .or. deltattop < Dipole_over_T)then
!        rol = -abs(rolbot) ! To have possiblity to track which 
!                        ! models are restricted by shell crossing
!       continue
!     else
        
!            ! Start Newton-Raphson
!            counter = 0
!            f = 1.d0
!            !rol = (roltop - rolbot) / (deltattop - deltatbot) * Dipole_over_T
!            rol = 0.5*(roltop + rolbot)
!            call VoidAsympGetDipole(rol,deltat)
!            do while (dabs(f)>1.d-5)

!               if(VoidFeedback>1)write(*,'(A)', advance="no")'.'
!               call VoidAsympGetDipole(rol,deltat)
!               f = deltat - Dipole_over_T

! !              write(*,'(40ES36.17)')rol,deltat,Dipole_over_T,f

!               if(f>0.d0)then
!                  roltop = rol
!                  deltattop = deltat
!                  df = (deltat - deltatbot) / (rol - rolbot) 
!                  drol = - dabs(f / df)
!                  if(rol+drol<rolbot)drol = (rolbot - rol)/2.d0
!               else
!                  rolbot = rol
!                  deltatbot = deltat
!                  df = (deltattop - deltat) / (roltop - rol) 
!                  drol = dabs(f / df)
!                  if(rol+drol>roltop)drol = (roltop - rol)/2.d0
!               end if
!               if(deltattop<deltatbot) then ! towards crossing zero, push forward
!                ! write(*,*)'yo'
!                 rolbot = rol
!                 roltop = rol*1.05d0
!                 rol = rol*1.025d0
!               end if
!               if(rol+drol==rol)then
!                 if(abs(f/Dipole_over_T)<1.d-1)then
!                   exit
!                 else
!                   counter = 1001
!                 end if
!               end if

!               rol = rol + drol
!               counter = counter + 1
!               if(counter > 1000 .or. rol+drol==rol)then
!                 success = .false.
!                 exit
!               end if

!            end do ! while (dabs(f)>1.d-5)
           
!     end if

!      if(success) exit ! j=1,4
!      ! else:
!      SetLBoost = (10.d0)**(-4-j)
!      if(VoidTesting)write(*,*)'increasing accuracy boost to',SetLBoost,i
!   end do ! j = 1, 4

!   if(VoidFeedback>1)write(*,*)'*'
!   if(.not.success)then
!     write(*,*)
!      call VoidErrorReport()
!       write(0,*)'f,rel f:',f,f/Dipole_over_T
!       write(0,*)'rol,drol:',rol,drol
!      stop 'VoidAsympDipole is hanging.'
!   end if
!   SetLBoost = SetLBoostmem

! !  write(*,*)'Dipole rol:',rol
! !  if(VoidFeedback>1)call VoidDipoleOutput(rol,WMAPDipole)
!   if(VoidFeedback>1)call VoidDipoleOutput(rol,deltat*COBETemp)
!   VP%r3p346 = VoidComovDistMPC(0.d0,rol*VP%L,VP%t0,VP%r0)
!   if( abs(rol - 1.d0) > 1.d-3) then
!     VP%Copernican = VoidVol(VP%t0,abs(rol*VP%L)) / VoidVol(VP%t0,VP%L)
!   else
!     VP%Copernican = 1.d0
!   end if
  

!   ! Fix damage:
!   VP%r0= r0mem
!   VP%t0 = tfmem
!   call init_flrw()
!   VP%tf = tfmem
!   VP%t0 = t0mem
!   call DoTheMagic()

! !stop 'Dipole estimation done.'
! end subroutine VoidAsympDipole


subroutine VoidDipoleOutput(rol,dipole)
  real(dl) :: rol, r, Lpaper,dipole
  real(dL) :: mtildetompc = 299792.458d0 ! = c [km / s]
  character(len=32) :: thisfmt
  character(len=128) :: mystring
  character(len=2) :: strlen
  ! r = rol * L
  ! r * Mtilde / Mtilde_paper 

  Lpaper = VP%L * VP%Mtilde / VP%Mtilde_paper * mtildetompc

  r = rol * Lpaper


  thisfmt = '(A20,F10.3)'
  strlen = '30'
  write(*,'(A'//strlen//',F10.3)')'zB = ',VP%zB
!  write(*,'(A'//strlen//',F10.3)')'L [Mpc] = ',Lpaper
!  Lpaper = VP%L * VP%Mtilde / VP%Mtilde_paper * mtildetompc / voida(r=VP%LTB_FLRW_r,t=vp%t0) * voida(r=0.d0,t=vp%t0)
  Lpaper = VoidComovDistMPC(0.d0,VP%L,VP%t0,VP%r0)
  write(*,'(A'//strlen//',F10.3)')'L [Mpc] observers scale = ',Lpaper
  Lpaper = VoidComovDistMPC(0.d0,rol*VP%L,VP%t0,VP%r0)
  write(*,'(A'//strlen//',F10.3)')'r_obs max [Mpc] obs scale = ', Lpaper
  write(*,'(A'//strlen//',F10.3)')'with a dipole [mK] of ', dipole*1.d3
  write(*,'(A'//strlen//',F10.5)')'rmax / L:',rol
  
  write(mystring,'(A,I4,A,F10.3,A,F10.3,A,F10.3,A,F10.5)')'Profile ',VoidProfile,'  & ',VP%zB,' & ',Lpaper,' & ',r,'& ',rol
  
  write(*,'(A)')trim(mystring)
  open(1024,file='VoidDipole.txt', status="old", position="APPEND")
  write(1024,'(A)')trim(mystring)// ' \\'
  close(1024)

end subroutine VoidDipoleOutput


subroutine VoidAsympGetDipole(roverL, deltat, sdeltat)
  real(dl), intent(in) :: roverl
  real(Dl), intent(out) :: deltat
  real(dl), intent(out), optional :: sdeltat ! signed deltat
  real(dl) :: dt1, dt2, r0mem, t0mem
!  
  real(dl), parameter :: zLSS = 1.1d3
  real(dl) :: tLSS
  logical :: alloutside
  
  ! changed scheme to:
  tLSS = nt(z=zLSS)

  if(voiderror)then
       deltat = 1.d30
       if(present(sdeltat))sdeltat = deltat
       if(VoidFeedback>1)call VoidRedMessage("GetAsympDipole did not pass initialization")
       return
  end if
  
  call VoidAsympDipoleOnCone(roverl * VP%L, VP%t0, tLSS, zLSS, deltat,alloutside)
  if(present(sdeltat))sdeltat=deltat
  deltat = abs(deltat)
  return
  !
!!$!VoidTestingIntegrationOutput = .false.
!
!  r0mem=VP%r0
!
!  VP%r0 = roverl * VP%L
!
!  call init_flrw
!  call dothemagic()
!  if(voiderror)then
!       deltat = 1.d30
!       if(present(sdeltat))sdeltat = deltat
!       if(VoidFeedback>1)call VoidRedMessage("GetAsympDipole did not pass initialization")
!       return
!  end if
!
!!$!  call VoidSettF()
!!$!  call dothemagic()
!
!  call VoidSetAsymp
!  dt1 = VoidTVOverTF()
!  if(voiderror)then
!       deltat = 1.d30
!       if(present(sdeltat))sdeltat = deltat
!       if(VoidFeedback>1)call VoidRedMessage("GetAsympDipole did not pass VoidSetAsymp")
!       return
!  end if
!
!
!  VP%r0 = - roverl * VP%L
!
!
!  call init_flrw
!
!
!  call dothemagic()
!!$!  call VoidSettF()
!!$!VoidTestingIntegrationOutput = .true.
!!$!  call dothemagic()
!!$!stop 'check output'
!  if(voiderror)then
!       deltat = 1.d30
!       if(present(sdeltat))sdeltat = deltat
!       if(VoidFeedback>1)call VoidRedMessage("GetAsympDipole did not pass re-initialization")
!       return
!  end if
!
!  call VoidSetAsymp
!  dt2 = VoidTVOverTF()
!  if(voiderror)then
!       deltat = 1.d30
!       if(present(sdeltat))sdeltat = deltat
!       if(VoidFeedback>1)call VoidRedMessage("GetAsympDipole did not pass 2nd VoidSetAsymp")
!       return
!  end if
!
!  deltat = dabs(dt1 - dt2)
!  if(present(sdeltat))sdeltat = dt1-dt2
!  VP%r0 = r0mem
!
!!$!  write(*,*)'#', dt1, dt2
!
end subroutine VoidAsympGetDipole

 subroutine VoidAsympDipole_vs_r()
   implicit none
   integer :: i, imin, imax
   real(dl) :: f, rol, df, deltat, t0mem, tfmem
   real(dl), parameter :: COBETemp = 2.726
!! Plot dipole vs distance from center for 0 < r < L
   write(*,*)'Entering VoidAsympDipole_vs_r'

  tfmem = VP%tf
  t0mem = VP%t0

   f = VoidComovDistMPC(0.d0,VP%L,VP%t0,VP%r0)
   imin=0
   imax=20
   open(12345,file="VoidDipole_vs_r.txt",Status="REPLACE")
!   do i = -imax+2,imax
   do i = -imax+2,imax
      rol=dble(i-1)/dble(imax)
!      SetLBoost = 1.d-4
      call VoidAsympGetDipole(rol,df,deltat)
      if(voiderror)then
         voiderror = .false.
         cycle
      end if
      write(12345,'(20ES24.15)')rol,f*rol,deltat,deltat*COBETemp
      write(*,'(A)', advance="no")'.'
  end do
  writE(*,*)
  close(12345)

  ! Fix damage:
  VP%t0 = tfmem
  call init_flrw()
  VP%tf = tfmem
  VP%t0 = t0mem
  call DoTheMagic()
 end subroutine VoidAsympDipole_vs_r


  subroutine voidasymp_integratebackwards(tl,rl, & ! input
                                              tf) ! output
    use void_dvode
    implicit none
    real(dl), intent(in) :: tl, rl
    real(dl), intent(out) :: tf
    
    ! DVODE:
    integer, parameter :: neq=3, nevents=2
    integer :: itask, istats(31), iout, ierror, istate
    real(dl) :: atol(neq), rstats(22), rtol, t_dvode, tout_dvode, dumdy(neq), g0(nevents), dummy, y(neq)
    TYPE (VODE_OPTS) :: OPTIONS
    ! END DVODE
    integer :: i, mystat
    logical :: repeat
    real(dl) :: memtdir

    if(voidfeedback>1)write(*,*)'integrating backwards to find CMB obs in FLRW.'

    memtdir = VP%tdir
    VP%tdir = 1.d0 ! Because r goes backwards -> all fine.
!write(0,*)VoidDistPrec, rl, tl, VP%t0
!stop
    RTOL = VoidDistPrec
    ATOL = VoidDistPrec*1.d-2
    ITASK = 1
    ISTATE = 1
    t_dvode = rl
    tout_dvode = 0.d0

    OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, &
        RELERR=RTOL,NEVENTS=nevents) ! Nevents is size of gfun's result.
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., ABSERR=ATOL(1), &
!        RELERR=RTOL,NEVENTS=1)
    repeat = .true.
    i = 1
    
    
  ! y(1) = t
  ! y(2) = dt/du = - (z + 1)
  ! y(3) = k
    y(1) = tl
    y(2) = -1.d1 ! irrelevant
    y(3) = voidkofr(rl) +1.d3! irrelevant too.
    
    do while (repeat)
      CALL DVODE_F90(voidasymp_derivswrap,NEQ,Y,t_dvode,tout_dvode,ITASK,ISTATE,OPTIONS,G_FCN=voidasymp_GFUN)
      ! DVODE done.
      if(ISTATE==1)then ! nothing done: t=tout
        repeat=.false.
      else if(ISTATE==2)then ! success: tout reached.
        repeat=.false.
      else if(ISTATE==3)then ! semi-success: root reached, continue integration.
        repeat=.true.

        call voidasymp_GFUN(NEQ, t_dvode, Y, Nevents, G0)
        if(abs(g0(1))<RTOL)then  
          repeat = .false.
        end if
        if(i>500)then
          call VoidErrorReport()
          write(0,*)'i>500 in voidasymp integrator'
          stop
        end if
        
      else if (ISTATE<1) then
        call VoidErrorReport()
        write(0,*)"DVODE backward background Integration did not converge in Voiddistances (asymp)."
        stop
      end if

      if(.not.repeat)exit
      i=i+1
    end do

    call RELEASE_ARRAYS
    if(associated(OPTIONS%rtol))deallocate(OPTIONS%rtol,stat=mystat)
    if(associated(OPTIONS%atol))deallocate(OPTIONS%atol,stat=mystat)

    tf = y(1)
!    write(*,*)'.... were here.. is it any good?',tf
!    write(*,*)'istate',istate
!    write(*,*)'r:',t_dvode

    VP%tdir = memtdir

  end subroutine voidasymp_integratebackwards

  subroutine voidasymp_derivswrap(neq,t,y,dy)
    implicit none
    integer, intent(in) :: neq
    real(dl), intent(in) :: t, y(neq)
    real(dl), intent(out) :: dy(neq)

    call ltb_aderivsdr(t,y,dy)

  end subroutine voidasymp_derivswrap
    
  SUBROUTINE voidasymp_GFUN(NEQ, T, Y, NG, GOUT)
    implicit none
    integer, intent(in) :: neq, NG
    real(dl), intent(in) :: t, y(NEQ)
    real(dl), intent(out) :: GOUT(ng)
    real(dl) :: dy(neq)
    real(dl) :: reladot, thissign
    
    if(ng/=2)then
      call VoidErrorReport()
      write(0,*)'NG /= Nevents in gfun in voidasymp.f90'
      stop
    end if

    gout(1) = t
    
    call voidasymp_derivswrap(neq,t,y,dy)

    if(any(abs(dy)>1.d29) .or. any(isnan(dy)) .or. notposdefsingeoy(y) .or. voiderror)then
      gout(2)=0.d0
!        write(*,*)dy
!        write(*,*)t,y
      VoidError=.true.
    else 
      gout(2) = 1.d0
    end if


  end subroutine voidasymp_GFUN


! in VoidAsympDipoleOnCone we calculate the CMB
! dipole on a point on the lightcone, for kSZ etc.
! Get z(tLSS) in both directions, take difference.
!subroutine VoidAsympDipoleOnCone(z, r, t, tLSS, zLSS, deltat,alloutside) ! why was z here?
subroutine VoidAsympDipoleOnCone(r, t, tLSS, zLSS, deltat,alloutside)
  real(dl), intent(in) :: r, t, tLSS, zLSS !  z,
  real(Dl), intent(out) :: deltat
  logical, intent(out) :: alloutside
  real(dl) :: dt1, dt2, r0mem, t0mem, zgoal, zmaxmem, zLSS1, zLSS2
  real(dl) :: rLSS1, rLSS2 ! for testing
  real(dl), parameter :: maxTempDiff = 1.d0
  real(dl) :: SetLBoostmem
    
  alloutside = .false.
  setlboostmem = SetLBoost
  SetLBoost = min(SetLBoost,1.d-4)
!  r0mem = VP%r0
!  t0mem = VP%t0

!  VP%r0 = r
!  VP%t0 = t
  
  ! Rescale, because dothemagic always starts with z=0. No problem.
!  zgoal = min((1.d0+zLSS)/(1.d0+z)*(1.d0+maxTempDiff) - 1.d0,zLSS) ! maxTempDiff is only a safety margin!
!  zmaxmem = void_zmax
!  void_zmax = zgoal

!  call dothemagic()
!  if(voiderror)then
!       deltat = 1.d30
!       if(VoidFeedback>1)call VoidRedMessage("VoidAsympDipoleOnCone did not pass initialization")
!       return
!  end if

!  zLSS1 = nz(t=tLSS)
!write(*,*)zLSS1,(1+zLSS)/(1+z),'==?'
!  rLSS1 = nr(zLSS1)
!  zLSS1=(1.d0+zLSS)/(1.d0+z)-1.d0
  ! duh...
  if(voiderror)then
       deltat = 1.d30
       if(VoidFeedback>1)call VoidRedMessage("VoidAsympDipoleOnCone did not pass VoidSetAsymp")
       return
  end if


!  VP%r0 = - r
  !write(*,'(A,10ES14.5)')'VoidAsympDipoleOnCone z, r, t, tLSS, zLSS:',z, r, t, tLSS, zLSS
  call voidasymp_integratekSZtoLSS(r,t,tLSS,zLSS1,rLSS1)
  call voidasymp_integratekSZtoLSS(-r,t,tLSS,zLSS2,rLSS2)
  
  SetLBoost = SetLBoostmem
  
  !call dothemagic()
  if(voiderror)then
       deltat = 1.d30
       if(VoidFeedback>1)call VoidRedMessage("VoidAsympDipoleOnCone did not pass re-initialization")
       return
  end if

!  zLSS2 = nz(t=tLSS)
!  rLSS2 = nr(zLSS2)
  if(voiderror)then
       deltat = 1.d30
       if(VoidFeedback>1)call VoidRedMessage("VoidAsympDipoleOnCone did not pass 2nd VoidSetAsymp")
       return
  end if

    dt1 = (1.d0+zLSS1)
    dt2 = (1.d0+zLSS2)
  deltat = (zLSS1 - zLSS2)/(dt1+dt2)*2

  ! Analytical zero:
  if(rLSS1>VP%LTB_FLRW_r .and. rLSS2<-VP%LTB_FLRW_r)then
  !if(rLSS2<-VP%LTB_FLRW_r)then
    deltat=0.d0
    alloutside = .true. ! for the calling routine to know when to stop
  end if

!  write(*,'(A,20ES14.5)')'# DipoleonCone ',r,t,deltat,zLSS1,rLSS1,zLSS2,rLSS2,tLSS,zLSS
end subroutine VoidAsympDipoleOnCone

subroutine VoidSetkSZZmax()
  implicit none
  type(vpnr_type) :: myVPNR!, myVPNR_FLRW
  real(dl) :: tLSS
!write(*,*)'Starting VoidSetkSZZmax'
  tLSS = VP%t1100
  ! Integrate from VP%L outward to t0:
  if(tLSS>nt(r=VP%LTB_FLRW_r))then
    VP%kSZzmax = void_zmax ! this only happens when LSS is inside LTB, which never happens unless 
    ! we have pure FLRW.
  else
    call voidasymp_integratekSZzmax(tLSS,VP%LTB_FLRW_r, & ! input
                                              VP%kSZzmax) ! output
  end if
  
!write(*,*)'Done VoidSetkSZZmax'

end subroutine VoidSetkSZZmax


! Integrate from outside starting at last scattering time, 
! towards larger r but also growing time. That is, integrate
! from LTB edge from last scattering outward to today.
! Stop when this geodesic crosses the observers
! lightcone: that defines zmax, the largest distance
! at which the LSS of a remote obsever is affected 
! by the LTB patch.
  subroutine voidasymp_integratekSZzmax(tl,rl, & ! input
                                              zmax) ! output
    use void_dvode
    implicit none
    real(dl), intent(in) :: tl, rl
    real(dl), intent(out) :: zmax
    
    ! DVODE:
    integer, parameter :: neq=3, nevents=5
    integer :: itask, istats(31), iout, ierror, istate
    real(dl) :: atol(neq), rstats(22), rtol, t_dvode, tout_dvode, dumdy(neq), g0(nevents), dummy, y(neq)
    TYPE (VODE_OPTS) :: OPTIONS
    ! END DVODE
    integer :: i, mystat
    logical :: repeat
    real(dl) :: memtdir


    memtdir = VP%tdir
    VP%tdir = -1.d0 ! Because r goes backwards -> all fine.
!write(0,*)VoidDistPrec, rl, tl, VP%t0
!stop
    RTOL = VoidDistPrec
    ATOL = VoidDistPrec*1.d-2
    ITASK = 1
    ISTATE = 1
    t_dvode = rl
    tout_dvode = 1.d4

    OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, &
        RELERR=RTOL,NEVENTS=nevents) ! Nevents is size of gfun's result.
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., ABSERR=ATOL(1), &
!        RELERR=RTOL,NEVENTS=1)
    repeat = .true.
    i = 1
    
    
  ! y(1) = t
  ! y(2) = dt/du = - (z + 1)
  ! y(3) = k
    y(1) = tl
    y(2) = -1.d1 ! irrelevant
    y(3) = voidkofr(rl) +1.d3! irrelevant too.
    call voidkSZzmax_GFUN(NEQ, t_dvode, Y, Nevents, G0)
    if(abs(g0(1))<RTOL)then  
      repeat = .false.
    end if
    if(any(g0(3:5)<0))then
      repeat = .false.
    end if

    do while (repeat)
      CALL DVODE_F90(voidasymp_derivswrap,NEQ,Y,t_dvode,tout_dvode,ITASK,ISTATE,OPTIONS,G_FCN=voidkSZzmax_GFUN)
      ! DVODE done.
      if(ISTATE==1)then ! nothing done: t=tout
        repeat=.false.
      else if(ISTATE==2)then ! success: tout reached.
        repeat=.false.
      else if(ISTATE==3)then ! semi-success: root reached, continue integration.
        repeat=.true.

        call voidkSZzmax_GFUN(NEQ, t_dvode, Y, Nevents, G0)
        if(abs(g0(1))<RTOL)then  
          repeat = .false.
        end if
        if(any(g0(3:5)<0))then
          repeat = .false.
        end if
        if(i>500)then
          call VoidErrorReport()
          write(0,*)'i>500 in voidasymp_integratekSZzmax integrator'
          stop
        end if
        
      else if (ISTATE<1) then
        call VoidErrorReport()
        write(0,*)"DVODE backward background Integration did not converge in Voiddistances voidasymp_integratekSZzmax (asymp)."
        stop
      end if

      if(.not.repeat)exit
      i=i+1
    end do

    call RELEASE_ARRAYS
    if(associated(OPTIONS%rtol))deallocate(OPTIONS%rtol,stat=mystat)
    if(associated(OPTIONS%atol))deallocate(OPTIONS%atol,stat=mystat)

    VP%tdir = memtdir

    ! Result:
    zmax = nz(r=t_dvode)

  end subroutine voidasymp_integratekSZzmax
    
  SUBROUTINE voidkSZzmax_GFUN(NEQ, T, Y, NG, GOUT)
    implicit none
    integer, intent(in) :: neq, NG
    real(dl), intent(in) :: t, y(NEQ)
    real(dl), intent(out) :: GOUT(ng)
    real(dl) :: dy(neq)
    real(dl) :: reladot, thissign
    real(dl) :: ntLightCone
    
    if(ng/=5)then
      call VoidErrorReport()
      write(0,*)'NG /= Nevents in gfun in voidkSZzmax_GFUN in voidasymp.f90'
      stop
    end if

    gout(1) = t
    
    call voidasymp_derivswrap(neq,t,y,dy)

    if(any(abs(dy)>1.d29) .or. any(isnan(dy)) .or. notposdefsingeoy(y) .or. voiderror)then
      gout(2)=0.d0
!        write(*,*)dy
!        write(*,*)t,y
      VoidError=.true.
    else 
      gout(2) = 1.d0
    end if

    ! Test when this CMB photon crosses the lightcone.
    ntLightCone = nt(r=t)
    if(ntLightCone<y(1))then
      gout(3)=-1.d0
    else
      gout(3)=1.d0
    end if

    if(y(1)>VP%t0)then
      gout(4)=-1.d0
    else
      gout(4)=1.d0
    end if

    if(ntLightCone<0)then
      gout(5)=-1.d0
    else
      gout(5)=1.d0
    end if

!write(*,'(100ES14.5)')t,y,gout,ntLightCone,VP%t0

  end subroutine voidkSZzmax_GFUN

!(-r,t,tLSS,zLSS2)
! Integration with t as time parameter, NOT r!!!
  subroutine voidasymp_integratekSZtoLSS(rini, tini, tgoal, & ! input
                                              zLSS,rLSS) ! output
    use void_dvode
    implicit none
    real(dl), intent(in) :: rini, tini, tgoal
    real(dl), intent(out) :: zLSS, rLSS
    
    ! DVODE:
    integer, parameter :: neq=3, nevents=2
    integer :: itask, istats(31), iout, ierror, istate
    real(dl) :: atol(neq), rstats(22), rtol, t_dvode, tout_dvode, dumdy(neq), g0(nevents), dummy, y(neq)
    TYPE (VODE_OPTS) :: OPTIONS
    ! END DVODE
    integer :: i, mystat
    logical :: repeat
    real(dl) :: memtdir

!    if(voidfeedback>1)write(*,*)'integrating backwards to find zmax for kSZ.'

    RTOL = VoidDistPrec * SetLBoost
    ATOL = VoidDistPrec*1.d-2 * SetLBoost
    ITASK = 1
    ISTATE = 1
    t_dvode = tini
    tout_dvode = tgoal

    OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, &
        RELERR=RTOL,NEVENTS=nevents) ! Nevents is size of gfun's result.
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., ABSERR=ATOL(1), &
!        RELERR=RTOL,NEVENTS=1)
    repeat = .true.
    i = 1
    
    
  ! y(1) = t
  ! y(2) = dt/du = - (z + 1)
  ! y(3) = k
    y(1) = rini
    y(2) = -1.d0
    y(3) = voidkofr(rini)
    
    do while (repeat)
      ! Integration with t as time parameter, NOT r!!!
      CALL DVODE_F90(voidasymp_integratekSZtoLSS_derivswrap,NEQ,Y,t_dvode,tout_dvode,ITASK,ISTATE,OPTIONS,G_FCN=voidasymp_integratekSZtoLSS_GFUN)
      ! DVODE done.
      if(ISTATE==1)then ! nothing done: t=tout
        repeat=.false.
      else if(ISTATE==2)then ! success: tout reached.
        repeat=.false.
      else if(ISTATE==3)then ! semi-success: root reached, continue integration.
        repeat=.true.
        if(VoidError)exit

        call voidasymp_integratekSZtoLSS_GFUN(NEQ, t_dvode, Y, Nevents, G0)

        if(abs(g0(1))<RTOL)then  
          repeat = .false.
        end if
        if(i>500)then
          call VoidErrorReport()
          write(0,*)'i>500 in voidasymp integrator'
          stop
        end if
        
        istate = 2
        
      else if (ISTATE<1) then
        if(VoidError)then
          exit
        else
          call VoidErrorReport()
          write(0,*)"DVODE backward background Integration did not converge in Voiddistances voidasymp_integratekSZtoLSS (asymp)."
          stop
        end if
      end if

      if(.not.repeat)exit
      i=i+1
    end do

    call RELEASE_ARRAYS
    if(associated(OPTIONS%rtol))deallocate(OPTIONS%rtol,stat=mystat)
    if(associated(OPTIONS%atol))deallocate(OPTIONS%atol,stat=mystat)

    zLSS = -y(2)-1.d0
    rLSS = y(1)

  end subroutine voidasymp_integratekSZtoLSS
    
  SUBROUTINE voidasymp_integratekSZtoLSS_GFUN(NEQ, T, Y, NG, GOUT)
    implicit none
    integer, intent(in) :: neq, NG
    real(dl), intent(in) :: t, y(NEQ)
    real(dl), intent(out) :: GOUT(ng)
    real(dl) :: dy(neq)
    real(dl) :: reladot, thissign
    real(dl) :: ntLightCone
    
    if(ng/=2)then
      call VoidErrorReport()
      write(0,*)'NG /= Nevents in gfun in voidkSZzmax_GFUN in voidasymp.f90'
      stop
    end if

    gout(1) = t
    
    call voidasymp_integratekSZtoLSS_derivswrap(neq,t,y,dy)

!    if(any(abs(dy)>1.d29) .or. any(isnan(dy)) .or. notposdefsingeoy(y) .or. voiderror)then
    if(any(abs(dy)>1.d29) .or. any(isnan(dy)) .or. voiderror)then ! notposdefsingeoy is for r-integration.
      gout(2)=0.d0
!        write(*,*)dy
!        write(*,*)t,y
      VoidError=.true.
    else 
      gout(2) = 1.d0
    end if
  end subroutine voidasymp_integratekSZtoLSS_GFUN

  ! Integration with t as time parameter, NOT r!!!
  subroutine voidasymp_integratekSZtoLSS_derivswrap(neq,t,y,dy)
    implicit none
    integer, intent(in) :: neq
    real(dl), intent(in) :: t, y(neq)
    real(dl), intent(out) :: dy(neq)
    real(dl) :: r, my(neq)
    real(dl) :: thisS
    ! y(1) = r
    ! y(2) = z
    r=abs(y(1))
    my(1)=t
    my(2:)=y(2:)
    call ltb_aderivsdr(r,my,dy)
    dy(1) = dy(1)**(-1)
    dy(2:) = dy(2:)*dy(1)

! Actually, now LLTB is working properly.
! The following problem can STILL occur, 
! and this is possible if the shell crossing
! radius is very narrow, such that it is 
! easily 'jumped over' by the original outward
! integration, and now we run into it when we 
! start the integration on a point close to the
! shell crossing: initial steps are small
! and only grow if the integration allows for it.
if(.false.)then
    if(VoidError .AND. VP%tdir*dy(1)>0.d0)then
        write(*,*)'dy:'
        write(*,*)dy(1:2),t
      ! delete cache:
      call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=0.d0,t=t,S=thisS)
      write(0,*)'enabling lltb_feedback'
      call lltb_feedback(100)
      call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,S=thisS)
        write(*,*)dy(1:2),t
        write(*,*)y(1:2),t
        write(*,*)thisS
    call VoidErrorReport()
    stop 'Problems for dipole on past lightcone. Should not happen if lightcone is healthy.'
    end if
end if
  end subroutine voidasymp_integratekSZtoLSS_derivswrap

