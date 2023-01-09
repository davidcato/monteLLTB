!> DoTheMagic performs the geodesics integration, for the currently loaded cosmology.
!> The results are always stored in the global type VNPR. 
  subroutine DoTheMagic!(Vp,Vpnr)
!    use NROdeInt
!    use ode_path
    use void_dvode
    implicit none
!    real(dl) :: inivars(3), vars(3) ! 1=t, 2=z, 3=trick to not miss peak
!    real(dl), allocatable :: inivars(:), vars(:) ! 1=t, 2=z, 3=trick to not miss peak, 4=A, 5 dA/dl
!    integer, parameter :: maxkount = 100, nvar = 5, finecoarsefixed=30
    integer, parameter :: maxkount = 300, nvar = 5, finecoarsefixed=100
    integer :: finecoarse
!    integer, parameter :: maxkount = 25, nvar = 5, finecoarse=18
    real(dl) :: inivars(nvar), vars(nvar) ! 1=t, 2=z, 3=trick to not miss peak, 4=A, 5 dA/dl
    real(dl) :: xp(maxkount), yp(nvar,maxkount)
    real(dl) :: v0, v1, nreps
    integer :: stopvar!, nvar
    integer :: kount
    real(dl) :: stopval, intvarstart, intvarexactlist(1)
    ! for testing:
    integer :: i, mystat!, testimax
!    real(dl) :: dtm_tr(1:2), testr
!    real(dl), allocatable :: ddyy(:,:), testy(:)
    ! DVODE:
    integer, parameter :: neq=nvar, nevents=6
    integer :: itask, istats(31), iout, ierror, istate
    real(dl) :: atol(neq), rstats(22), rtol, t_dvode, tout_dvode, dumdy(neq), g0(nevents), dummy, y(neq)
    TYPE (VODE_OPTS) :: OPTIONS
    ! END DVODE
    real(dl) :: stopvarinival, thisstopval, ystopvarfinecoarse
    logical, parameter :: nrsuccess = .true. ! deprecated
    integer :: thispower

    v0=VP%r0 !0.d0 ! r

    if(.not. VoidError) then
        
    !    dxsav=0.d0 ! simply save all steps in integration.
    !    save_steps=.true.

        intvarstart=VP%r0  !0.d0 ! r
        finecoarse = finecoarsefixed
        if(v0>VP%L )finecoarse=0
    !    if(VP%r0 /= 0.d0 .and. .not.(Do_bestfit_output.or.VoidTesting))stop 'Off-center observer while not testing. Do not do that.'

        ! new off-centerode_path
    !    if(VP%r0/=0.d0)then
    !       nvar=5 ! PARAMETER 
    !    else
    !       nvar=3
    !    end if
    !    allocate(vars(nvar),inivars(nvar))
        if(nvar>4)then
           inivars(4)=0._dl
           inivars(5)=1._dl ! only valid of inivars(2)=-1.
        end if
        ! end off center

        v0=intvarstart !0.d0 ! r
        v1=1.d10

        inivars(1)= VP%t0
        inivars(2)= - 1.d0 !  -(z +1)
        inivars(3)= VP%kmax + VP%kb_om

        if(inivars(1)<0.)then
           write(*,*)'DoTheMagic called with negative initial t. Must be a bug.'
           call VoidErrorReport()
           stop 'DoTheMagic called with negative initial t. Must be a bug.'
        end if

        intvarexactlist(1)=VP%deadon

        vars=inivars

        stopvar = 2 !
        stopval = - (1.d0 + void_zmax) ! now including FLRW outside void. !-(1.d0+VP%zB)  ! stop at 1+z = 1+zB

        if(vp%tdir<0.)then
           stopval = -1._dl!-1.d-2!-1.d0/1100.d0
           vars(2) = -1101.d0
           inivars(2) = -1101.d0
        end if

    !    if(voidfeedback>1)write(*,'(A,F20.11,A,F20.11,A,F20.11,A,F20.11)')' Stopping beyond z =',void_zmax, ', for zB =',VP%zB, ', L =',VP%L!,', r_deadon =',VP%deadon
!        if(voidfeedback>1)write(*,'(A,F20.11,A,F20.11,A,F20.11,A,F20.11)')' Stopping beyond y(2) =',stopval, ', for zB =',VP%zB, ', L =',VP%L,' precision:',VoidDistPrec*SetLBoost!,', r_deadon =',VP%deadon
        if(voidfeedback>1)write(*,'(A,F20.11,A,F20.11,A,F20.11,A,ES11.2)')' Stopping beyond y(2) =',stopval, ', for zB =',VP%zB, ', L =',VP%L,' precision:',VoidDistPrec*SetLBoost!,', r_deadon =',VP%deadon
        

        nreps = VoidDistPrec * SetLBoost
    !    nrsuccess=.true.
        RTOL = nreps
        ATOL = RTOL*1.d-2
        ITASK = 1
        ISTATE = 1
        t_dvode = v0
        tout_dvode = 1.d4+v0
        
        OPTIONS = SET_NORMAL_OPTS(ABSERR_VECTOR=ATOL, &
            RELERR=RTOL,NEVENTS=nevents) ! Nevents is size of gfun's result.
    !    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., ABSERR=ATOL(1), &
    !        RELERR=RTOL,NEVENTS=1)

    !    call odeint(vars,intvarstart, intvarexactlist, stopvar,stopval,nreps,nreps,0.d0,ltb_nderivsdr,rkqs)
    !    call odeint(vars,intvarstart, intvarexactlist,stopvar,stopval,nreps,nreps,0.d0,ltb_aderivsdr,rkqs)
        xp=0.d0
        yp=0.d0
        xp(1) = v0
        yp(:,1)=inivars
        stopvarinival = inivars(stopvar)
        ystopvarfinecoarse = inivars(stopvar)
    !write(*,*)vars,inivars
    !write(*,*)abs(vars(1))*rtol+atol(1),abs(vars(2))*rtol+atol(2),abs(vars(3))*rtol+atol(3),abs(vars(4))*rtol+atol(4),abs(vars(5))*rtol+atol(5)
        do i = 1, maxkount
!if(VP%tdir < 0 .and. i>29)VoidFeedback=100
          kount = i
    !      thisstopval = stopvarinival + dble(i)/dble(maxkount)*(stopval - stopvarinival)
          thisstopval = stopval
          if(i==1+finecoarse)ystopvarfinecoarse = vars(stopvar)
          ! Fine gridding for r < L
          if(i<=finecoarse)then ! fine graining in k(r), coarse outside.
            if(VP%L>0.d0)then 
              tout_dvode = (VP%L-v0) * dble(i)/dble(finecoarse) + v0
            else          
              tout_dvode = 1.d0/VP%H0 * dble(i)/dble(maxkount-finecoarse)+ v0
            end if
            ! overrule for DE_fitting:
            if(VoidIncludeFLRWDEFit/=0 .and. vp%tdir>0.)then
              tout_dvode =1.d4/VP%H0
              thisstopval = - ( 1.d0 + min(void_zmax,VoidIncludeFLRWDEFitzmax) * dble(i)/dble(finecoarse))
            end if
          else ! then gridding in z, first fine, then bigger steps for higher z.
    !        tout_dvode = tout_dvode + 4.d0/VP%H0 * dble(i-finecoarse)/dble(maxkount-finecoarse)
            tout_dvode = 1.d4/VP%H0
            if(tout_dvode<t_dvode)then
                if(voidfeedback>2)call VoidRedMessage('adjusting tout_dvode')
                tout_dvode=abs(t_dvode*1.d2)
            end if
            if(vp%tdir>0 .and. finecoarse>0 )then
              thispower = 5
            else
              thispower = 1
            end if
    !        thisstopval = vars(stopvar) + (dble(i-finecoarse) / dble(maxkount - finecoarse))**5 * (stopval - vars(stopvar))

            thisstopval = ystopvarfinecoarse +  (dble(i-finecoarse) / dble(maxkount - finecoarse))**thispower * (stopval - ystopvarfinecoarse)
          end if
          ! write(*,*)v0,t_dvode,tout_dvode, VP%H0,-vars(2)-1
          ! write(*,*)-vars(2)-1, i, thisstopval, ystopvarfinecoarse,t_dvode,tout_dvode
          ! Do integration
      if(voidfeedback>2)then
        call VoidRedMessage('Before dvode')
          write(0,'(A,3ES15.6)')'t_dvode,tout_dvode, v0:',t_dvode,tout_dvode,v0
      end if
          CALL DVODE_F90(DoTheMagic_derivswrap,nvar,vars,t_dvode,tout_dvode,ITASK,ISTATE,OPTIONS,G_FCN=dothemagic_GFUN)
          ! If integration didn't work: exit
      if(voidfeedback>2)then
        if(VoidError)call VoidRedMessage('voiderror after dvode')
          write(0,*)'istate:',istate
      end if
          if(all(ISTATE/=(/1,2,3/)))voiderror=.true.
          if(voiderror)exit
          ! store results of this step
          yp(:,i)=vars
          xp(i)=t_dvode
          if(abs(vars(stopvar)/stopval-1.d0)<nreps)exit
          if(ISTATE/=1)ISTATE=2
   
        end do

    ! dvode cleanup:
    call RELEASE_ARRAYS
    if(associated(OPTIONS%rtol))deallocate(OPTIONS%rtol,stat=mystat)
    if(associated(OPTIONS%atol))deallocate(OPTIONS%atol,stat=mystat)


        if(voidfeedback>1)write(*,*)'Finished geodesic integration. Error?',VoidError
        !write(*,*)'zmax dothemagic', yp(stopvar,kount),stopval, kount, nrsuccess, voiderror
        if(.not.nrsuccess)VoidError = .true.

    !  VoidTestingIntegrationOutput = .true.
        if(VoidTestingIntegrationOutput)then
           if(voidfeedback>2)write(*,*)'Writing integration results to test_ltb.txt.'
           open(123,file='test_ltb.txt')
           do i = 1, kount
              write(123,'(20ES14.5)')xp(i),yp(:,i), avoidS(xp(i),yp(1,i)), avoidSdot(xp(i),yp(1,i)), avoidRp(xp(i),yp(1,i)), avoidS(xp(kount),yp(1,i)), avoidSdot(xp(kount),yp(1,i)), voidR(xp(i),yp(1,i))
           end do
           close(123)
           if(voidfeedback>2)write(*,*)'Closed test_ltb.txt.'
        end if
    !write(*,*)i,kount,voiderror, ISTATE
    !if(.not.voiderror)stop 'i''m here'

    !    if(VoidError)return
    ! moved to after (de)-allocation
        
    end if ! end if (.not. voiderror)

    if(VoidError)then
      inivars= 0.d0

      yp=0.d0
      xp=0.d0
      kount=3
      xp(1)=1.d0
      xp(2)=2.d0
      xp(3)=3.d0
      xp(4)=4.d0
      
      yp(:,1)=1.d0      
      yp(:,2)=2.d0
      yp(:,3)=3.d0
      yp(:,4)=4.d0
    end if

!    if(associated(VPNR%ztr))deallocate(VPNR%ztr,stat=mystat) ! mystat, such that we can just call this when %ztr had never been allocated.
!    if(associated(VPNR%ddztr))deallocate(VPNR%ddztr,stat=mystat)
!    if(associated(VPNR%ztrerr))deallocate(VPNR%ztrerr,stat=mystat)
!    if(associated(VPNR%rtz))deallocate(VPNR%rtz,stat=mystat) ! mystat, such that we can just call this when %ztr had never been allocated.
!    if(associated(VPNR%ddrtz))deallocate(VPNR%ddrtz,stat=mystat)
    call deallocateThisVPNR(VPNR)

    ! NEW: add 1, to include initial conditions in the arrays.
    allocate(VPNR%ztr(nvar,kount+1)) 
    allocate(VPNR%ddztr(nvar,kount+1))
    allocate(VPNR%ztrerr(nvar,kount+1))
    allocate(VPNR%rtz(nvar,kount+1)) 
    allocate(VPNR%ddrtz(nvar,kount+1))

! Nope, do not return yet: first go through GetLocalMinMax_in_z in order
! to allocate some pointers!
!!$    if(VoidError)then
!!$       deallocate(vars,inivars)
!!$       return
!!$    end if
    VPNR%ztr(1,2:kount+1)=-yp(2,1:kount)-1.d0
    VPNR%ztr(2,2:kount+1)=yp(1,1:kount)
    VPNR%ztr(3,2:kount+1)=xp(1:kount)

    VPNR%ztr(1,1)= - inivars(2) - 1.d0
    VPNR%ztr(2,1)= inivars(1)
    VPNR%ztr(3,1)= v0 

    VPNR%rtz(1,:)=VPNR%ztr(3,:)
    VPNR%rtz(2,:)=VPNR%ztr(2,:)
    VPNR%rtz(3,:)=VPNR%ztr(1,:)
    
    if(nvar>4)then
       VPNR%ztr(4,2:kount+1)=yp(4,1:kount)
       VPNR%ztr(5,2:kount+1)=yp(5,1:kount)
       
       VPNR%ztr(4,1)= inivars(4)
       VPNR%ztr(5,1)= inivars(5)
       
       VPNR%rtz(4,:)=VPNR%ztr(4,:)
       VPNR%rtz(5,:)=VPNR%ztr(5,:)       
    end if

    ! --------------------------------------- !
    ! Depricated lines......
    VPNR%ztrerr(:,1) = 0.d0
    VPNR%ztrerr(1,2:kount+1)=0.d0!yperr(2,1:kount)
    VPNR%ztrerr(2,2:kount+1)=0.d0!yperr(1,1:kount)
    VPNR%ztrerr(3,2:kount+1)=0.d0!yperr(3,1:kount)
    ! -------------------------------------- !

    ! -------------------------------------- !
    ! Get number of mins and maxs, for multiple
    ! outcomes for DA(z)
    

    call GetLocalMinMax_in_z()
    
    
    ! -------------------------------------- !

!    if(VoidError)return

!    call aspline(VPNR%ztr(1,:),VPNR%ztr(2:3,:),VPNR%ddztr(2:3,:))
    ! new 23/02/2011:
    if(VPNR%monotonic_in_z)then
       call aspline(VPNR%ztr(1,:),VPNR%ztr(2:nvar,:),VPNR%ddztr(2:nvar,:))
       !call maspline(VPNR%ztr(1,:),VPNR%ztr(2:3,:),VPNR%ddztr(2:3,:),VPNR%kmaxz,VPNR%kminz)
    else
       VPNR%ddztr=0.d0
    end if
    call aspline(VPNR%rtz(1,:),VPNR%rtz(2:nvar,:),VPNR%ddrtz(2:nvar,:))
    ! end new

!!$    if(VoidTestingIntegrationOutput)then
!!$       if(voidfeedback>2)write(*,*)'Re-Writing integration results to test_ltb.txt with ddy.'
!!$       open(123,file='test_ltb.txt')
!!$       do i = 1, kount
!!$          write(123,'(22ES14.5)')xp(i),yp(:,i), avoidS(xp(i),yp(1,i)), avoidSdot(xp(i),yp(1,i)), avoidRp(xp(i),yp(1,i)), avoidS(xp(kount),yp(1,i)), avoidSdot(xp(kount),yp(1,i)), voidR(xp(i),yp(1,i)),VPNR%ztr(:,i),VPNR%ddztr(2:3,i)
!!$       end do
!!$       close(123)
!!$       if(voidfeedback>2)write(*,*)'Closed test_ltb.txt.'
!!$    end if

    if(VoidTestingIntegrationOutput)then
!    if(.false.)then
       if(voidfeedback>2)write(*,*)'Writing DA(z) to test_ltb_DA.txt.'
       open(123,file='test_ltb_DA.txt')
       do i=1,500
!          write(123,*)void_zmax*dble(i-1)/499.,DA(void_zmax*dble(i-1)/499.)!,z_FLRW(NR(VP%zB*dble(i-1)/99.)*VP%Mtilde)
          write(123,*)min(void_zmax,3._dl)*dble(i-1)/499.,DA(min(void_zmax,3._dl)*dble(i-1)/499.)!,z_FLRW(NR(VP%zB*dble(i-1)/99.)*VP%Mtilde)
       end do
! High resolution (or change 500) and degenerate z:
!!$       allocate(testy(size(yp(:,1))), ddyy(size(yp(:,1)),kount))
!!$       call aspline(xp(1:kount), yp(:,1:kount), ddyy)
!!$       testimax = 500*kount
!!$       do i = 1, testimax
!!$          testr = (xp(kount) - xp(1))/(dble(testimax+100))*dble(i)
!!$
!!$          call asplint(xp(1:kount),yp(:,1:kount),testr,testy,ddyy)
!!$
!!$          dtm_tr(1) = testy(1)
!!$          dtm_tr(2) = testr
!!$          write(123,'(7E14.5)')testr,testy(1),-testy(2)-1.d0, DA(0.d0, dtm_tr)
!!$       end do
!!$       deallocate(testy,ddyy)
       close(123)
       if(voidfeedback>2)write(*,*)'Closed test_ltb_DA.txt.'
    end if

!    deallocate(vars,inivars)

  contains
  
    SUBROUTINE DoTheMagic_GFUN(NEQ, T, Y, NG, GOUT)
      implicit none
      integer, intent(in) :: neq, NG
      real(dl), intent(in) :: t, y(NEQ)
      real(dl), intent(out) :: GOUT(ng)
      real(dl) :: dy(neq)
!      real(dl) :: reladot, thissign
      
      if(ng/=Nevents)then
        call VoidErrorReport()
        write(0,*)'NG /= Nevents in gfun in voidasymp.f90'
        stop
      end if

      ! if a voiderror is issued at this point, it is to reject a model. Leave it.
      if(voiderror)then
        gout=0.d0
!      write(*,*)'gout',gout
!      voiderror = .false.
        return ! don't assess different cases.
      end if

      gout(1) = Y(stopvar) - thisstopval


      gout(2) = t - intvarexactlist(1)
      
      call DoTheMagic_derivswrap(neq,t,y,dy)

      ! If there is a recoverable error, voiderror is false, but dy=1.d30.
      ! continue but with small steps.
      if(any(abs(dy)>1.d29) .or. any(isnan(dy)) .or. notposdefsingeoy(y))then
        gout(3)=0.d0
!        write(*,*)dy
!        write(*,*)t,y
!        VoidError=.true.
      else 
        gout(3) = 1.d0
      end if

      if(y(1)<0.d0)then
        gout(4) = 0.d0
      else
        gout(4)=1.d0
      end if
      
      ! added 02/2012: E(r) can cross zero in case of nonzero kb:
      gout(5) = 1+2*t**2*VP%Mtilde**2*voidkofr(t)
      if(gout(5)<0.d0+1.d-8)then
        voiderror=.true.
        gout(5)=0.d0
        ! rejected.
      end if

      gout(6) = abs(t) - VP%L

    end subroutine DoTheMagic_GFUN

  end subroutine DoTheMagic

  function notposdefsingeoy(y)
    implicit none
    real(dl) :: y(:)
    logical :: notposdefsingeoy
    
    notposdefsingeoy = .false.
    if(y(1)<0.d0)notposdefsingeoy=.true.

  end function notposdefsingeoy

 subroutine DoTheMagic_derivswrap(neq,t,y,dy)
    implicit none
    integer, intent(in) :: neq
    real(dl), intent(in) :: t, y(neq)
    real(dl), intent(out) :: dy(neq)

    call ltb_aderivsdr(t,y,dy)
    if(VoidFeedback>99)then
      write(*,*)'r,y,dy:',t,y,dy
      write(*,*)'k,dk:',voidkofr(t),voiddkdr(t)
    end if
  end subroutine DoTheMagic_derivswrap


!> VoidSetL numerically inverts the relations between a certain coordinate size L and a redshift size zb.
!> It does so by calling DoTheMagic iteratively. Therefore, when VoidSetL succesfully finishes,
!> the final resulting distance integration is already stored in VPNR.
  subroutine VoidsetL
!    use AMLutils
    implicit none
    real(dl), save :: zBgoal, thisr, kmax
    real(dl) :: ltop, lbot, ThiszB, L
    real(dl) :: Lprecision=1.d-20!16!3
    integer, parameter :: nmem = 4!2
!    real(dl) :: prevR(nmem), prevL(nmem), prevZb(nmem), df, f, dLb, prevdf(nmem), ddf, ddldf2(nmem), sortedprevf(nmem), sortedprevL(nmem), prevf(nmem), GuessedL
    real(dl) :: prevR(nmem), prevL(nmem), df, f, dLb, prevdf(nmem), prevf(nmem), GuessedL
    integer :: i, j
    integer :: NewtRaphStart = 3
    logical :: LFound, CanUseSplint
! New size in mpc's:
  real(dL) :: mtildetompc = 299792.458d0 ! = c [km / s]
  real(dl) :: drdl(2,2)
  real(dl) :: aflrw0, alltb0, inzb, myacc, Sltb_L, Sltb_r0, GA_r_over_L, GA_rMpc_over_L_xSout
  logical :: CheckKofr
! end new
!    integer :: thetime(3)
!    call system_clock(count=thetime(1),count_rate=thetime(3))    

    if(VoidError)return


    inzb=VP%zb
    ! If VP%zb < 0, then VP%zb = -L in Mpc
    if(VP%zb<0)then

       VP%L=0.d0
       call VoidDistFuncs(obs="FLRW",r=0.d0,t=VP%t0,a=aflrw0)

!       aflrw0=voida(r=1.d-20,t=VP%t0)
       VP%L=-VP%zb
       alltb0=voida(r=0.d0,t=VP%t0)
       VP%L=VP%L / mtildetompc * VP%Mtilde_paper / VP%Mtilde !/ alltb0 * aflrw0! inversion of eq. in subroutine voiddipoleout in voidasymp.f90
        if(alltb0==0.d0 .or. VP%L>1.d300)then
             if(VOidFEedback>1)write(*,*)'a_lltb(r=0,t0) is zero. This indicates that there is something sick in wLLTBBackground initial time. Error.'
             VoidError=.true.
             return
        end if
       i=0
       do while (VoidCheckKofr())
          VP%L = VP%L *0.9d0
          i=i+1
          if(i>100)then
             if(VOidFEedback>1)write(*,*)'Cannot find healthy initial L to start VoidSetL.'
             VoidError=.true.
             return
          end if
       end do

       call setrtopini


       L = VP%L
       VP%deadon = L
       
       ! So far it is wrong. Now invert real distance:
       ltop = VP%rtopini

       lbot = 0.d0
       if(VoidError)return
       LFound = .false.
       CanUseSplint = .false.
       
       void_zmax = void_zmax_fixed
       thisr = 1.d10
       i=0
       dlb=1.d10
       drdl(:,1)=0.d0
       drdl(:,2)=1.d10

       Sltb_L = avoids(L,vp%t0)
       VP%r0 = L


       
       if(VoidZMaxType==0)then

          GA_r_over_L=VoidGARadius()/VP%L

       else if(VoidZMaxType==2)then

          GA_r_over_L=1.d0
          thisr = -inzb
          VP%L = thisr / mtildetompc * VP%Mtilde_paper / VP%Mtilde
       else

          GA_r_over_L=1.d0
          write(*,*)'L = LTB patch size.'
       end if

        ! deprecated function call.
       !GA_rMpc_over_L_xSout=VoidComovDistMPC(0.d0,GA_r_over_L * VP%L,VP%t0,1.1*VP%L)*Sltb_L/VP%L
       !            thisr=VoidComovDistMPC(0.d0,L,VP%t0,VP%r0)


       do while (min(dabs(Thisr/inzb+1.d0),dabs(dlb/L)).gt.1.d-15)
          ! Adaptive accuracy in RObs: first rough, but better when we approach the answer.
          myacc=1.d-1 * min(dabs(Thisr/inzb+1.d0),dabs(dlb/L),1.d0)
          myacc=max(myacc,voiddistprec)


!            thisr=VoidComovDistMPC(VP%r0,L,VP%t0)
          CheckKofr = VoidCheckKofr()
          if(.not.CheckKofr)then ! false if L is too big - checks that 1+2r^2Mt^2k(r) > 0 for all r

             if(nint(VoidRObs)==-1)then
                call VoidSetRObs(inzb,integration=.false.,acc=myacc)
             else if (nint(VoidRObs)==-2)then
                VP%r0=VP%L
             else if (nint(VoidRObs)==-3)then
                VP%r0=GA_r_over_L*VP%L
             else if (nint(VoidRObs)==-4)then
                VP%r0=-VP%L
             else if (nint(VoidRObs)==-5)then
                VP%r0=-GA_r_over_L*VP%L
             else if ((VoidRObs)==0.d0)then
                VP%r0=0.d0
             else if ((VoidRObs)>0.d0)then
                call VoidSetRObs(inzb,integration=.false.,acc=myacc)
             else
                stop 'Do not understand your value of VoidRObs'
             end if
             Sltb_r0 = avoids(vp%r0,vp%t0)
             ! Slow:
             ! thisr=VoidComovDistMPC(0.d0,VoidGARadius(),VP%t0,VP%r0)
             ! faster:
              thisr=VoidComovDistMPC(0.d0,GA_r_over_L * VP%L,VP%t0,VP%r0)
             ! fastest:
             ! thisr = GA_rMpc_over_L_xSout*VP%L/Sltb_r0
             ! accurate:
!!$             if(abs(thisr/inzb+1.d0)>1.d-6)then
!!$                thisr = GA_rMpc_over_L_xSout*VP%L/Sltb_r0
!!$             else
!!$                thisr=VoidComovDistMPC(0.d0,GA_r_over_L * VP%L,VP%t0,VP%r0)
!!$             end if

             drdl(1,2)=drdl(1,1)
             drdl(1,1)=thisr
             drdl(2,2)=drdl(2,1)
             drdl(2,1)=L

          else
             Sltb_r0 = 1.d0
             thisr = 1.d30
          end if

          ! Approximate Newton-Raphson works for this situation:
          if(i>newtraphstart)dlb=(-inzb - thisr)/(drdl(1,1)-drdl(1,2))*(drdl(2,1)-drdl(2,2))

          if(thisr>-inzb)then
             !if(i<newtraphstart .or. i > 500)
             ltop = L
             dlb = -abs(dlb)
             if(l+dlb<lbot .or. i < newtraphstart .or. i > 100 .or.isnan(dlb) )dlb=-0.5_dl*(ltop-lbot)
          else
             !if(i<newtraphstart .or. i > 500)
             lbot = L
             dlb=abs(dlb)
             if(l+dlb>ltop .or. i < newtraphstart .or. i > 100 .or.isnan(dlb))dlb=0.5_dl*(ltop-lbot)
          end if
          L = L + dlb
          VP%L = L
          i=i+1
          if(i>LMaxSteps)exit
          
            ! new 28/03/2011
            if(VoidFeedback>1)write(6,'(A)', advance="no")'.'
            write(*,'(20ES24.15)')thisr,inzb,l,ltop,lbot

            if(abs(lbot/ltop-1.d0)<1.d-10 .and. abs(abs(dlb)/abs(0.5_dl*(ltop-lbot))-1.d0)<1.d-16)then
               if ( dabs(Thisr/inzb+1.d0) > 1.d-2 )then
!                  write(*,*)'Fuck.'
                  if(abs(thisr)<abs(inzb))then
                     ltop = ltop * 2.d0
                  else
                     lbot = lbot / 2.d0
                  end if
               end if
            end if
         end Do
         If(VoidFeedback>1)write(*,*)'Found Robs and L in', i, 'steps.'
!         call VoidSetRObs(inzb,integration=.false.)
!            thisr=VoidComovDistMPC(0.d0,L,VP%t0,VP%r0)
         !write(*,*)'L:',VP%L
       if(SetLBoost .gt. 1.d-2)SetLBoost = SetLBoost / 1.d2 ! At higher precisionx
       If(VoidFeedback>1)write(*,*)'Repeating at higher precision.' 
!       if(VoidCheckKofr())stop 'this cannot be the case.'
        if(VoidCheckKofr())then ! actually it can, for some alpha we can simply never find L, and end up here with an error.
          voidError = .true.
          return
        end if
        call VoidSetRObs(inzb,integration=.true.)

!       call system_clock(count=thetime(2))
!       write(*,*)'Set L_GA in ',dble(thetime(2)-thetime(1))/dble(thetime(3)),'seconds.'

       ! Integration is done in voidsetrobs
!       call DoTheMagic
       
       ! but first get z_FLRW at L, exactly:
!       VP%deadon = L
!       call init_flrw
!       if(VoidError)return
!       call dothemagic
       if(VoidError)return

       if(maxval(VPNR%ztr(3,:))<VP%L)then
          ! L is not in integration domain
          ! which means that zmax = 1100 occurs
          ! before r=L. Last scattering then
          ! is inside the LTB metric, in which
          ! case we cannot calculate the cmb.
          VoidError=.true.
          return
       end if
       
       thiszB = NZ(VP%L)
       VP%zb=thiszb
       If(VoidFeedback>1)write(*,*)'New zb:',VP%zb

!       if((thiszb<0).or.(thiszb>void_zmax))then
       if((thiszb>void_zmax))then
          ! Probably L is not in integration domain
          ! which means that zmax = 1100 occurs
          ! before r=L. Last scattering then
          ! is inside the LTB metric, in which
          ! case we cannot calculate the cmb.
          VoidError=.true.
          call VoidErrorReport()
!          call wlltb_error('This error should not happen in VoidSetL')
          write(0,*)'This error should not happen in VoidSetL'
          stop
!          return
       end if

       zbgoal=VP%zb
       
       if (thiszb<0.) then ! observer is outside. Not interesting
          if(voidfeedback>1)write(*,*)'Rejecting models with observer outside: zb<0'
          VoidError=.true.
          return
       end if

    else
       ! If vp%zb>0:
      
       kmax = VP%kmax
       df = 0.d0

       call setrtopini

       if(VoidError)return
       LFound = .false.
       CanUseSplint = .false.

       ! For finding L we donot have to integrate up to void_zmax_fixed
       ! but we do need some more points, otherwise the interpolation is too inaccurate
       void_zmax = VP%zB*5.d0 

       Do while (.not. LFound)

          ltop = VP%rtopini
          lbot = 0.d0
          L = 1.d-1 * VP%rtopini
          L = 5.d-1 * VP%rtopini
          VP%L=L
          !    ThiszB = 0.d0
          thisr=1.d10
          zBgoal = VP%zB
          i=0


          do while (dabs(Thisr/L-1.d0).gt.Lprecision)
             CanUseSplint=.false.
             !Newton-Raphson:

             ! f = r(L) - L
             ! df = dr/dL - 1
             ! x = L
             ! dx = - f / df
             f = thisr - L
             if(i.ge.NewtRaphStart)df = (prevR(1)-thisr)/(prevL(2)-prevL(1)) - 1.d0 ! can we improve this?
             !          ddf = (prevdf(1) - df)/(prevL(2)-prevL(1))
             !          ddf = 0.d0
             !didn't help.

             if(thisr.lt.0.d0)then
                dLb = dLb! * 0.5d0
                exit ! force higher precision
             else if (Thisr.lt.L) then

                ltop = L
                !             if(i.gt.nmem .and. .not.(any(prevR.gt.1.d29)))then
                if(CanUseSplint)then ! Never used.
                   dLb = GuessedL - L
                else if(i.gt.NewtRaphStart)then
                   dLb = -dabs(f/df)
                   !                dLb = -dabs(f/(df+dLb*ddf))
                   if(L+dLb.lt.lbot )dLb = - 0.5d0*(ltop - lbot)
                else
                   dLb = -0.5d0*(ltop - lbot)             
                end if
                !          L = (ltop + lbot) * 0.5d0
                !          if(L.lt.lbot)lbot=L
             else if(thisr.gt.1.d29) then ! voiderror occured
                dLb = dLb!*0.5d0
                exit ! force higher precision
             else

                lbot = L

                !             if(i.gt.nmem .and. .not.(any(prevR.gt.1.d29)))then
                if(CanUseSplint)then ! Never used
                   dLb = GuessedL - L
                else if(i.gt.NewtRaphStart)then
                   dLb = dabs(f/df)
                   !               dLb = dabs(f/(df+dLb*ddf))
                   if(L+dLb.gt.ltop )dLb = 0.5d0*(ltop - lbot)
                else
                   dLb = 0.5d0*(ltop - lbot)             
                end if
                !          L = (ltop + lbot) * 0.5d0
             end if
             if(isnan(dLb))write(*,*)f,df
             if(L+dLb.lt.lbot )dLb = - 0.5d0*(ltop - lbot)
             !          if(L+dLb.gt.ltop )dLb = - 0.5d0*(ltop - lbot)
             ! WV new for collapsing scenarios
             !          if(L+dLb.gt.VP%rtopini )dLb = - 0.5d0*(VP%rtopini - lbot)

             L = L + dLb

             if(prevL(1)==L.and. i .gt.NewtRaphStart)then ! df = infty
                exit
             end if
             VP%L=L
             VP%deadon = L

             call DoTheMagic


             if(VoidError)then
                thisr=1.d30
                VoidError=.false.
             else
                thisr=NR(zBgoal)
             end if

             if(dabs(ltop/lbot-1.d0).lt.Lprecision)exit
             if(ltop .eq. lbot)then
                write(*,'(6E14.5)')lbot,l,ltop,VP%kmax,VP%zB,VP%H0
                write(*,*) 'Voiddistances cannot find L(zB).'
                if(dabs(Thisr/L-1.d0).lt.Lprecision)then
                   exit
                else
                   VoidError=.true.
                   return
                end if
             end if
             i=i+1
             if(i.gt.Lmaxsteps) then
                if(voidfeedback>2)write(*,*) 'Voiddistances got stalled while searching L for r.'
                !             VoidError=.true.
                !             return
                exit ! force higher precision
             end if

             do j = nmem, 2, -1
                prevR(j) = PrevR(j-1)
                prevL(j) = PrevL(j-1)
                prevf(j) = Prevf(j-1)
                prevdf(j) = Prevdf(j-1)
             end do
             prevR(1) = thisr
             prevL(1) = L
             prevf(1) = f
             prevdf(1) = df

             !          CanUseSplint = ((i.lt.NewtRaphStart).and.(i.gt.nmem) .and. .not.(any(prevR.gt.1.d29)) .and. ( (maxval(sortedprevf).gt.0.d0) .and. (minval(sortedprevf).lt.0.d0)))
             CanUseSplint = .false. ! It doesn't help at all.
!!$          if(CanUseSplint)then
!!$             sortedprevL = prevL
!!$             sortedprevf = prevf
!!$             call voidsort(sortedprevf,sortedprevL)
!!$             
!!$             call spline_double(sortedprevf,sortedprevL,nmem,ddldf2)
!!$             
!!$             GuessedL= voidsplint(sortedprevf,sortedprevL,ddldf2,0.d0)
!!$          end if

             ! Wv new: difficulty with collapsing scenarios:
             thiszB = NZ(VP%L)
             if((dabs(thiszB/zBgoal - 1.d0).gt.VoidZBPrecision*1.d2))then
                if(abs(L/Lbot-1.d0)<1.d-12)lbot=lbot*0.1d0
                if(abs(L/Ltop-1.d0)<1.d-12)ltop=0.5_dl*(VP%rtopini+ltop)

             end if
             !          writE(*,*)'bot:',abs(L/Lbot-1.d0)
             !          writE(*,*)'top:',abs(L/Ltop-1.d0)

          end do

          if(voidfeedback.gt.0)write(*,*)'Found L in ',i,'steps.'



          VP%L = L
          if(L.eq.0.d0)then
             VoidError=.true.
             return
          end if

          if(VoidError)return

          ! Check if this result is actually good, if not, repeat at higher precision
          thiszB = NZ(VP%L)
          if(dabs(thiszB/zBgoal - 1.d0).gt.VoidZBPrecision)then
             LFound = .false.
             SetLBoost = SetLBoost / 1.d1
             if(voiderror)then
                if(voidfeedback.gt.0)write(*,*)'Voiderror and inaccurate L(zB): giving up.'
                return
             end if
             if(SetLBoost*VoidDistPrec.lt.VoidMinPrec)then
                voiderror = .true.
                exit
             end if
             if(voidfeedback.gt.0)then
                write(*,*)'Searching L at higher precision:', SetLBoost
                write(*,*)'So far we had:',dabs(thiszB/zBgoal - 1.d0)  
                write(*,*)'Found z(L), desired z(L):', thiszB, zBgoal
                write(*,*)'Lbot,L,Ltop:', 0.d0,VP%L,VP%rtopini
                write(*,*)'L(zB):',NR(zBgoal)
             end if
          else
             LFound = .true.
             exit
          end if
       end Do


       ! For finding L we donot have to integrate up to void_zmax_fixed
       ! Now we do need that.
       ! Repeat now for full range:
       void_zmax = void_zmax_fixed
       if(SetLBoost .gt. 1.d-2)SetLBoost = SetLBoost / 1.d2 ! At higher precisionx

       If(VoidFeedback>1)write(*,*)'Repeating init_flrw.'
       ! but first get z_FLRW at L, exactly:
       VP%deadon = L
       call init_flrw
       if(VoidError)return


       If(VoidFeedback>1)write(*,*)'Repeating dothemagic.'
       call dothemagic
       if(VoidError)return

       thiszB = NZ(VP%L)
       !    write(*,*)'At last higher accuracy run:',dabs(thiszB/zBgoal - 1.d0)

       If(VoidFeedback>1)write(*,*)'Final call for robs.'
       ! new 28/03/2011
       call VoidSetRObs(inzb)

    end if

    if(voidfeedback.gt.1)write(*,*)'L:',VP%L
    VP%L_paper=VP%L * VP%t0bg / VP%t0
    if(voidfeedback.gt.1)write(*,*)'L (paper dimensions):',VP%L_paper!VP%L * VP%t0bg / VP%t0
    if(voidfeedback.gt.1)write(*,*)'z(L)/zB - 1:',thiszB/zBgoal - 1.d0

!    if(VoidWantErrors)call EstimateErrors



!!$    ! DEBUGGING:
!!$    do i=1,100
!!$       VP%L=L*(1+dble(i-50)/100.d0/5.d0)
!!$       call DoTheMagic!(VP,VPNR)
!!$!       if(VoidError)return
!!$       if(VoidError)then
!!$          thisr=1.d30
!!$          VoidError=.false.
!!$       else
!!$          thisr=NR(zBgoal)
!!$       end if
!!$       write(*,*)VP%L, thisr, zBgoal
!!$    end do
!!$stop 'Debugged?'
    !       write(*,*)'finding L. L, r:',L, thisr

    !stop ' set l'


  end subroutine VoidsetL


!> setrtopini: find turn-over point for z(r)
  subroutine setrtopini
    implicit none
    real(dl) :: z,z1,z2!, realzB
    real(dl) :: r,rtopini
    integer :: j
    
    if(VoidError)return

    call init_FLRW()
    
    
    ! first find rtopini
    z2=1.d-30
    z1=0
    rtopini=1.d0
    r=1.d-2/VP%Mtilde
    j=0

    if(3.d0*VP%zB.gt.void_zmax)stop 'Increase void_zmax to deal with this z_B. function z_FLRW.'
    do while ((z2.lt.3.d0*VP%zB))
       rtopini=r
       r=r*1.1
       
       z = z_FLRW(r)
       if(z.gt.1.d29)then
          VoidError=.false.
          rtopini=rtopini/1.1/1.1 ! two steps back. Should be fine, right?
          exit
       end if
       z1=z2
       z2=z
       j=j+1
       
       if(j.gt.maxsteps)then
          write(*,*)'z,r, VP%H0, VP%Hlocal, VP%kmax, VP%zB:'
          write(*,'(9E14.5)')z,r, VP%H0, VP%Hlocal, VP%kmax, VP%zB
          write(*,*) 'Find rtopini got stalled.'
          VoidError=.true.
          exit
       end if

    end do

    VP%rtopini = rtopini
    

  end subroutine setrtopini
  
  subroutine init_FLRW
    implicit none
    integer :: mystat
    real(dl) :: Lmem

    if(voiderror)return

    VP%ini_flrw = .true.
    
    Lmem = VP%L

    VP%L= 0.d0
    if(Lmem.le.0.d0)then
!!$       VP%deadon = Lmem
!!$    else
       VP%deadon = 0.d0
    end if

    call dothemagic

!    if(VoidWantErrors)call EstimateErrors


!    if(associated(VPNR_FLRW%rtz))deallocate(VPNR_FLRW%rtz,stat=mystat)
!    if(associated(VPNR_FLRW%ddrtz))deallocate(VPNR_FLRW%ddrtz,stat=mystat)

!    if(associated(VPNR_FLRW%ztr))deallocate(VPNR_FLRW%ztr,stat=mystat)
!    if(associated(VPNR_FLRW%ddztr))deallocate(VPNR_FLRW%ddztr,stat=mystat)
!    if(associated(VPNR_FLRW%ztrerr))deallocate(VPNR_FLRW%ztrerr,stat=mystat)
!    if(associated(VPNR_FLRW%kmaxz))deallocate(VPNR_FLRW%kmaxz,stat=mystat)
!    if(associated(VPNR_FLRW%kminz))deallocate(VPNR_FLRW%kminz,stat=mystat)
    call deallocateThisVPNR(VPNR_FLRW)
  
!    allocate(VPNR_FLRW%rtz(size(VPNR%ztr(:,1)),size(VPNR%ztr(1,:))))
!    allocate(VPNR_FLRW%ddrtz(size(VPNR%ddztr(:,1)),size(VPNR%ddztr(1,:))))

!    allocate(VPNR_FLRW%ztr(size(VPNR%ztr(:,1)),size(VPNR%ztr(1,:))))
!    allocate(VPNR_FLRW%ddztr(size(VPNR%ddztr(:,1)),size(VPNR%ddztr(1,:))))
!    allocate(VPNR_FLRW%ztrerr(size(VPNR%ztr(:,1)),size(VPNR%ztr(1,:))))
!    allocate(VPNR_FLRW%kmaxz(VPNR%nmaxz))
!    allocate(VPNR_FLRW%kminz(VPNR%nminz))
    
!    VPNR_FLRW%rtz = VPNR%rtz
!    VPNR_FLRW%ddrtz = VPNR%ddrtz

!    VPNR_FLRW%ztr = VPNR%ztr
!    VPNR_FLRW%ddztr = VPNR%ddztr
!    VPNR_FLRW%ztrerr = VPNR%ztrerr
!    VPNR_FLRW%kmaxz = VPNR%kmaxz
!    VPNR_FLRW%kminz = VPNR%kminz
!    VPNR_FLRW%nmaxz = VPNR%nmaxz
!    VPNR_FLRW%nminz = VPNR%nminz

!    VPNR_FLRW%monotonic_in_z = VPNR%monotonic_in_z
    call copyVPNRToThisVPNR(VPNR,VPNR_FLRW)

    VP%L = Lmem


    VP%ini_flrw=.false.

 end subroutine init_FLRW

!> GetLocalMinMax_in_z determines the indices in VPNR's results, for which the redshift has a local minimum or maximum.
!> These indices are used when calling the multivalued splines later on.
  subroutine GetLocalMinMax_in_z()
    implicit none
    integer, allocatable :: kminima(:),kmaxima(:)
    real(dl), allocatable :: x(:)
    integer :: n,mystat

    n=size(VPNR%ztr(1,:))
    allocate(kminima(n),kmaxima(n),x(n))

    x=VPNR%ztr(1,:)


    call GetLocalMinMax(x,kmaxima,kminima,VPNR%nmaxz,VPNR%nminz)



    if(associated(VPNR%kmaxz))deallocate(VPNR%kmaxz,stat=mystat)
    if(associated(VPNR%kminz))deallocate(VPNR%kminz,stat=mystat)

    allocate(VPNR%kmaxz(VPNR%nmaxz))
    allocate(VPNR%kminz(VPNR%nminz))
    
    VPNR%kmaxz = kmaxima(1:VPNR%nmaxz)
    VPNR%kminz = kminima(1:VPNR%nminz)
    
    if(VPNR%nmaxz>1.or.VPNR%nminz>1)then
       VPNR%monotonic_in_z = .false.
    else
       VPNR%monotonic_in_z = .true.
    end if

    deallocate(kminima,kmaxima,x)

  end subroutine GetLocalMinMax_in_z

!> VoidSetRobs is used in combination with VoidSetL. For each iteration of VoidSetL, VoidSetRObs is called to 
!> place the observer at the right coordinate, corresponding to a certain comoving distance in the observer's frame.
  subroutine VoidSetRobs(inzb,integration,acc)
    use NROdeInt
    use ode_path
    implicit none
    real(dl) :: inzb
    logical, optional :: integration
    real(dl), optional :: acc
    logical :: dointegration
    integer :: i
    real(dl) :: thisroverl, myobst
!    real(dl), external :: rombint
    real(dl) :: testl
    real(dl) :: thisr,dlb,aflrw0,alltb0, rtopbot(2),ltopbot(2)
!    integer, parameter :: NewtRaphStart = 1, raphoff=20
    integer  :: NewtRaphStart = 1, raphoff=50 ! raphson doesn't work
    real(dl) :: RObsAccuracy = VoidDistPrec
    real(dl), parameter :: lfrac = 1.d-1 ! Tuned: 0.1 seems fastest.
    real(dl) :: thisrobs, thisrobs_c, GoalRobs, Sobs, Sprime, thisrmax
    real(dl) :: sofr(1), thissofr(2), myacc
    real(dl), allocatable :: sofrres(:,:), ddsofr(:,:)
    integer :: mymin, mymax
! integer :: thetime(3)
! call system_clock(count=thetime(1),count_rate=thetime(3))    

    myobst=VP%t0
    if(present(acc))then
       RObsAccuracy = acc
    else
       RObsAccuracy = VoidDistPrec
    end if
    

    if(VoidRObs==0.d0)then
       thisr=0.d0
       VP%r0 = thisr
       GoalRobs=-1.
       
    else if(nint(VoidRObs)==-1)then
       if(VP%d0_2/=0.)then
          VP%r0 = VP%d0_2 * VP%L
          thisr=VoidComovDistMPC(0.d0,VP%r0,VP%t0,VP%r0)
          GoalRobs=-1.
       else if(VP%d0_2==0.d0)then
          VP%r0 = 0.d0
          thisr=0.d0
          GoalRobs=-1.
       end if
    else if(nint(VoidRObs)==-2)then
       VP%r0 = VP%L
       thisr=VoidComovDistMPC(0.d0,VP%r0,VP%t0,VP%r0)
       GoalRobs=-1.
    else if(nint(VoidRObs)==-3)then
       VP%r0 = VoidGARadius()
       thisr=VoidComovDistMPC(0.d0,VP%r0,VP%t0,VP%r0)
       GoalRobs=-1.
    else if(nint(VoidRObs)==-4)then
       VP%r0 = -VP%L
       thisr=VoidComovDistMPC(0.d0,VP%r0,VP%t0,VP%r0)
       GoalRobs=-1.
    else if(nint(VoidRObs)==-5)then
       VP%r0 = -VoidGARadius()
       thisr=VoidComovDistMPC(0.d0,VP%r0,VP%t0,VP%r0)
       GoalRobs=-1.
    else if (VoidRObs>0.d0)then
       GoalRobs = VoidRObs
    else
       write(*,*)'VoidRobs:',voidrobs
       write(0,*)'VoidRobs:',voidrobs
       stop 'What the fuck is VoidROBs?'
    end if
    
    if(GoalRobs>0)then
       !    if(inzb<0.d0)stop 'RObs/=0 only works with zb>0.'
       
!       if(.true.)then
          ! new 28/04/2011: try odeint:
          ! Nope, is slower, but much more accurate.
          sofr=0.
          thisrmax = 1.d1 * GoalRobs / (VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde)
          !       myacc=RObsaccuracy/VoidDistPrec * 1.d-10
          myacc=VoidDistPrec
          nrsuccess=.true.
          call odeint(sofr,0.d0, (/thisrmax/),0,thisrmax,myacc,myacc,0.d0,RObsDerivs,rkqs)
          if(.not.nrsuccess .or. kount<2)then
             voiderror=.true.
             return
          end if

          allocate(sofrres(2,1+kount))
          allocate(ddsofr(2,1+kount))
          
          sofrres(1,2:kount+1)=xp(1:kount)
          sofrres(2,2:kount+1)=yp(1,1:kount)
          
          sofrres(:,1)=0.d0
          
          ltopbot(1) = 1.d10
          ltopbot(2) = 0.d0
          do i = 1, kount+1
             sofrres(2,i) = sofrres(2,i) * VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde / avoidS(sofrres(1,i),vp%t0)
             if(sofrres(2,i)<GoalRobs)then
                ltopbot(2) = sofrres(1,i)
                rtopbot(2) = sofrres(2,i)
             else
                if(ltopbot(1)>sofrres(1,i))then
                   ltopbot(1)=sofrres(1,i)
                   rtopbot(1)=sofrres(2,i)
                end if
             end if
          end do

          ! The rescaling doesn't do any good to the accuracy of the results.
          ! Zeros may appear in the differences. Aspline cannot deal with that.
          ! Solution:
          mymin = 1
          mymax = kount+1
          if (size(sofrres(1,:))<mymax)then
             voiderror=.true.
             return
          end if
          do while (any(sofrres(:,mymin+1)-sofrres(:,mymin)==0.d0))
             mymin = mymin+1
          end do

          do while (any(sofrres(:,mymax)-sofrres(:,mymax-1)==0.d0))
             mymax = mymax-1
          end do
          ! End solution.

          
          call aspline(sofrres(2,mymin:mymax),sofrres,ddsofr(:,mymin:mymax))
          call asplint(sofrres(2,mymin:mymax),sofrres,GoalRobs,thissofr,ddsofr(:,mymin:mymax))
          VP%r0 = thissofr(1)

          deallocate(sofrres)
          deallocate(ddsofr)
          
!!$       else
!!$          ! new 28/03/2011
!!$          aflrw0=voida(r=1.d-20,t=VP%t0)
!!$          alltb0=voida(r=0.d0,t=VP%t0)
!!$          VP%r0=GoalRObs / VPc/1.d-3 * VP%Mtilde_paper / VP%Mtilde / alltb0 * aflrw0
!!$          
!!$          VP%r0 = 0.d0
!!$          thisrobs = 0.d0
!!$          thisrobs_c = 0.d0
!!$          i=0
!!$          do while ( thisrobs < GoalRobs)
!!$             ltopbot(2) = VP%r0
!!$             rtopbot(2) = thisrobs
!!$             VP%r0 = VP%r0 + lfrac * VP%L
!!$             thisrobs_c = thisrobs_c + lfrac * VP%L * myintegrand(VP%r0) * VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde
!!$             thisrobs = thisrobs_c / avoidS(VP%r0,vp%t0)
!!$             if(thisrobs>GoalRobs)exit       
!!$             i=i+1
!!$             if(i>LMaxSteps)then
!!$                voiderror = .true.
!!$                if(voidfeedback>0)write(*,*)'Cannot find RObs.'
!!$                exit
!!$             end if
!!$          end do
!!$          ltopbot(1) = VP%r0
!!$          rtopbot(1) = thisrobs
!!$          !write(*,*)thisrobs,ltop,lbot
!!$          !    ltop = VP%rtopini
!!$          !    lbot = 0.d0
!!$          if(VoidError)return
!!$          
!!$          !    thisr = 1.d10
!!$          thisr = thisrobs
!!$          
!!$!    ltopbot(1)=ltop
!!$!    ltopbot(2)=lbot
!!$!    rtopbot(1)=thisr
!!$!    rtopbot(2)=0.d0
!!$          
!!$       end if
       
       i=0
       dlb=1.d10    
       thisr=VoidComovDistMPC(0.d0,VP%r0,VP%t0,VP%r0)
       
!!$       if(.false.) then     
!!$          ltopbot(1) = max(VP%L,5*VP%r0)
!!$          do while (min(dabs(Thisr/GoalRObs-1.d0),dabs(dlb/VP%r0)).gt.RObsAccuracy)
!!$             Sprime = (avoids(VP%r0*1.01,VP%t0)-avoids(VP%r0*0.99,VP%t0))/(0.02*VP%r0)
!!$             ! get sobs last: then this value is cached. Will be called again straight after this 
!!$             ! routine in voidsetl. Safes a fraction of a fraction of a second.. (about 0.05 sec.)
!!$             Sobs = avoids(VP%r0,VP%t0)
!!$             
!!$             dlb = (VoidRObs - thisr)/ VPc/1.d-3 * VP%Mtilde_paper / VP%Mtilde &
!!$                  / ( 1._dl - Sprime / Sobs * VoidRobs * ( VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde) )
!!$             !          if(ltopbot(1)==VP%rtopini .or. ltopbot(2)==0.d0)NewtRaphstart=NewtRaphStart+1
!!$             if(thisr>GoalRObs)then
!!$                
!!$                !   dlb = (GoalRObs - thisr) * (VP%r0-ltopbot(2)) / (thisr-rtopbot(2))
!!$                
!!$                rtopbot(1)=thisr
!!$                ltopbot(1)=VP%r0
!!$                if(i<newtraphstart .or. i > raphoff)ltopbot(1) = VP%r0
!!$                dlb = -abs(dlb)
!!$                if(thisr+dlb<ltopbot(2) .or. i < newtraphstart .or. i > raphoff)then
!!$                   dlb=-0.5_dl*(ltopbot(1)-ltopbot(2))
!!$                end if
!!$                !             if(i < newtraphstart .or. i > raphoff)dlb=-0.5_dl*(ltop-lbot)
!!$             else
!!$                
!!$                !  dlb = (GoalRObs - thisr) * (ltopbot(1)-VP%r0) / (rtopbot(1)-thisr)
!!$                
!!$                rtopbot(2)=thisr
!!$                ltopbot(2)=VP%r0
!!$                dlb=abs(dlb)
!!$                if(i<newtraphstart .or. i > raphoff)ltopbot(2) = VP%r0
!!$                if(thisr+dlb>ltopbot(1) .or. i < newtraphstart .or. i > raphoff)then
!!$                   dlb=0.5_dl*(ltopbot(1)-ltopbot(2))
!!$                end if
!!$                !             if(i < newtraphstart .or. i > raphoff)dlb=0.5_dl*(ltop-lbot)
!!$             end if
!!$             VP%r0 = VP%r0 + dlb
!!$             i=i+1
!!$             if(i>LMaxSteps)then
!!$                voiderror = .true.
!!$                if(voidfeedback>0)write(*,*)'Cannot fin RObs.'
!!$                exit
!!$             end if
!!$             !write(*,'(10ES24.15)')thisr,voidrobs,vp%r0,vp%l,ltop,lbot
!!$             thisr=VoidComovDistMPC(0.d0,VP%r0,VP%t0,VP%r0)
!!$          end Do
!!$       end if
    end if ! Now r0 is set for all cases.

    dointegration = .true.
    if(present(integration))then
       if(.not.integration)dointegration=.false.
    end if
    
!!$    if(voidfeedback>1)then
!!$       write(*,*)'r0:',VP%r0
!!$       write(*,*)'L:',VP%L
!!$       write(*,*)'kmax:',VP%kmax
!!$       write(*,*)'precision',VoidDistPrec
!!$       write(*,*)'boost:',setlboost
!!$    end if


    if(dointegration)then
       
       
       VP%deadon = VP%L
! why do this twice?
!       call init_flrw
       if(VoidError)return
       
    
       call dothemagic
       if(VoidError)return
       if(voidfeedback>1)then
          write(*,*)'RObs/=0. RObs:',thisr
          testl=nz(vp%L)
          write(*,*)'New zb:',testl
          testl=VoidComovDistMPC(VP%r0,VP%L,VP%t0,VP%r0)
          write(*,*)'Nearest border (negative means we''re outside):',testl
          testl=VoidComovDistMPC(-VP%r0,VP%L,VP%t0,VP%r0)
          write(*,*)'Furthest border:',testl
          testl=VoidComovDistMPC(0.d0,VP%L,VP%t0,VP%r0)
          write(*,*)'New L(Mpc):',testl

          testl=VoidComovDistMPC(0.d0,VoidGARadius(),VP%t0,VP%r0)

          write(*,*)'New L_GA(Mpc):',testl
          testl = nz(r=VoidGARadius())
          write(*,*)'New z_GA(Mpc):', testl
       end if
       
    end if

!    call system_clock(count=thetime(2))
!    write(*,*)'Set Robs in ',dble(thetime(2)-thetime(1))/dble(thetime(3)),'seconds.'

!!$  contains
!!$    
!!$    function MyIntegrand(r)
!!$      real(dl) :: r
!!$      real(dl) :: MyIntegrand
!!$      
!!$      MyIntegrand = aVoidS(r,vp%t0) 
!!$      
!!$    end function MyIntegrand
!!$    
    

  end subroutine VoidSetRobs

!> RObsDerivs is used for the integration performed in VoidSetRobs.
  subroutine RObsDerivs(x,y,dy)
    implicit none
    real(dl), intent(in) :: x, y(:)
    real(dl), intent(out) :: dy(:)
    integer :: n
    
    n = size(y)
    if(n/=1)stop 'wrong size in MyDerivs in VoidSetRobs.'
    
    dy(1) = aVoidS(x,VP%t0)
    
    
  end subroutine RObsDerivs
