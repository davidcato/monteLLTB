
   subroutine VoidBaoTofZ(zini,zrec, trbao, trrec)
     use NROdeint
     real(dl) :: zini, zrec, z
!     real(dl) :: tr(1:2), trbao(1:2), trrec(1:2), y(1), intvar, GetExactList(1)
     real(dl) :: tr(1:2), trbao(1:2), trrec(1:2), y(1), GetExactList(1)
     ! For volume element:
     real(dl) :: tmin, tmax, t, r, thisz, u(3), myS, myR, delt, iniS, iniR, tmem, zmem
     integer :: bcounter
     integer, parameter :: NewtRaphStart = 2
     ! To choose between choices:
     logical, parameter :: ImaMonkey = .false., YoureaMonkey=.false.!.true.
     logical :: isan=.true.
     z = zini

     call asplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),z,tr,VPNR%ddztr(2:3,:))

     trbao = tr

     if(ImaMonkey)then ! then I'm wrong.
        
        call VoidBaoRememberThisR(trbao(2),"put")
        
        y(1) = trbao(1)
!!$     GetExactList(1) =  zrec
!!$     call odeint(y, zini, GetExactList, 0, GetExactList(1), VoidDistPrec*1.d-4, VoidDistPrec, 0.d0, VoidBaodtdz, rkqs)
        GetExactList(1) =  -(1.d0 + zrec)
        call odeint(y, -(1.d0 + zini), GetExactList, 0, GetExactList(1), VoidDistPrec*1.d-4, VoidDistPrec, 0.d0, VoidBaodtdz, rkqs)
        
        
        trrec(1) = y(1)
        trrec(2) = trbao(2)
        
        
        call VoidBaoRememberThisR(z,"rst")
        
     else if(YoureaMonkey)then ! then you're wrong.
        void_zmax = zrec + 5.d0
        call DoTheMagic

        call asplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),zrec,trrec,VPNR%ddztr(2:3,:))
        trrec(2) = trbao(2)
     else ! else we're right.
        ! Or do we just want the growth of the volume element?
        bcounter = 0
        tmax = VP%t0
        tmin = 1.d-10
        !t = trrec(1)
        t = (tmax + tmin) * 0.5d0
        iniS = avoidS(trbao(2),trbao(1))
        iniR = voidR(trbao(2),trbao(1))
        isan=.true.
        do while (0 .ne. 1)
           r = trbao(2)
!           u = voidu(r,t)
           myS = avoidS(r,t)
           myR = voidR(r,t)
           thisz = ((iniS * iniR * iniR)/(myS * myR * myR))**(1.d0/3.d0)*(1+zini) - 1.d0

           if(thisz.gt.zrec)then ! t to small
              tmin = t
              delt = (tmax - tmin) * 0.5d0
              if(isan)then
                 if(bcounter>NewtRaphStart)then
                    delt = dabs( (thisz - zrec) / (thisz-zmem) * (t-tmem) )
                    if(isnan(delt))isan=.false. ! No more Newton-Raphson if we encounter thisz==zmem.
                    if((t + delt .gt. tmax).or.isnan(delt))  delt = (tmax - tmin) * 0.5d0
                 end if
              end if
           else ! t to large
              tmax = t
              delt = (tmin - tmax) * 0.5d0
              if(isan)then
                 if(bcounter>NewtRaphStart)then
                    delt = - dabs( (thisz - zrec) / (thisz-zmem) * (t-tmem))
                    if(isnan(delt))isan=.false. ! No more Newton-Raphson if we encounter thisz==zmem.
                    if((t + delt .lt. tmin).or.isnan(delt))  delt = (tmin - tmax) * 0.5d0
                 end if
              end if
           end if

           if(dabs(thisz / zrec - 1.d0) .lt. VoidDistPrec)exit
           if(.not.isan.and.(abs(thisz-zmem)<1.d-30))exit ! Because then it is the best we can do..

           tmem = t
           zmem = thisz
           t = t + delt
           bcounter = bcounter + 1
!write(*,'(4ES15.6,L3,10ES15.6)')t,tmax,tmin,delt,isan,dabs(thisz / zrec - 1.d0),thisz,zmem,thisz-zmem
           if(bcounter>maxsteps)then
              VoidError = .true.
              t = -1.d30
              exit
           end if
        end do
        trrec(1) = t
        trrec(2) = trbao(2)
     end if

 
   end subroutine VoidBaoTofZ


   subroutine VoidBaoRememberThisR(r,action)
     character(len=3) :: action
     real(dl) :: r
     real, save :: TheR
     logical, save :: isset = .false.

     if (action=='put') then
        TheR = r
        isset = .true.
     else if (action=='get') then
        if(isset)then
           r = TheR
        else
           stop 'VoidBaoRememberThisR is called without r being set.'
        end if
     else if(action=='rst') then
        TheR = -1.d30
        isset = .false.
     else
        write(*,*)'Syntax error in VoidBaoRememberThisR: '//action
        stop 'Syntax error in VoidBaoRememberThisR'
     end if


   end subroutine VoidBaoRememberThisR
     
   
   ! v = -(1+z)
   ! y(1) = t
   ! dy(1) = dt/dv
   subroutine VoidBaodtdz(v,y,dy)
     real(dl), intent(in) :: v
     real(dl), intent(in) :: y(:)
     real(dl), intent(out) :: dy(:)
     real(dl) :: r, t, u(3)

     call VoidBaoRememberThisR(r,"get")

     t = y(1)
     u = voidu(r,t)

     dy(1) = - avoidS(r,t,u) / avoidSdot(r,t,u) / v
!     dy(1) = -avoidS(r,t,u) / avoidSdot(r,t,u) / (1.d0+v)

   end subroutine VoidBaodtdz


   subroutine VoidDistFuncs(z,t,r,Rltb,Rp,Rpd,S,Sd,a,H,Mtilde,Mtilde_paper, obs)
     use wlltb
     implicit none
     
     real(dl), optional :: z,t,r,Rltb,Rp,Rpd,S,Sd,H,Mtilde,Mtilde_paper,a
     real(dl) :: u(3), tr(1:2), z1, k, Lmem
     character(len=*), optional :: obs
     logical :: FLRW

     FLRW = .false.
     if(present(obs))then
        if(obs=='FLRW')then
           FLRW=.true.
           Lmem = VP%L
           VP%L = 0.d0
        end if
     end if

     if(present(z))then
        z1 = z
        if(.not.FLRW)then
           call asplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),z,tr,VPNR%ddztr(2:3,:))
        else
           call asplint(VPNR_FLRW%ztr(1,:),VPNR_FLRW%ztr(2:3,:),z,tr,VPNR_FLRW%ddztr(2:3,:))
        end if
     else if (present(t) .and. present(r))then
        tr(1) = t
        tr(2) = r
     else
        stop 'Need to specify either z or t and r in VoidDistFunc.'
     end if

!     u = voidu(tr(2),tr(1))

     k = Voidkofr(tr(2))

     if(present(Rltb))Rltb = voidR(tr(2),tr(1))!,u)
     if(present(Rp))Rp = avoidRp(tr(2),tr(1))!,u)
!     if(present(Rpd))Rpd = avoidSdot(tr(2),tr(1))  * dsqrt( 1.d0 + 2.d0 * tr(2)**2 * k * VP%Mtilde**2  )
      if(present(Rpd))  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=tr(2),t=tr(1),Rpd=Rpd)

     if(present(S))S = avoidS(tr(2),tr(1))!,u)
     if(present(H))H = voidH(tr(2),tr(1))!,u)
     if(present(a))a = voida(tr(2),tr(1))!,u)
     if(present(Sd))Sd = avoidSdot(tr(2),tr(1))!,u)
     if(present(Mtilde))Mtilde = VP%Mtilde
     if(present(Mtilde_paper))Mtilde_paper = VP%Mtilde_paper

     if(FLRW)then
        VP%L = Lmem
     end if


   end subroutine VoidDistFuncs
   
   
   ! vdi_bao_Dv returns Dv to a point in redshift space, given the LTB metric / observer.
   function vdi_bao_Dv(z)
    implicit none
    real(dl) :: vdi_bao_Dv
    real(dl), intent(in) :: z
    ! oops: Y -> R
    real(dl) :: thisY, thisYp, thisYpd
    real(dl), parameter :: third = 1.d0/3.d0
    real(dl) :: tr(2), properDv, comovDv, properDv3
    
    if(VoidError)then
      vdi_bao_Dv=1.d-30
      return
    end if
    
    thisY = DA(z)
    
    if(abs(VP%zb)>1.d-16)then
       call asplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),z,tr,VPNR%ddztr(2:3,:))
    else
       call asplint(VPNR_FLRW%ztr(1,:),VPNR_FLRW%ztr(2:3,:),z,tr,VPNR_FLRW%ddztr(2:3,:))
    end if

    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
          r=tr(2),t=tr(1),&
          Rp=thisYp,Rpd=thisYpd)

    properDv3 = (VPc/1.d3)* z/(1.d0+z)*thisYp*thisY**2/thisYpd
    properDv = properDv3**third

    comovDv = properDv * vdi_bao_proptocomov3D(tr(1),tr(2))

    vdi_bao_Dv = comovDv

   end function vdi_bao_Dv  
   
    ! vdi_bao_proptocomov3D gives the effective scalefactor that 
    ! relates a proper length at some time t, positioned at r, to 
    ! a comoving size, given that this length actually defines then side
    ! of a 3D volume L^3. FOr example for the sound horizon that appears
    ! in the 3D Bao: 1 radial scaling, 2 angular scalings.
    function vdi_bao_proptocomov3D(t,r)
      implicit none
      real(dl) :: vdi_bao_proptocomov3D
      real(dl), intent(in) :: t,r
      
      real(dl) :: a0,Yp0,az,Ypz

      call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
          r=r,t=t,&
          a=az,Rp=Ypz)
      call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
          r=VP%r0,t=VP%t0,&
          a=a0,Rp=Yp0)

!      vdi_bao_proptocomov3D = ((a0**2*Yp0)/(az**2*Ypz))**(1.d0/3.d0)
! Oops: comoving observer has scale factor (volume) = 1...
      vdi_bao_proptocomov3D = (1.d0/(az**2*Ypz))**(1.d0/3.d0)

    end function vdi_bao_proptocomov3D
    
    ! vdi_bao_proptocomovFLRWDrag calls vdi_bao_proptocomov3D
    ! for time and radius that correspond to zdrag,
    ! taking into account that this zdrag must come from CAMB
    ! which only knows about the effective FLRW observer. Hence,
    ! t(z) and r(z) are those of the FLRW observer.
    function vdi_bao_proptocomovFLRWDrag(zdrag)
      implicit none
      real(dl) :: vdi_bao_proptocomovFLRWDrag
      real(dl), intent(in) :: zdrag
      !
      real(dl) :: tr(2)
      
      if(zdrag > void_zmax_fixed .or. voiderror)then

        VoidError = .true.
        vdi_bao_proptocomovFLRWDrag = 1.d30
      
      else
      
        call asplint(VPNR_FLRW%ztr(1,:),VPNR_FLRW%ztr(2:3,:),zdrag,tr,VPNR_FLRW%ddztr(2:3,:))
        vdi_bao_proptocomovFLRWDrag = vdi_bao_proptocomov3D(tr(1),tr(2))

      end if

    end function vdi_bao_proptocomovFLRWDrag
    
    function timeFLRWDrag(zdrag)
      implicit none
      real(dl) :: timeFLRWDrag
      real(dl), intent(in) :: zdrag
      !
      real(dl) :: tr(2)
      
      if(zdrag > void_zmax_fixed .or. voiderror)then

        VoidError = .true.
        timeFLRWDrag = 1.d30
      
      else
      
        call asplint(VPNR_FLRW%ztr(1,:),VPNR_FLRW%ztr(2:3,:),zdrag,tr,VPNR_FLRW%ddztr(2:3,:))
        timeFLRWDrag = tr(1)

      end if

    end function timeFLRWDrag
    
    function vdi_bao_omegamofz(z)
      implicit none
      real(dl) :: vdi_bao_omegamofz
      real(dl), intent(in) :: z
      
      
    end function vdi_bao_omegamofz
    
