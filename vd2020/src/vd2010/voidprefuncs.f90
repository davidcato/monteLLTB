  
  subroutine sett0bg_Mt
		use lltb, only: oldlltb => lltb_normalization
    implicit none
    real(dl) :: x,x2, ub, u0(2)
!    integer :: i
    real(dl), parameter :: tau0 = (1.d0 / 6.d0 / pi)**(1.d0/6.d0)
!    logical, save :: warned = .false.

    call lltb_normalization(H0=VP%H0_paper,kb=VP%kb_om,kmax=VP%kmax,Lambda=VP%Lambda_paper,w=VP%w, & ! input
         Mtilde=VP%Mtilde_paper,Hlocal=VP%Hlocal,delta0=VP%delta0,t0=VP%t0bg,tbar=VP%tbar) ! output
    if(VoidDeltaMatterOnly)then
       if(abs(VP%w+1.d0)>1.d-10)stop 'cannot use DeltaMatterOnly with w/=-1.'
       call oldlltb(H0=VP%H0_paper,kb=VP%kb_om,kmax=VP%kmax,Lambda=VP%Lambda_paper,& !w=VP%w, & ! input
         Mtilde=VP%Mtilde_paper,Hlocal=VP%Hlocal,delta0=VP%delta0,t0=VP%t0bg)!,tbar=VP%tbar) ! output

    end if

    VP%Mtilde = VP%Mtilde_paper ! Back to unrescaled stuff.


    if(.false.) then
       !-----------------------------------------------------!
       
       ! Normalization such that a(L,t0)=1.d0:
       VP%Mtilde_paper=dsqrt(3.d0 * VP%H0**2/8.d0/pi/(1.d0+0.75d0 * VP%kb_om/pi))
       ! Pre 09/02/2010:
       ! VP%Mtilde=5.d-2/min(VP%zB,5.d-2) ! Just an 'empirical' choice, for numerical (mis-)behaviour.
       ! Actually, it seems to make no difference anymore.
       ! Or:
       ! VP%Mtilde=VP%Mtilde_paper
       ! Or: 
       VP%Mtilde = VP%Mtilde_paper / 10.d0 ! Seems to be fast with this choice..

!!$    if(voidtesting)then
!!$!       if(.not.warned)then
!!$       !write(*,*)'Danger: fixing Mtilde.'
!!$       call VoidBigAlert('Danger: fixing Mtilde.')
!!$!          warned = .true.
!!$!       end if
!!$       VP%Mtilde = 0.00006922622950819671
!!$    end if



       x=1.d0 + 1.5d0 * VP%kb_om / pi
       x2 = x

       if(x.eq.1.d0)then
          VP%t0bg = 2.d0/3.d0 / VP%H0 ! lim ub -> 0
       else if(dabs(x).ge.1.d0)then
          ub = acosh(x)
          VP%t0bg =  dabs(dsqrt(2.d0) / VP%H0 * dsqrt(dabs(1.d0+0.75d0*VP%kb_om/pi)) &
               * (dsinh(ub) - ub)&
               / dabs(dcosh(ub) - 1.d0)**(1.5d0))
       else if(dabs(x).lt.1.d0) then
          ub = acos(x)
          VP%t0bg = dabs(dsqrt(2.d0) / VP%H0 * dsqrt(dabs(1.d0+0.75d0*VP%kb_om/pi)) &
               * (dsin(ub) - ub)&
               / dabs(dcos(ub) - 1.d0)**(1.5d0))
       else
          write(*,*) 'weird value of x? sett0bg_Mt.'
          VoidError=.true.
          return
       end if

       !-----------------------------------------------------!

       if(voidfeedback.gt.1)write(*,*)'Mtilde:', VP%Mtilde


       !-----------------------------------------------------!

       if(dabs(x).gt.1.d0)then
          x2= (dsinh(ub) - ub)
       else if(dabs(x).lt.1.d0) then
          x2= - (dsin(ub) - ub)
       else
          ub=0.d0
          x2=0.d0
       end if

       if(VP%kb_om.eq.0.d0)then
          x=VP%kmax**(1.5d0) * dsqrt(3.d0)/2.d0/pi**(1.5d0)
          u0 = uofsinhu_u(x)

          VP%Hlocal = 3.d0 * VP%kmax**(3.d0/2.d0) * VP%Mtilde / dsqrt(2.d0)/pi * dsinh(u0(1))/(dcosh(u0(1))-1.d0)**2

          ! now we've got u0 and ub:
          VP%delta0 = (3.d0 / 2.d0 / pi * (VP%kmax + VP%kb_om) / (dcosh(u0(1)) - 1) )**3 - 1.d0

       else if (VP%kmax + VP%kb_om .eq. 0.d0) then

          ! lucky you..
          ! x = 0 -> q = 0 / u = 0
          VP%Hlocal = dsqrt(2.d0)/pi * VP%Mtilde * dabs(VP%kb_om)**(3.d0/2.d0)*x2! (dsinh(ub) - ub)

          VP%delta0 = 0.d0

       else if (VP%kmax+VP%kb_om .gt. 0.d0) then ! real

          x=dabs(dabs(VP%kmax/VP%kb_om + 1.d0)**(1.5d0)* x2)! x2 = (dsinh(ub) - ub)
          u0 = uofsinhu_u(x)
          VP%Hlocal = 3.d0 * dabs(VP%kmax + VP%kb_om)**(3.d0/2.d0) * VP%Mtilde / dsqrt(2.d0)/pi * dsinh(u0(1))/(dcosh(u0(1))-1.d0)**2


          ! now we've got u0 and ub:
          VP%delta0 = (3.d0 / 2.d0 / pi * (VP%kmax + VP%kb_om) / (dcosh(u0(1)) - 1) )**3 - 1.d0

       else if (VP%kmax+VP%kb_om .lt. 0.d0) then ! imaginary

          x=dabs(dabs(VP%kmax/VP%kb_om + 1.d0)**(1.5d0)*x2)! x2 = (dsinh(ub) - ub)
          u0 = qofsinq_q(x)
          VP%Hlocal = 3.d0 * dabs(VP%kmax + VP%kb_om)**(3.d0/2.d0) * VP%Mtilde / dsqrt(2.d0)/pi * dsin(u0(1))/(dcos(u0(1))-1.d0)**2

          ! now we've got u0 and ub:
          VP%delta0 = (3.d0 / 2.d0 / pi * (VP%kmax + VP%kb_om) / (dcos(u0(1)) - 1) )**3 - 1.d0

       else
          write(*,*) 'weird value of k? sett0bg_Mt.'
          VoidError=.true.
          return
       end if

    end if ! End skipping


    if(VP%Mtilde/=VP%Mtilde_paper)then
       VP%H0 = VP%H0_paper *  VP%Mtilde /  VP%Mtilde_paper
       VP%Lambda = VP%Lambda_paper *  (VP%Mtilde /  VP%Mtilde_paper)**2
    else
       VP%H0 = VP%H0_paper 
       VP%Lambda = VP%Lambda_paper 
    end if

    VP%Hlocal_paper = VP%Hlocal * VP%Mtilde_paper / VP%Mtilde

    VP%t0=VP%t0bg * VP%Mtilde_paper / VP%Mtilde

    VP%t0_Gyr = VP%t0bg * 977.8d0 ! Gyr

    !----------------------------------------------------------!
    ! Transformation such that delta0 satisfies -1 < delta0 < 1
    ! hence, delta0 = (rho_largest - rho_smallest)/rho_largest
    if(VoidDelta0Type == 0)then
       if(VP%delta0>0._dl)VP%delta0=1._dl - 1._dl / (VP%delta0 + 1)
    else
       write(*,*)'delta0 is normal.'
    end if
    !----------------------------------------------------------!

    if(vp%t0<0.)then ! then the bakcground cosmology is probably contracting as well. Sick.
       write(*,*)'Negative t0??'
       voiderror=.true.
       return
    end if


    if(voidfeedback.gt.2)write(*,*)'tau0/tau0_paper (curvature dependent):',(VP%t0bg*VP%Mtilde_paper)**(1.d0/3.d0)/tau0
    if(voidfeedback.gt.2)write(*,*)'t0:',VP%t0
    if(voidfeedback.gt.2)write(*,*)'t0bg:',VP%t0bg
    if(voidfeedback.gt.2)write(*,*)'t0 (Gyr):',VP%t0_Gyr

    if(voidfeedback.gt.2)write(*,*)'Hlocal:',VP%Hlocal
    if(voidfeedback.gt.2)write(*,*)'Hlocal_paper:',VP%Hlocal_paper
    if(voidfeedback.gt.2)write(*,*)'delta0:',VP%delta0

    if((voidfeedback.gt.1).and. voiderror)write(*,*)'VoidError in sett0bg_Mt'

!!$write(*,*)'Mt: ',testMt,VP%Mtilde_paper
!!$write(*,*)'Hl: ',testHl,VP%Hlocal_paper
!!$write(*,*)'d0: ',testd0,VP%delta0
!!$write(*,*)'t0: ',testt0,VP%t0bg

  end subroutine sett0bg_Mt

  subroutine VoidKmaxofCMBkmaxTRUE(inputkmax1)
    implicit none
!    real(dl) :: VoidKmaxofCMBkmax
    real(dl), intent(in) :: inputkmax1
    real(dl) :: inputkmax
    integer, parameter :: nmem = 2
    real(dl) :: k, ktop,kbot, d0, memf(nmem),memk(nmem)
    real(dl) :: f, df, dk, d0acc
    
    integer :: i, j
    integer, parameter :: NewtRaphStart = 4
    logical :: foundd0
!    integer :: thetime(3)
!    call system_clock(count=thetime(1),count_rate=thetime(3))
    real(dl) :: thistturn

    inputkmax = dble(inputkmax1)
    ! wv new: log prior on delta0>0
    if(inputkmax1>1.)then
       inputkmax = 1.d0 - 10.d0**(1-inputkmax1)
    end if

    if(inputkmax.lt.-1.d0)then
        call VoidErrorReport
        stop 'delta0 cannot be smaller than -1.'
     end if


!    if(dabs(VP%zB)<1.d-20)then
    if(VP%is_pure_flrw)then
        VP%kmax = 0.d0
       call sett0bg_Mt
       return
    end if

!!$    if(inputkmax.ge.0.d0)then
!!$        VP%kmax = inputkmax
!!$       call sett0bg_Mt
!!$       return
!!$    end if
    ! else it is delta0.

    ! Numerical inversion
    ktop = kmax_max
!    kbot = 0.d0
    kbot = -ktop
    k = (ktop + kbot) * 0.5d0

    VP%kmax = k
    call sett0bg_Mt
    if(VP%t0 < 0.)return ! Because then even the background is collapsing.
    d0 = VP%delta0
    
    foundd0 = .false.
    f=1.d0


    i=1
    do while (.not. foundd0)
       if(d0 .gt. inputkmax) then ! kmax was too small
          if(i.le. NewtRaphStart)then
             dk = (ktop - k)*0.5d0
              if (i==1) dk = dk/1.d1
          else
             dk = dabs(f / df)
             if(k+dk.ge. ktop)dk = (ktop - k)*0.5d0
          end if
          kbot = k
       else
          if(i.le. NewtRaphStart)then
             dk = (kbot - k)*0.5d0
              if (i==1) dk = dk/2.5d1
          else
             dk = -dabs(f / df)
             if(k+dk.le. kbot)dk = (kbot - k)*0.5d0
          end if
          ktop = k
       end if
       
       if(nmem.ge.2)then
          do j = nmem, 2, 1
             memf(j) = memf(j-1)
             memk(j) = memk(j-1)
          end do
       end if

       memf(1) = f
       memk(1) = k

       k = k + dk
       VP%kmax = k
       call sett0bg_Mt
       d0 = VP%delta0

       f = inputkmax - d0
!       if(i.ge. NewtRaphStart)df = (memf(1) - f) / (memk(1) - k)
       if(i.ge. NewtRaphStart)df = (memf(1) - f) / dk
       if(dabs(f/inputkmax).lt.u_precision)exit
       if(isnan(df))exit
       if(dabs(ktop-kbot)<1.d-30)exit
       !write(*,'(7ES15.4,I8)')d0,inputkmax,f,df,ktop,kbot,k,i
       i = i + 1
       if(i.gt.maxsteps)then
          if(voidfeedback.gt.0)write(*,*)'Problem in kmax(delta0). Hanging -> rejecting model.'
          voiderror=.true.
          return
       end if
!if(voidtesting)write(*,*)i
      if(inputkmax==0.d0)then
        d0=0.d0
        VP%kmax=0.d0
        exit
      end if
   end do
!if(voidtesting)write(*,*)'kmax found.'

   d0acc = dabs(d0/inputkmax-1.d0)
   if(inputkmax==0.d0)d0acc=0.d0
!    call system_clock(count=thetime(2))
!    write(*,*)'Delta0 precision:',dabs(d0/inputkmax-1.d0)
!    write(*,*)'Found in ',dble(thetime(2)-thetime(1))/dble(thetime(3)),'seconds.'
   if(d0acc.gt.u_precision)then
      !if(voidfeedback.gt.0)
      write(*,'(A)')' kmax(delta0) is too inaccurate. Are the input parameters physical? delta0 < delta0_min? Rejecting model.'
      voiderror=.true.
      return
   end if
      
    ! Now check that this universe has not already suffered a crunch somewhere:

   VP%L=1._dl
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=0._dl,t=1.d-10*VP%t0bg,tturn=vp%tturn)
    if(vp%t0bg.ge.2._dl*vp%tturn)then
       if(voidfeedback>1)write(*,*)'Centre is beyond crunch.'
       voiderror=.true.
    end if
    if(voidfeedback>1)write(*,*)'tturn_LTB:',vp%tturn


   VP%L=0._dl
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=0._dl,t=1.d-10*VP%t0bg,tturn=vp%bgtturn)
    if(vp%t0bg.ge.2._dl*vp%bgtturn)then
       if(voidfeedback>1)write(*,*)'Background is beyond crunch.'
       voiderror=.true.
    end if

    if(voidfeedback>1)write(*,*)'tturn_FLRW:',vp%bgtturn

    if(voiderror)then
       if(voidfeedback>1)write(*,*)'VoidError in VoidKmaxofCMBkmaxTRUE'       
       return
    end if

!    if(.not.voiderror)voiderror = VoidCheckKofr()

 end subroutine VoidKmaxofCMBkmaxTRUE


 subroutine VoidKmaxofCMBkmax(VCP)
!    use VDII_Params
   type(voidcmbparams) :: VCP
!   type(VoidParams) :: ThisVP
   real(dl) :: thiscmbkmax!, kmaxmem

   thiscmbkmax = VCP%delta0
   
   ! New: need to set dummy, because korf is fully checked in next routine
   VP%kmax_2 = 0.d0
   VP%L2_1 = 0.d0
   call VoidKmaxofCMBkmaxTRUE(thiscmbkmax)
   call sett0bg_Mt
   ! kmax has to be k(r) - k_b at r=0. So, set kmax_2 to present kmax, 
   ! continue to get kmax(d0_2).
   if(any(VoidProfile==(/3,7,8,9/)))then
      VP%kmax_2 = VP%kmax
      VP%L2_1 = VCP%L2_1
!      VP%d0_2 = VCP%d0_2
! New cosmomc10_Voids:
      VP%d0_2 = (ThisCMBkmax + 1.d0) * VCP%d0_2 - 1.d0
! end new
      call VoidKmaxofCMBkmaxTRUE(VP%d0_2)
      call sett0bg_Mt
      if(voidfeedback>2)then
         write(*,*)'kmax1:', VP%kmax
         write(*,*)'kmax2:', VP%kmax_2
         write(*,*)'d01:',VCP%delta0
         write(*,*)'d02:',VCP%d0_2
         write(*,*)'zB1:',VP%zB
         write(*,*)'L2/L1:',VP%L2_1
      end if
   end if
   
   ! kmax has to be k(r) - k_b at r=0. So, set kmax_2 to present kmax, 
   ! continue to get kmax(d0_2).
   if(any(VoidProfile==(/10,11/)))then
      VP%kmax_2 = VP%kmax
      VP%L2_1 = VCP%L2_1
!      VP%d0_2 = VCP%d0_2
! New cosmomc10_Voids:
      VP%d0_2 = (ThisCMBkmax + 1.d0) * VCP%d0_2 - 1.d0
! end new
      call VoidKmaxofCMBkmaxTRUE(VP%d0_2)
      call sett0bg_Mt
      VP%kmax_2 = VP%kmax_2 - VP%kmax
      if(voidfeedback>2)then
         write(*,*)'kmax1:', VP%kmax
         write(*,*)'kmax2:', VP%kmax_2
         write(*,*)'d01:',VCP%delta0
         write(*,*)'d02:',VCP%d0_2
         write(*,*)'zB1:',VP%zB
         write(*,*)'L2/L1:',VP%L2_1
      end if
   end if
   
   ! kmax has to be k(r) - k_b at r=0. So, set kmax_2 to present kmax, 
   ! continue to get kmax(d0_2).
   if(any(VoidProfile==(/15/)))then
      VP%kmax_2 = VP%kmax
      VP%L2_1 = VCP%L2_1
!      VP%d0_2 = VCP%d0_2
! New cosmomc10_Voids:
      VP%d0_2 = VCP%d0_2
! end new
      call VoidKmaxofCMBkmaxTRUE(VP%d0_2)
      call sett0bg_Mt
      VP%kmax_2 = VP%kmax_2 - VP%kmax
      if(voidfeedback>2)then
         write(*,*)'kmax1:', VP%kmax
         write(*,*)'kmax2:', VP%kmax_2
         write(*,*)'d01:',VCP%delta0
         write(*,*)'d02:',VCP%d0_2
         write(*,*)'zB1:',VP%zB
         write(*,*)'L2/L1:',VP%L2_1
      end if
      if(vp%zb<0)setlboost=1.d-3
   end if
   
   if(any(VoidProfile==(/16/)))then
      VP%kmax_2 = 0._dl
      VP%L2_1 = VCP%L2_1
      VP%d0_2 = VCP%d0_2
      if(VP%d0_2<1._dl)stop 'd0_2 must be > 1 for profile 16.'
      if(vp%zb<0)setlboost=1.d-2
   end if

   if(any(VoidProfile==(/17/)))then
      ! With this one, k(r=0)=kb
      VP%kmax_2 = VP%kmax
      VP%kmax=0._dl
!      call VoidKmaxofCMBkmaxTRUE(VP%d0_2)
      call sett0bg_Mt
      VP%L2_1 = VCP%L2_1
      VP%d0_2 = VCP%d0_2
      if(VP%d0_2<1._dl)stop 'd0_2 must be > 1 for profile 17.'
      if(vp%zb<0)setlboost=1.d-2
   end if

   if(any(VoidProfile==(/18,19/)))then
      VP%kmax_2 = 0._dl
      VP%L2_1 = VCP%L2_1
      VP%d0_2 = VCP%d0_2
      if((VP%L2_1<0._dl).or.(VP%L2_1>1._dl))stop 'L1_2 must be 0 < L1_2 < 1 for profile 18.'
      if(vp%zb<0)setlboost=1.d-2
   end if

   if(.not.any(VoidProfile==(/1,2,3,4,5,6,12,13,14,18/)))then
      if(VoidRObs<0)stop 'RObs as MC variable is only possible for selected profiles.'
   end if
   if(VoidRObs<0)then
      VP%d0_2 = VCP%d0_2 ! re-use unused parameter.
   end if

    ! Enforce collapse here:
    if(VoidForceCollapse)then
       if(VP%t0<vp%tturn)then
          call VoidBigAlert('You chose to force collapse at r=0. Starting to reject models now.')
          voiderror=.true.
       end if
    end if
    
    if(voidfeedback.gt.1)write(*,*)'t0 (Gyr):',VP%t0_Gyr, VP%t0
    if(voidfeedback.gt.1)write(*,*)'delta0:',VP%delta0

    
    if(voiderror)return

 end subroutine VoidKmaxofCMBkmax

! VoidCheckKofr test that 1+2r^2k(r)Mt^2 > 0 for 0<r<L.
 function VoidCheckKofr()
   logical :: VoidCheckKofr


  VoidCheckKofr = VoidCheckPositiveDefiniteness(0.d0,VP%L,myf,mydf)

  contains

  function mydf(r)
    implicit none
    real(dl) :: mydf
    real(Dl), intent(in) :: r
    
    mydf = 2*r*VoidKofr(r) + r**2*VoiddKdr(r)
    
  end function mydf

  function myf(r)
    implicit none
    real(dl) :: myf
    real(Dl), intent(in) :: r
    
    myf = 1._dl + 2._dl * r**2 * VoidKofr(r) * VP%Mtilde**2
    
  end function myf

 end function VoidCheckKofr
 
 ! VoidCheckPositiveDefiniteness checks that func(x) is positive definite on range
 ! xmin < x < xmax, returns false (!) if it is positive definite. So, true if there is an error.
 ! Function dfunc(x) is f'(x), is optional, speeds up process.
! after a coarse probing, it finds the exact minimum, to be sure we didn't miss it. 
 function VoidCheckPositiveDefiniteness(xmin,xmax,func,dfunc)
   interface
      function afunc(x)
        use precision
        implicit none
        real(dl) :: afunc
        real(dl), intent(in) :: x
      end function afunc
   end interface
   procedure(afunc) :: func
   procedure(afunc), optional :: dfunc
   real(dl), intent(in) :: xmin,xmax ! range over which to probe func
   logical :: VoidCheckPositiveDefiniteness
   real(dl) :: r, Eofr, minE, minr
   real(dl) :: rt, rb, dkt, dkb
   real(dl) :: fmem(3,2) ! steps times ( x, f(x) )
   integer :: i, imax

   imax = 100!0!1000
  
   do i = 1, imax
      r = xmin + dble(i-1)/dble(imax-1)*(xmax-xmin)
!      Eofr = 1._dl + 2._dl * r**2 * VoidKofr(r) * VP%Mtilde**2
      Eofr = func(r)
      VoidCheckPositiveDefiniteness = 0._dl >= Eofr

      if(i==1)then
        minE = Eofr
        minr = r
      end if
      if(Eofr<minE)then
        minE=Eofr
        minr = r
      end if

!        write(123,'(50ES14.5)')r/VP%L,Eofr,myf()
      if(VoidCheckPositiveDefiniteness)return
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
        VoidCheckPositiveDefiniteness = 0._dl >= Eofr
        if(VoidCheckPositiveDefiniteness)return
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
      myf = (fmem(3,2)-fmem(1,2))/(fmem(3,1)-fmem(1,1))
    end if
    
  end function myf
   
 end function VoidCheckPositiveDefiniteness
