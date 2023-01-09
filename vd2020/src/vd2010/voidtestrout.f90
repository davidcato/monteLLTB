   subroutine PlotKofr()
     implicit none
     integer :: i, imax
     real(dl) :: r
     open(1230,file='test_kofr.txt')
     write(0,*)'Writing this profile to kofr.txt.'
     write(*,*)'1: r/L, 2: k(r/L), 3: dkdr(r/L), 4: Numerical dkdr'
     write(1230,*)'# 1: r/L, 2: k(r/L), 3: dkdr(r/L), 4: Numerical dkdr'

     imax = 1000
     do i = 1, imax
        r=dble(i-1)/dble(imax-1)*VP%L
        write(1230,'(20ES30.16E5)')r,voidkofr(r),voiddkdr(r),(voidkofr(r)-voidkofr(r*0.999d0))/(r*0.001d0)
     end do
     close(1230)

     

   end subroutine PlotKofr


   subroutine PlotDensity()
     use wlltb
     implicit none
     integer :: i, imax, j, jmax
     real(dl) :: r,t, u(3), rho, mR, mRp, mS, rmin, rmax, thisk, thisint(2), lastr, dr, rhout!, mte
     real(dl) :: rg,tg, ug(3), rhog, mRg, mRpg, mSg, RF, RpF, rhorat
     real(dl) :: nextr, a_flrw, r_flrw, tstart, tturn
     real(dl) :: testap, testrp, testa, try_ap, tests, testsd, mSdot, testh, testrpd, testapd
     real(dl) :: testhp, mRpd, numHp,r_flrw_new, aF, ComovLMpc, ComocRobsMPC, RhoAtRobs, ma, rhout0
     character (len = 1024) :: mstring
     real(dl) :: rmpc, testhofz, garadius, For3,For3out
     
     write(*,*)'VoidError:',voiderror
     voiderror = .false.
     write(*,*)'kmax:',VP%kmax
     write(*,*)'Lbar:',VP%L*VP%Mtilde
     write(*,*)'Writing profile functions to file.'
!!$     write(*,*)'H0_out:',VP%H0
!!$     write(*,*)'H0_in:',VP%Hlocal_paper
!!$     write(*,*)'(1+z_FLRW) / (1+z_in):',(VP%H0/VP%Hlocal_paper)**(2.d0/3.d0)
!!$     mte=VoidTVoverTF()
!!$!     write(*,*)'VoidTVoverTF:',VoidTVoverTF()

     imax = 1000!00
!     imax = 1000
     jmax = 1!2!75!1!8!75
     rmin = max(0.d0+1.d-10,VP%r0)
     
     garadius = VoidGARadius()

     ! You can try to set rmin=0, but you will crash:
     ! the geodesics integration starts as r=r0, and the following
     ! loop will ask for t_geodesic(r). So if rmin < r0, this t will  
     ! be an extrapolated spline, going to crazy values very fast.
     if(abs(VP%r0)>abs(rmin))rmin=VP%r0

!     if (VoidTesting)then
!        rmin = 0.d0+1.d-2*VP%L
!     end if

     if(voidprofile==2)then
        rmax = 5.d0*VP%L
     else
        rmax = 1.5d0*VP%L
     end if
!rmax = nr(z=void_zmax)
!     rmax = 10.d0*VP%alpha
!     rmin = VP%L*0.99d0
!     rmax = VP%L*1.01d0
     t = VP%t0
     thisint = 0.d0
     lastr = rmin
     if(.not.Do_bestfit_output)open(222,file='test_LTB_density.txt')
     mstring = '# 1: r/L, 2: rho(t0), 3: R(t0), 4: Rprime(t0), 5: k(r), 6: z_ltb, 7: z_flrw, 8: S. Along geodesic: 9: t, 10: rho, 11: R, 12: Rp, 13: S, 14 rho/rhoFLRW(r,t)... 17: rho(t)/rhoout(t), 18:r_flrw = Rltb/a_FLRW(t)/L, 19: t(Gyr), 20: tturn - t, 21: testrp, 22: testa, 23: testap, 24: r, 25: num_ap, 26: tests. 27: testsdot, 28:Sdot, 29: adot, 30: num_adotprime, 31: Rdotprime, 32: testHp, 33: derived HP, 34: numHp, 35: apdot, 36: r_flrw_new, 37: k''(r), 38: N(k''(r)), 39: DA(z_ltb), 40: r^2 k(r) Mtilde^2, 41: R_obs (as col 36), 42: rho(robs)/rho_flrw, 43: rho(r,t) / rho_flrw(t0), 44: r(Mpc), 45: H(r,t0), 46: H(z), 47: L(Mpc)'
     write(222,'(A)')trim(mstring)
     write(*,'(A)')trim(mstring)     

     ma = voida(rmax,VP%t0)
     mRp = avoidRp(rmax,VP%t0)
     rhout0 = VP%Mtilde**2 / mRp / ma**2

!tstart = 0.1_dl!1._dl!0.1_dl !0.3_dl
tstart = NT(z=1100.d0)
if(jmax==1)tstart=VP%t0
do j = 1,jmax

     thisint = 0.d0

!t=VP%t0*(tstart+(1._dl-tstart)*(dble(j)/dble(jmax)))
if(jmax/=1)then
  t=tstart + (VP%t0-tstart)*dble(j-1)/dble(jmax-1)
else
  t=tstart
end if
ComovLMpc=VoidComovDistMPC(0.d0,VP%L,t,VP%r0)
ComocRobsMPC=VoidComovDistMPC(0.d0,VP%r0,t,VP%r0)
mR = voidR(VP%r0,t)
ma = voida(VP%r0,t)
mRp = avoidRp(VP%r0,t)
mRpd = avoidRpdot(VP%r0,t)
mS = avoidS(VP%r0,t)
mSdot = avoidSdot(VP%r0,t)
RhoAtRobs = VP%Mtilde**2 / mRp / ma**2

     u = voidu(rmax,t)
     mR = voidR(rmax,t,u)
     mRp = avoidRp(rmax,t,u)
     rhout = VP%Mtilde**2 * rmax**2 / mRp / mR**2
     ! r_flrw_new = rmin
     r_flrw_new = VoidComovDistMPC(0.d0,rmin, t, VP%r0)/ComovLMpc
        
     do i = 1, imax

        r = rmin + (rmax-rmin)* dble(i-1)/dble(imax-1) 
        nextr =  rmin + (rmax-rmin)* dble(i)/dble(imax-1) 
        if(r<VP%L .and. nextr>VP%L)r=VP%L

!        u = voidu(r,t)
        mR = voidR(r,t)
        ma = voida(r,t)
        mRp = avoidRp(r,t)
        mRpd = avoidRpdot(r,t)
        mS = avoidS(r,t)
        mSdot = avoidSdot(r,t)
        rho = VP%Mtilde**2  / mRp / ma**2

        call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
             dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=r,t=t,tturn=tturn)


        call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
             dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=r,t=t,a=testa,ap=testap,apd=testapd,s=tests,sd=testsd, h=testh,rpd=testrpd,hp=testhp,For3=For3)
testrp=testa+r*testap
        call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
             dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=2*VP%L,t=t,For3=For3out)



        call VoidDistFuncs(r=r,t=t,a=a_flrw,obs='FLRW')
        r_flrw = mR / a_flrw / VP%L
        
        r_flrw_new = r_flrw_new + (rmax-rmin)/dble(imax-1) * mS / a_flrw /VP%L
if(ms<0.d0)call VoidBigAlert('Warning: test routine has S(r,t) < 0. Unphysical.')
        rg = r
        if(rg<VP%r0)stop 'You are asking for a radius which is not inside the range of geodesics integration. Stop plotdensity, voidtestrout.f90'
        tg = nt(r=rg)
        if(isnan(tg))stop 'isnan tg=vp%t0'


!        ug = voidu(rg,tg)
        mRg = voidR(rg,tg)
        mRpg = avoidRp(rg,tg)
        mSg = avoidS(rg,tg)
        rhog = VP%Mtilde**2 * rg**2 / mRpg / mRg**2
        testhofz=voidH(rg,tg)


        call VoidDistFuncs(r=rg,t=tg,Rltb=RF,Rp=RpF,obs='FLRW')
        rhorat =  (RpF * RF**2) / (mRpg * mRg**2)  

        thisk = voidkofr(r)
        dr = r - lastr 
        if(r>0.d0 .and. r<GARadius)thisint = thisint + dr * (/rho,rhout/) * mR**2 * mS * 4.d0 * pi
        lastr = r
numhp=(voidH(r,t)-voidH(r*0.999d0,t))
numhp=numhp/(r*0.001d0)
!     (voidH(r,t)-voidH(r*0.999d0,t))/(r*0.001d0)

        rMpc=r_flrw_new * ComovLMpc

        write(mstring,'(18ES21.10E5,F10.3,28ES21.10E5)')r/VP%L,rho,mR,mRp,thisk,NZ(r),z_FLRW(r), mS, tg, rhog,mRg,mRpg, mSg, rhorat, (voidkofr(r)-voidkofr(r*0.999d0))/(r*0.001d0),voiddkdr(r),rho/rhout, r_flrw, t*977.8, tturn-t,testrp, testa, testap, r,(voida(r,t)-voida(r*0.999d0,t))/(r*0.001d0),tests, testsd, mSdot, testa*testh,(voida(r,t)*voidH(r,t)-voida(r*0.999d0,t)*voidH(r*0.999d0,t))/(r*0.001d0), testrpd, testhp, (mRpd-mRp*voidH(r,t))/mR, numHp, testapd, r_flrw_new, voiddkdr(r),(voidkofr(r)-voidkofr(r*0.999d0))/(r*0.001d0), DA(0.d0,(/tg,rg/)),r**2 * VP%Mtilde**2 * thisk,ComocRobsMPC/ComovLMpc,RhoAtRobs/Rhout, rho/rhout0, rMpc, testh, testhofz, ComovLMpc
!        if(modulo(i,imax/10)==1)write(*,'(A)')trim(mstring)
        if(modulo(i,imax/10)==1)then
           ! This may be intel fortran only: '(A,$)' for noadvance
           write(6,'(A)', advance="no")'.'
        end if
        write(222,'(A)')trim(mstring)
        if(voiderror)then
           write(*,*)'Exiting with error. R too large? S imaginary?'
!           exit
           voiderror=.false.
        end if

     end do
write(222,*)
!     write(*,*)
!     write(*,*)'t/tlss,t/t0:',t/NT(z=1100.d0),t/VP%t0
!     write(*,*)'Integral:',thisint
!     write(*,*)'\delta M / M:',thisint(1)/thisint(2)-1.d0

!     write(*,*)
end do
     if(.not.Do_bestfit_output)close(125)

     
     !stop 'density written.'
        
     
   end subroutine PlotDensity
   
   subroutine CompareWithMath
     implicit none
     integer, parameter :: fnum = 125
     !     character(len=128), parameter :: infile = "Math_LTB_largek.tsv", outfile = "voiddiffs.txt"
     character(len=128), parameter :: infile = "Math_FLRW_largek.tsv", outfile = "voiddiffsF.txt"
     character(len=128) :: dummy
     real(dl), allocatable :: mathfile(:,:)
!     real(dl), allocatable :: oneline(:)
     real(dl) :: myS, mySdot, hisS,hisSdot
!     integer :: nitems, nlines, i, j, mystat, lengthwidth(2)
     integer :: nitems, nlines, i, mystat, lengthwidth(2)
     integer :: z, r, t, S, Sdot
     
     VP%kmax = 15.d0
     VP%H0 = 72.d0
     call sett0bg_Mt
     VP%L = 0.d0!14.9365
     call sett0bg_Mt
     
     !write(*,*)'u0:',voidu(0.d0,VP%t0)
     !write(*,*)'t0:',VP%t0
     !write(*,*)'Mtilde:',VP%Mtilde
     !stop
     !0.005 / VP%Mtilde !3539.1209623946206
     
     lengthwidth = GetLengthWidth()
     nlines = lengthwidth(1)
     nitems = lengthwidth(2)
     
     allocate(mathfile(nlines,nitems))
     
     open(fnum,file=infile)
     
     do i = 1,nlines
        read(fnum,'(A)')dummy
        read(dummy,*)mathfile(i,1:nitems) ! weird, this handles tabs properly....
     end do
     
     close(fnum)
     
     
     ! assuming items are
     ! z  r  t  S  Sdot
     z = 1
     r = 2
     t = 3
     S = 4
     Sdot = 5
     
     
     open(fnum,file=outfile, status="REPLACE")
     
     
     do i = 1,nlines
        myS = avoidS(mathfile(i,r),mathfile(i,t))
        mySdot = avoidSdot(mathfile(i,r),mathfile(i,t))
        hisS = mathfile(i,S)
        hisSdot = mathfile(i,Sdot)
        write(dummy,'(8ES16.7)')mathfile(i,r),mathfile(i,t),myS,hisS,myS-hisS,mySdot,hisSdot,mySdot-hisSdot
        write(fnum,'(A)')dummy
        !        write(*,'(A)')dummy
        
     end do
     
     close(fnum)
     
     
     deallocate(mathfile)
     
     stop 'Compared.'
     
     return
   contains
     
     function GetLengthWidth()
       implicit none
       integer getlengthwidth(2), gli
       open(fnum,file=infile)
       mystat = 0
       gli=0
       do while (mystat ==0)
          read(fnum,'(A)',iostat=mystat)dummy
          if(mystat.ne.0)exit
          gli=gli+1
          if(gli.eq.1) GetLengthWidth(2)=CountItemsIn(trim(adjustl(dummy)))
       end do
       close(fnum)
       GetLengthWidth(1) = gli
     end function GetLengthWidth
     
     
     
   end subroutine CompareWithMath


   subroutine VoidDensAtCMB()
     implicit none
     real(dl) :: t,r, Lmem, altb, aflrw, delta
     
     t=nt(z=1100.d0)
     r=nr(1100.d0)
     Lmem = VP%L
     ! rho = VP%Mtilde**2 * r**2 / mRp / mR**2 = VP%Mtilde**2 / a(r,t)^3

     call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
             dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=0.d0,t=t,a=altb)

     VP%L=0.d0
     call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
             dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=0.d0,t=t,a=aflrw)

     delta = (aflrw/altb)**3 - 1.d0

     write(*,*)'(rho-rhobar)/rhobar at t_dec:',delta



     
     

   end subroutine VoidDensAtCMB

    subroutine VoidWOFZ()
      implicit none
      integer, parameter :: nn=200!1000
      real(dl) :: z(nn),w(nn)
      integer :: i, newn
      write(*,*)'hello VoidWOfZ. VoidError:',VoidError

      newn = min(size(VPNR%ztr(1,:)),nn)

      do i = 1, newn
!        if(VPNR%ztr(1,i)>2*VP%zb)exit
        if(VPNR%ztr(1,i)>max(2.d0,2*VP%zb))exit
      end do
      newn = min(newn,i)

      w=VP%w
      z=0.d0
      call VoidMakeWofZ(w(1:newn),z(1:newn),2*VP%zb,VPNR%ztr(:,1:newn))

      open(1234,file="Void_wofz.txt",status="REPLACE")

      do i = 1, newn!nn

        write(1234,'(10ES14.5)')z(i),w(i)+1.d0, VP%w,VP%zb, VP%w

      end do

      close(1234)

    end subroutine VoidWOFZ

!
! eq (3) from http://arxiv.org/pdf/astro-ph/0702670.pdf
!
    subroutine VoidMakeWOFZ(w_inout,z_inout,zmax,ztr)
      use vdii_tools
      implicit none
      real(dl), intent(out) :: z_inout(:), w_inout(:)
      real(dl), parameter :: zmin=0.01d0
      real(dl), intent(in) :: zmax!=2.d0
      real(dl), intent(in), optional :: ztr(:,:)
!      integer, parameter :: imax = 1000
      integer :: imax
      real(dl) :: z
      integer :: i
      real(dl) :: thederiv(2), omm, omk, theW, theH0, thisHz, realH
      real(dl) :: Rltb, Rp, Rd, Rpd, Sd, a,H, Rpp, Rpdd, thisk, thisdkdr, M2, Rppd, Rdd, add
      real(dl) :: r, t, numW, Dlum, myW, flrwdaz, flrwzofr, prevr, prevz
      character(len=1024) :: outstr
      type(FakeCMBParams) :: CMBInside
    real(dl), external :: rombint

      ! vdii_D1D2 returns the first and second derivative of f(x) in x,
      ! accurate up to order h^4.
      !  subroutine vdii_D1D2(g,x,d, h)

      if(VoidError)then
        z_inout = 0.d0
        w_inout = VP%w
        return
      end if
!      zmax = 2*VP%zb
      imax = size(z_inout)
      if(size(w_inout)/=imax)then
        write(0,*)"Sizes of z(:) and w(:) do not match in MakeVoidWofZ."
        stop
      end if
      if(present(ztr))then
        if(size(ztr(1,:))/=imax)then
          write(0,*)"Sizes of z(:) and ztr(1,:) do not match in MakeVoidWofZ."
          stop
        end if
        if(size(ztr(:,1))<3)then
          write(0,*)"Size ztr(:,1) is less than 3 in MakeVoidWofZ:",size(ztr(:,1))
          stop
        end if
      end if

      M2 = VP%Mtilde**2
      if(VoidWofZObsParams==0)then

        omm = 1 - VP%omv - VP%omk
        omk = VP%omk
        theH0 = vdii_H0avv(0.023d0,0.1d0)/3d5 ! Riess' HST range, in 1/Mpc

        theH0 = VP%H0
        
        omm= omm * VP%H0**2 / theH0**2
        omk= omk * VP%H0**2 / theH0**2

      else if(VoidWofZObsParams==1)then
        call VoidMakeCMBin(VP%CMBParams,CMBInside)

        theH0 = CMBInside%H0
        omm = (CMBInside%omdmh2 + CMBInside%ombh2) / (1.d-4*theH0**2)
        omk = CMBInside%omk

      else
        write(0,*)'Do not know with this value of VoidWofZObsParams:',VoidWofZObsParams
        stop 
      end if

!      write(*,*)'theH0 here:',theH0
!      write(*,*)'omm here:',omm
!      write(*,*)'omv here:',VP%omv
!      write(*,*)'omk here:',VP%omk,vp%kmax
      r=0
      flrwzofr=0
      flrwdaz=0
      z=0
      do i = 1, imax
        prevr = r
        prevz=z
!write(*,'(I5)', advance="no")i
        if(present(ztr))then
          z = ztr(1,i)
          t = ztr(2,i)
          r = ztr(3,i)
        else
          z=zmin+dble(i-1)/dble(imax-1)*(zmax-zmin)
          r=NR(z)
          t=NT(z)
        end if
        z_inout(i) = z
!write(*,*)z,t,r
        call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
                   dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=r,t=t, &
                   Rltb=Rltb,a=a,add=add,H=H,Rp=Rp,Rpd=Rpd,Sd=Sd,Rpdd=Rpdd)
        Rd = a*H*r
        Rdd = r*add
        thisk = VoidKofr(r)
        thisdkdr = VoiddKdr(r)

        call vdii_D1D2(myRpfunc,r,thederiv,min(1.d-4,VP%L*1.d-1))
        Rpp = thederiv(1)
        call vdii_D1D2(myRpdfunc,r,thederiv,min(1.d-4,VP%L*1.d-1))
        Rppd = thederiv(1)

!        Rltb = Rltb * theH0
!        Rp = Rp * theH0
!        Rpd = Rpd * theH0
!        Rpp = Rpp * theH0
!        Rppd = Rppd * theH0
!        Rpdd = Rpdd * theH0
!        Rd = Rd * theH0
!        Rdd = Rdd * theH0


        theW = ((1 + z)*(2*(Rdd*Rp**2*Rpd - 2*Rd*Rp*Rpd**2 + 2*Rltb*Rpd**3 - &
              Rd*Rp**2*Rpdd + Rpd*Rpp - Rp*Rppd + M2*r**2*Rp*Rpd*thisdkdr + &
              2*M2*r*(r*Rpd*Rpp + Rp*(Rpd - r*Rppd))*thisk + Rp*Rpd**2*Sqrt(1 + &
              2*M2*r**2*thisk) + Rp**2*Rpdd*Sqrt(1 + 2*M2*r**2*thisk) - &
              Rd*Rpd*Rpp*Sqrt(1 + 2*M2*r**2*thisk) + Rd*Rp*Rppd*Sqrt(1 + &
              2*M2*r**2*thisk))*(1 + omk*Rltb**2*theH0**2*(1 + z)**2) - &
              Rpd**2*(Rltb*Rpd + Rp*(-Rd + Sqrt(1 + 2*M2*r**2*thisk)))*(1 + &
              (omk*theH0**2*(2*Rltb*Rpd + Rp*(-Rd + Sqrt(1 + &
              2*M2*r**2*thisk)))**2*(1 + z)**2)/Rpd**2)))/(3.*Rpd**2*(Rltb*Rpd + &
              Rp*(-Rd + Sqrt(1 + 2*M2*r**2*thisk)))*(-1 - z + &
              omm*Rltb**2*theH0**2*(1 + z)**4 - (2*Rltb*theH0**2*(2*Rltb*Rpd + &
              Rp*(-Rd + Sqrt(1 + 2*M2*r**2*thisk)))*(1 + z)**3*(omk + omm + &
              omm*z))/Rpd + (theH0**2*(2*Rltb*Rpd + Rp*(-Rd + Sqrt(1 + &
              2*M2*r**2*thisk)))**2*(1 + z)**3*(omk + omm + omm*z))/Rpd**2))

        myW =-1 + (-2*omk*(1 + z) - 3*omm*(1 + z)**2 + (2*Sqrt((Rpd**2*(1 + &
              omk*Rltb**2*theH0**2*(1 + z)**2))/(theH0**2*(-(Rd*Rp) + Rltb*Rpd + &
              Rp*Sqrt(1 + 2*M2*r**2*thisk))**2))*(omk*Rltb*theH0*(1 + z) - &
              ((Rdd*Rp**2*Rpd + Rpd*Rpp - Rp*Rppd + M2*r**2*Rp*Rpd*thisdkdr + &
              2*M2*r*Rp*Rpd*thisk + 2*M2*r**2*Rpd*Rpp*thisk - &
              2*M2*r**2*Rp*Rppd*thisk - Rp*Rpd**2*Sqrt(1 + 2*M2*r**2*thisk) + &
              Rp**2*Rpdd*Sqrt(1 + 2*M2*r**2*thisk) - Rd*(Rp**2*Rpdd + &
              Rpd*Rpp*Sqrt(1 + 2*M2*r**2*thisk) - Rp*Rppd*Sqrt(1 + &
              2*M2*r**2*thisk)))*(1 + omk*Rltb**2*theH0**2*(1 + &
              z)**2))/(Rpd*theH0*(-(Rd*Rp) + Rltb*Rpd + Rp*Sqrt(1 + &
              2*M2*r**2*thisk))**2*(1 + z))))/Sqrt(1 + omk*Rltb**2*theH0**2*(1 + &
              z)**2))/(-3*omk*(1 + z) - 3*omm*(1 + z)**2 + (3*Rpd**2*(1 + &
              omk*Rltb**2*theH0**2*(1 + z)**2))/(theH0**2*(-(Rd*Rp) + Rltb*Rpd + &
              Rp*Sqrt(1 + 2*M2*r**2*thisk))**2*(1 + z)))

!        call vdii_D1D2(Dlumfunc,z,thederiv,1.d-2)
!        Dlum = theH0*DA(z) * (1+z)**2
!        numW = (2*(1 + z)*(-((-Dlum + (1 + z)*thederiv(1))*(1 + &
!            omk*thederiv(1)**2))/2. + (Dlum**2*omk + (1 + &
!            z)**2)*thederiv(2)))/(3.*(-Dlum + (1 + z)*thederiv(1))*(-1 + &
!            Dlum**2*omm - z - 2*Dlum*(omk + omm*(1 + z))*thederiv(1) + (1 + &
!            z)*(omk + omm*(1 + z))*thederiv(1)**2))

!        flrwdaz = flrwdaz + rombint(ddadz,prevz,z_inout(i),1.d-8)
!        flrwzofr = flrwzofr + rombint(flrwhz,prevr,r,1.d-12)
!        write(outstr,'(10ES14.5)')z,Rltb,flrwdaz / (1+z),flrwzofr,VP%r0,DA_unnorm(z)
!        write(*,'(A)')trim(outstr)
!        write(1234,'(A)')trim(outstr)

        w_inout(i)=theW
      end do
!write(*,*)'check output 1234'

    contains
    
      function myRpfunc(rr)
        implicit none
        real(dl) :: myRpfunc
        real(dl), intent(in) :: rr

        call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
           dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=abs(rr),t=t, &
           Rp=myRpfunc)
        
      end function myRpfunc

      function myRpdfunc(rr)
        implicit none
        real(dl) :: myRpdfunc
        real(dl), intent(in) :: rr
        
        call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
           dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=abs(rr),t=t, &
           Rpd=myRpdfunc)
        
      end function myRpdfunc
      
      function Dlumfunc(zzz)
        implicit none
        real(dl) :: Dlumfunc
        real(dl), intent(in) :: zzz
        
        Dlumfunc = theH0*DA(zzz) * (1+zzz)**2

      end function Dlumfunc

      function flrwhz(zzz)
        implicit none
        real(dL) :: flrwhz, zzz

        flrwhz = theH0 * ( (1-VP%omv-VP%omk)*(1+zzz)**3 + VP%omk*(1+zzz)**2 + VP%omv*(1+zzz)**(3*(1+VP%w)))**0.5d0

      end function flrwhz

      function ddadz(zzz)
        implicit none
        real(dL) :: ddadz, zzz, omv

        ddadz = 1.d0/flrwhz(zzz)
      end function ddadz



    end subroutine VoidMakeWOFZ


    subroutine VoidWofZDispersion(wdisp)
      implicit none
      real(dl), intent(out) :: wdisp(:)
      integer :: nn
      integer :: ni
      integer :: i, binpoints, j, next
      integer, allocatable :: bincount(:)
      real(dl), allocatable :: wz(:),zz(:)
      real(dl) :: binsize, wavv(size(wdisp))

      if(.not.VP%defit_set)then
!        write(0,*)"You need to call vdi_GetIgnObsDEPars before calling VoidWofZDispersion."
        wdisp = 0.d0 ! if you don't de DEFit.
      end if

      binpoints=10
      binsize = 0.1

      nn = size(wdisp)
      if(nn<1)then
        write(0,*)'Size of wdisp is wrong. Smaller than one:',nn
        stop
      end if

!      ni=binpoints*nn
      do i = 1, size(VPNR%ztr(1,:))
        if(VPNR%ztr(1,i)>nn*binsize)exit
      end do

      ni=min(i,size(VPNR%ztr(1,:)))

      allocate(wz(ni))
      allocate(zz(ni))
      allocate(bincount(nn))
      wz=VP%w
      zz=0.d0

      call VoidMakeWOFZ(wz,zz,VPNR%ztr(1,ni),VPNR%ztr(1:3,1:ni))

      wdisp=0
      wavv=0
      next=1
      bincount = 0
      do j = 1, nn
!        do i = (j-1)*binpoints+1,j*binpoints
!          wavv(j) = wavv(j) + wz(i)
!        end do
        do i = next, ni
          if(zz(i)>j*binsize)exit
          wavv(j)=wavv(j)+wz(i)
          bincount(j)=bincount(j)+1
        end do
        next=i
        if(bincount(j)/=0)then
          wdisp(j) = wavv(j)/dble(bincount(j))
        else
          if(VoidIncludeFLRWDEFit/=0)write(0,*)"We have an empty bin in VoidWofZDispersion. Set VoidIncludeFLRWDEFit /= 0 and VoidIncludeFLRWDEFitzmax equal to max[zb] in params.ini."
          wdisp=-10
          exit
!          call VoidErrorReport()
!          stop
        end if
      end do

      deallocate(wz)
      deallocate(zz)
      deallocate(bincount)

    end subroutine VoidWofZDispersion


    subroutine testVoidTopHatGaussianField()
      implicit none
      real(dl) :: GARadius,t,For3,a,For3out, aout
  
      t=NT(z=1100.d0)
      garadius = VoidGARadius(t)
      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=2*VP%L,t=t,For3=For3out,a=aout)
        

      write(*,*)'at z=1100:'
      write(*,*)'Comoving Radius:',VoidComovDistMPC(0.d0,GARadius,VP%t0,0.d0)
      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=0.d0,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=0:',For3 * aout**3/ (For3out * a**3  )-1

      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=GARadius,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=r_GA:',For3 * aout**3/ (For3out * a**3  )-1

      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=VP%L,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=L:',For3 * aout**3/ (For3out * a**3  )-1
  
      t=NT(z=200.d0)
      garadius = VoidGARadius(t)
      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=2*VP%L,t=t,For3=For3out,a=aout)
        

      write(*,*)'at z=200:'
      write(*,*)'Comoving Radius:',VoidComovDistMPC(0.d0,GARadius,VP%t0,0.d0)
      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=0.d0,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=0:',For3 * aout**3/ (For3out * a**3  )-1

      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=GARadius,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=r_GA:',For3 * aout**3/ (For3out * a**3  )-1

      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=VP%L,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=L:',For3 * aout**3/ (For3out * a**3  )-1

      t=VP%t0
      garadius = VoidGARadius(t)
      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=2*VP%L,t=t,For3=For3out,a=aout)
        
      write(*,*)'at z=0:'
      write(*,*)'Comoving Radius:',VoidComovDistMPC(0.d0,GARadius,VP%t0,0.d0)
      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=0.d0,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=0:',For3 * aout**3/ (For3out * a**3  )-1

      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=GARadius,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=r_GA:',For3 * aout**3/ (For3out * a**3  )-1

      call lltb_functions(H0_inf=vp%h0,lambda=vp%lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
               dtbbdr=voiddtbbdr,Mtilde=VP%mtilde,r=VP%L,t=t,For3=For3,a=a)
      write(*,*)'delta M / M for r=L:',For3 * aout**3/ (For3out * a**3  )-1

  end subroutine testVoidTopHatGaussianField



  subroutine testVoidDumpDAofz()
    implicit none
    integer :: i
    logical :: exist, opened
    character(len=32), parameter :: fname = "VoidDAofZ.txt", form="(4ES14.5)"
    integer, parameter :: myfnum = 666
    real(dl) :: z, da, dl
    if(VOidError)return
    opened = .false.
    inquire(myfnum,opened=opened)
    if(.not.opened)then
       ! initialize file
       ! open file overwrite:
       open(myfnum,file=fname,status="replace")
       write(myfnum,'(A)')"# z  d_A(z)[Mpc]   d_L(z)[Mpc"
    end if

    write(myfnum,'("# ",10(A,":",ES14.5,", "))')"H0[km/s/Mpc]:",VP%H0,"omk",VP%omk,"omDE",VP%omv,"w",VP%w,"delta0",VP%delta0,"zB",VP%zB,"L(Mpc)", VoidComovDistMPC(0.d0,VP%L,VP%t0,VP%r0),"alpha",VP%L2_1
    do i = 1, size(VPNR%ztr(1,:))
       z = VPNR%ztr(1,i)
       if(z>1.1d3)exit
       da = VoidDA(z)
       dl = (1+z)**2 * da
       write(myfnum,form) z, da, dl
    end do
    write(myfnum,*)""

  end subroutine testVoidDumpDAofz




