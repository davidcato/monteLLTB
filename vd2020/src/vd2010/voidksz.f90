!*
! *  voidksz.f90
! *  VoidDistancesII
! *
! *  Created by Wessel Valkenburg on 12/04/2012.
! *  Copyright 2012 Wessel Valkenburg. All rights reserved.
! *
! */

! 


! kflrw returns the effective wavenumber k for which to
! return the matter power spectrum, as a function of 
! the actual wavenumber kltb and redshift z.
! Initime is set as a parameter, the time at which
! LTB ~= FLRW.
! function kflrw(kltb,z)
function kflrw(z)
  implicit none
  real(dl) :: kflrw
  ! real(dL), intent(in) :: kltb,z
  real(dL), intent(in) :: z
  real(dL) :: initime
  real(dl), save :: cache(5)
  real(dl) :: aFz, aFi, dAoverY, HFz, Yi, Yz, Ypi, Ypz, Ypdz, ai, az, dAF, dALTB
  real(dl) :: rLTB, tFLRW, tLTB, rmpc, kltb
  real(dl) :: Lmem
  real(dl), parameter :: third = 1.d0/3.d0
  
  initime = VP%t0*1.d-5 ! Test this?

  ! Get initial time quantities. This can be cached.
  ! FLRW:
  Lmem = VP%L
  VP%L = 0.d0
  
  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
      r=1.d-2*VP%L,t=initime,&
      a=aFi)

  tFLRW = nt_FLRW(z=z)

  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
      r=1.d-2*VP%L,t=tFLRW,&
      a=aFz,H=HFz,Rltb=dAF) ! because for FLRW r_obs is irrelevant: dA = Y.
  
  VP%L = Lmem
  ! END FLRW
  
  !LTB
  rLTB = nr(z)

  tLTB = nt(z=z)
  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
      r=rLTB,t=initime,&
      a=ai,Rltb=Yi,Rp=Ypi)

  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
      r=rLTB,t=tLTB,&
      a=az,Rltb=Yz,Rp=Ypz,Rpd=Ypdz)

  dALTB = DA_unnorm(z)

  ! rmpc = nr(z) * VPc*1.d-3
  rmpc = rLTB * VPc*1.d-3 
  kltb = (2*dble(3000.) + 1.d0)/(2*rmpc)
  
!  kFLRW = kLTB * aFZ/aFi*( (Yi/Yz*dALTB/dAF)**2 * Ypi * HFz / Ypdz )**(third)
  kFLRW = kLTB * aFZ/aFi*( (Yi/Yz)**2 * Ypi / Ypz )**(third)
  
!  write(*,*)kFLRW,kLTB,kFLRW/kLTB
  
end function kflrw

subroutine vdi_kSZSetRedshifts(numz, zarray)
  implicit none
  integer, intent(inout) :: numz
  real(dl), intent(inout) :: zarray(:)
  integer :: inz, i, imax, lowb, maxnz
  real(dl) :: inzarray(size(zarray))
  
  if(numz==0)then ! This should actually never happen..
    numz=1
    zarray(1) = 0.d0
  end if
  maxnz=size(zarray)

  inz = numz
  inzarray = zarray
  
  ! just a guess to start, 10 in ltb, 10 out of ltb:
  imax = 10
  numz = 0
  do i = 1, imax
    ! don't include zero: that should already be in inzarray
    ! **0.7 because we expect more change close to border, more samples there.
    zarray(i) = VP%zb * (dble(i)/dble(imax))!**(0.9d0) ! ^.9 to have more samples close to edge = collapsing shell
    numz = numz + 1
  end do

  lowb = numz
  imax = maxnz - lowb - inz
  do i = 1, imax
!    zarray(lowb+i) = VP%zb + (1.1d3 - VP%zB) * (dble(i)/dble(imax))**5 ! just a nice number to have small steps first.
    zarray(lowb+i) = VP%zb + (VP%kSZzmax*1.01 - VP%zB) * (dble(i)/dble(imax))**0.7 ! 0.7 to have more samples close to where LTB patch is on LSS of observer.
    numz = numz + 1
  end do
  
  lowb = numz + 1
  numz = numz + inz ! because inz >= 1
  zarray(lowb:numz) = inzarray(1:inz) ! add initially present redshifts
  
  call vdi_kSZQsort(zarray(numz:1:-1)) ! inverse sort redshifts
  
end subroutine vdi_kSZSetRedshifts

! QuickSort for double real array.
RECURSIVE SUBROUTINE vdi_kSZQsort(a)

  real(dl), INTENT(IN OUT) :: a(:)
  INTEGER :: split

  IF(size(a) > 1) THEN
     CALL vdi_kSZPartition(a, split)
     CALL vdi_kSZQsort(a(:split-1))
     CALL vdi_kSZQsort(a(split:))
  END IF

END SUBROUTINE vdi_kSZQsort

SUBROUTINE vdi_kSZPartition(a, marker)

  real(dl), INTENT(IN OUT) :: a(:)
  INTEGER, INTENT(OUT) :: marker
  INTEGER :: left, right
  real(dl) :: pivot
  real(dl) :: temp

  pivot = (a(1) + a(size(a))) / 2  ! Average of first and last elements to prevent quadratic 
  left = 0                         ! behavior with sorted or reverse sorted data
  right = size(a) + 1

  DO WHILE (left < right)
     right = right - 1
     DO WHILE (a(right) > pivot)
        right = right-1
     END DO
     left = left + 1
     DO WHILE (a(left) < pivot)
        left = left + 1
     END DO
     IF (left < right) THEN 
        temp = a(left)
        a(left) = a(right)
        a(right) = temp
     END IF
  END DO

  IF (left == right) THEN
     marker = left + 1
  ELSE
     marker = left
  END IF

END SUBROUTINE vdi_kSZPartition

! subroutine vdi_kSZAddkSZToClScalar(clscal,YHe,Omb,Omc,CData)
!   use PRELTB()camb
!   implicit none
!   real(dl), intent(inout) :: clscal(2:)
!   real(dl), intent(in) :: YHe, Omb, Omc
!   type (PRELTB()CAMBdata) :: CData
!   real(dl) :: redshifts(100), fb
!   integer :: numz
!   real(dl), allocatable :: integrand(:), dipoleOfZ(:), rofz(:), tOfZ(:), FOfz(:), betaOfZ(:)
!   real(dl), allocatable :: dtaudr(:), y2(:), yint(:)
!   integer, parameter :: numl = 7, lmin_vdi = 2
!   real(dl) :: alll(numl)
!   real(dl) :: kSZCl(numl),d2kSZCl(numl), thisCl
!   integer :: il, lmax, iz, thisl
!   real(dl) :: r, deltat ! absolute deltat and signed deltat
!   real(Dl) :: tLSS, dbl, simplint
! !  type(VPNR_type) :: myVPNR, myVPNR_FLRW
!   logical :: alloutside
!   logical, parameter :: kSZFeedback=.false.

!   integer :: thetime(3), ndig
!   character(len=32) :: timeformat
! !
!   if(voidError)return
!   if(kSZFeedback.or.VoidFeedback>0)  write(*,*)'Calculating kSZ Spectrum'
!   if(kSZFeedback.or.VoidFeedback>1)call system_clock(count=thetime(1),count_rate=thetime(3))

!   if(.not.associated(VPNR%ztr))return

!   redshifts = 0.d0

!   numz = CP%transfer%num_redshifts
!   redshifts(1:numz) = CP%transfer%redshifts(1:numz)
  
!   allocate(integrand(numz))
!   allocate(dipoleOfZ(numz))
!   allocate(rOfZ(numz))
!   allocate(tOfZ(numz))
!   allocate(FOfZ(numz))
!   allocate(betaOfZ(numz))
!   allocate(dtaudr(numz))
!   allocate(y2(numz))
!   allocate(yint(numz))

!   ! pre-fetch these arrays, because we are going to mess with VPNR for the dipoles.
!   do iz = 1, numz
!     rOfZ(iz) = nr(z=redshifts(iz)) ! re-use this in MakeIntegrand
!     tOfZ(iz) = nt(z=redshifts(iz))
!   end do
!     ! May be healthy for observer at r<0. CHeck was only for debugging.
!     !  if(any(rOfZ<0))then
!     !    write(*,*)vpnr%monotonic_in_z
!     !    call VoidErrorReport()
!     !    stop 'r<0'
!     !  end if

!   tLSS = nt(z=1100.d0) ! Last Scattering Surface

!   ! Get dipole of z for all redshifts(:)
!   if(kSZFeedback.or.VoidFeedback>1)write(*,*)'Getting dipole(z)'
  
!   ! Backup current results
! ! NO LONGER necessary; not using dothemagic in the dipoleoncone routine.
! !  call nullifyThisVPNR(myVPNR)
! !  call nullifyThisVPNR(myVPNR_FLRW)
! !  call associateVPNRToThisVPNR(VPNR,myVPNR)
! !  call associateVPNRToThisVPNR(VPNR_FLRW,myVPNR_FLRW)
! !  call nullifyThisVPNR(VPNR)
! !  call nullifyThisVPNR(VPNR_FLRW)
  
!   dipoleOfZ=0.d0
!   betaOfZ=0.d0
!   FOfZ=0.d0
  
!   fb = omb / (omb+omc) ! Baryon fraction
!   if(VoidTesting)SetLBoost = 1.d-4
!   ! Prepare the r-dependent functions
!   if(VoidTesting)then
!      open(12345,file="VoidDipoleOnCone_kSZ.txt",status="REPLACE")
!      write(12345,*)0.,0.,'a'
!   end if
!   do iz = numz, 1, -1
!     if(rOfZ(iz)==0.d0)cycle
!     if(VoidTesting.and.iz<numz)then ! double spacing
!        dbl=0.5*(redshifts(iz)+redshifts(iz+1))
! !       call VoidAsympDipoleOnCone(dbl, nr(z=dbl), nt(z=dbl), tLSS, 1100.d0, simplint,alloutside)! absolute deltat and signed deltat
!        call VoidAsympDipoleOnCone(nr(z=dbl), nt(z=dbl), tLSS, 1100.d0, simplint,alloutside)! absolute deltat and signed deltat ! z removed
!        write(12345,*)dbl,simplint,'a'
!     end if
! !    call VoidAsympDipoleOnCone(redshifts(iz), rOfZ(iz), tOfZ(iz), tLSS, 1100.d0, dipoleOfZ(iz),alloutside)! absolute deltat and signed deltat
!     call VoidAsympDipoleOnCone(rOfZ(iz), tOfZ(iz), tLSS, 1100.d0, dipoleOfZ(iz),alloutside)! absolute deltat and signed deltat ! z removed
!     ! write(*,*)redshifts(iz),dipoleOfZ(iz)
!     if(VoidTesting)write(12345,*)redshifts(iz),dipoleOfZ(iz)
!     if(alloutside .or. VoidError)then ! if alloutside Then both directions never cross the LTB patch anymore (r = large), and dipole is identically zero.
!       ! But now we pre-calculate VP%kSZzmax, so maxval(redshifts) should be very close to this redshifts(iz).
!       dipoleOfZ(1:iz)=0.d0
!       exit
!     end if
!     betaOfZ(iz) = dipoleOfZ(iz)*0.5d0
!     dtaudr(iz) = vdi_kSZdtaudr(redshifts(iz),tOfZ(iz),rOfZ(iz),YHe,fb)
!   end do
! !  betaOfZ(:) = dipoleOfZ(:)*0.5d0
!   FOfZ = betaOfZ * dtaudr
  
!   if(VoidTesting)then
!      write(*,*)'Dipole on cone written to VoidDipoleOnCone_kSZ.txt'
!      close(12345)
!   end if

!   call vdi_SetComptonY(numz,CData%Params%Reion%redshift,dtaudr,betaOfZ,rOfZ)

!   ! restore backup:
! !  call deallocateThisVPNR(VPNR)
! !  call deallocateThisVPNR(VPNR_FLRW)
! !  call associateVPNRToThisVPNR(myVPNR,VPNR)
! !  call associateVPNRToThisVPNR(myVPNR_FLRW,VPNR_FLRW)
! !  call nullifyThisVPNR(myVPNR)
! !  call nullifyThisVPNR(myVPNR_FLRW)
  
!   if(.not.VoidError)then

!     lmax = ubound(clscal,1)

!     ! For each l
!     ! Get list of integrand values, then spline-integrate.
!     if(kSZFeedback.or.VoidFeedback>1)  writE(*,*)'Calculating kSZCl for selected multipoles'!,VoidError,associated(VPNR_FLRW%ztr),associated(VPNR%ztr),nt_FLRW(z=0.1d0)
!     do il = 1, numl
!       thisl = nint( dble(lmax-lmin_vdi)*dble(il-1)/dble(numl-1) + lmin_vdi )
!       alll(il) = thisl
!       dbl=dble(thisl)
!       call vdi_kSZMakeIntegrandForl(thisl,numz,redshifts(1:numz),FOfZ,rOfZ,CData%MTrans,CData%Params%H0/1.d2,integrand)
!   !    kSZCl(il) = splintegrate..
!   !      subroutine spline_integrate(x,y,y2,yint,n)
!       call spline(rOfZ(numz:1:-1),integrand(numz:1:-1),numz,0.d0,0.d0,y2) ! growing r
!       call spline_integrate(rOfZ(numz:1:-1),integrand(numz:1:-1),y2,yint,numz) ! growing r

!       simplint = sum((integrand(2:numz)+integrand(1:numz-1))*0.5*(rOfZ(1:numz-1)-rOfZ(2:numz)))
!       if(abs(yint(numz)/simplint-1.d0)>1.d-1)then 
!         if(VoidFeedback>1)write(*,'(A,ES14.5,A,ES14.5)')'oscillatory behaviour in kSZ: taking simple integral ',simplint,'instead of',yint(numz)
!         yint = simplint
!       end if

!       kSZCl(il) = 16.d0*pi**2 / (2*dbl+1.d0)**3 * yint(numz)
      
      

!   !    write(*,*)thisl,kSZCl(il) * dbl * (dbl + 1.d0) * 2.273d6**2 / (2*pi),clscal(thisl)  * 2.273d6**2
!   !    write(1234,*)thisl,kSZCl(il) * dbl * (dbl + 1.d0)/ (2*pi) * 2.273d6**2 ,clscal(thisl)  * 2.273d6**2 
!       ! clscal from camb is normalized to clscal = C_l * l(l+1)/(2pi)
!       kSZCl(il) = kSZCl(il) * dbl * (dbl + 1.d0)/ (2*pi)
!     end do
  
!   else
!     kSZCl(:)=1.d10
!   end if


!   deallocate(integrand,dipoleOfZ,rOfZ,tOfZ,FOfZ,betaOfZ,dtaudr,y2,yint)


!   if(kSZFeedback.or.VoidFeedback>1)  write(*,*)'Finished calculating kSZ Spectrum'


!   if(.not.VoidError)then
!     ! fill the array with interpolated results
!     call spline(alll,kSZCl,numl,0.d0,0.d0,d2kSZCl) ! growing r  

!     do il = lmin_vdi, lmax
!       thisCl = voidsplint(alll,kSZCl,d2kSZCl,dble(il))
! !      if(modulo(il,100)==0) write(*,'(I7,20ES14.5)')il,clscal(il),thisCl,clscal(il) + thisCl
! ! write(1234,'(I7,20ES14.5)')il,clscal(il),thisCl,clscal(il) + thisCl
!       if(il==3000)then 
!         VP%kSZ3000 = thisCl
!         VP%kSZIsSet = .true.
!       end if
!       clscal(il) = clscal(il) + thisCl
!     end do
! ! write(1234,*)
!     if(any(isnan(clscal)))then
!       call VoidErrorReport()
!       stop 'kanker in je clscal'
!     end if
!   else
!     clscal = clscal + 1.d10
!   end if
  
!   if(kSZFeedback.or.VoidFeedback>1)then
!     call system_clock(count=thetime(2),count_rate=thetime(3))
!     ndig = int(log10(dble(thetime(2)-thetime(1))/dble(thetime(3))))
!     write(timeformat,*)ndig+6
!     if(ndig<0 .or. ndig>20)ndig=20
!     timeformat='(A,F'//trim(adjustl(timeformat))//'.3,A)' !'(A,F6.3,A)'
!     write(*,timeformat)' Done kSZ in ',dble(thetime(2)-thetime(1))/dble(thetime(3)),' seconds.'
!   end if
! !stop 'ksz ended'
  
!   if(voidtesting)then
!      write(*,*)'kSZ Cl (l=3000):',VP%kSZ3000
!      stop 'Done kSZ.'
!   end if
! end subroutine vdi_kSZAddkSZToClScalar

! subroutine vdi_kSZMakeIntegrandForl(thisl,numz,redshifts,FOfZ,rOfZ,MTrans,h,integrand)
!   use PRELTB()camb
!   implicit none
!   integer, intent(in) :: thisl, numz
!   real(dl), intent(in) :: redshifts(:), FOfZ(:), rOfZ(:),h
!   real(dl), intent(out) :: integrand(:)
!   Type(MatterTransferData), intent(in) :: MTrans
  
!   integer :: iitf, cambitf
!   integer, parameter ::  in = 1, npoints = 4
!   real, parameter ::     dlnkh = 1.d-2
!   real :: minkh, outpower(npoints), khmax
!   real(dl) :: kk, kkh, kltb,z, rmpc, rmpch, cambkh
  
!   call vdi_kSZCheckSizes(redshifts,numz,'redshifts','vdi_kSZMakeIntegrandForl')
!   call vdi_kSZCheckSizes(FOfZ,numz,'FOfZ','vdi_kSZMakeIntegrandForl')
!   call vdi_kSZCheckSizes(rOfZ,numz,'rOfZ','vdi_kSZMakeIntegrandForl')
!   call vdi_kSZCheckSizes(integrand,numz,'integrand','vdi_kSZMakeIntegrandForl')
    
    
!   integrand=0.d0

!   ! outpower: 1D array with power at kh(:) ranging from
!   ! minkh: with steps
!   ! dlnkh: to maximum log(kh)=log(minkh) + npoints * dlnkh
!   ! npoints: size of outpower
!   ! itf: iterates over redshifts(:)
!   ! in: iterates over spectra, always one in our simple world.
  
!   !write(*,*)'integrand'
!   do iitf = 1, numz
!     if(redshifts(iitf)<=1.d-10)cycle ! redshift zero.
!     rmpc = rOfZ(iitf) * VPc*1.d-3 ! c[kms]
!     rmpch = rmpc*h ! r in Mpc / h
!     if(iitf>numz)stop 'kankerdekanker'
!     kltb = (2*dble(thisl) + 1.d0)/(2*rmpch)

! !   write(*,*)kltb
!     kkh = kflrw(kltb,redshifts(iitf))
!     outpower=0.d0


!     khmax = MTrans%TransferData(1,MTrans%num_q_trans,iitf)   / exp(npoints*dlnkh)

!     cambkh = min(kkh,real(khmax,dl)) ! prefer constant over extrapolation or zero.
!     minkh = cambkh

!     if(iitf>numz)stop 'kankerdekanker'
!     if(minkh<0)then
!       call VoidErrorReport()
!       stop 'minkh<0'
!     end if

!     cambitf = iitf
!     ! Bug in my code.. cambitf=iitf no longer necessary.
!     call PRELTB()Transfer_GetMatterPower(MTrans,outpower, cambitf, in, minkh, dlnkh, npoints)
!     ! if (kkh>khmax)outpower=0.d0

!     integrand(iitf) = rOfZ(iitf) * FOfZ(iitf)**2 * outpower(1) * kkh**3 / (2 * pi**2)
!     !    write(*,'(20ES14.5)')redshifts(iitf),integrand(iitf),kkh,cambkh,rOfZ(iitf) , FOfZ(iitf)**2 ,outpower(1)!, iitf
!   end do
  
  
! end subroutine vdi_kSZMakeIntegrandForl

function dtaudrz(zi, YHe, omb, omc, H0)
  implicit none
  !real(dl), intent(inout) :: clscal(2:)
  real(dl) :: zi, h, YHe
  real omb, omc, H0
  !type (PRELTB()CAMBdata) :: CData
  real(dl) :: fb, Omegab, Omegac
  real(dl), allocatable :: rofZ, tOfZ
  real(dl), allocatable :: dtaudrz

  !YHe = 0.2454d0
  h = H0 / 100.d0

  Omegab = omb/h**2
  Omegac = omc/h**2

  rOfZ = nr(z=zi) ! re-use this in MakeIntegrand
  tOfZ = nt(z=zi)

  fb = Omegab / (Omegab+Omegac) ! Baryon fraction

  dtaudrz = vdi_kSZdtaudr(zi,tOfZ,rOfZ,YHe,fb)

end function dtaudrz

function betaz(zi,zreio)
  implicit none
  !real(dl), intent(inout) :: clscal(2:)
  !type (PRELTB()CAMBdata) :: CData
  real(dl) :: zi, zp, zreio
  !type (PRELTB()CAMBdata) :: CData
  logical :: alloutside
  real(dl) :: tLSS
  real(dl), allocatable :: rofZ, tOfZ
  real(dl), allocatable :: betaz, dipoleOfZ

  dipoleOfZ=0.d0
  betaz=0.d0
  tLSS = nt(z=1100.d0) ! Last Scattering Surface

  rOfZ = nr(z=zi) ! re-use this in MakeIntegrand
  tOfZ = nt(z=zi)
  
  zp = zreio+4.d0
  ! if(zi<VP%zB)then
    if(zreio.lt.zp)then ! computing for of all redshifts
      call VoidAsympDipoleOnCone(rofZ, tOfZ, tLSS, 1100.d0, dipoleOfZ,alloutside)! absolute deltat and signed deltat ! z removed
      if(alloutside)then
          betaz=0
      else
          betaz = dipoleOfZ*0.5d0
      end if
    ! else
    !   betaz = 0.d0
    end if
  ! call VoidAsympDipoleOnCone(rofZ, tOfZ, tLSS, 1100.d0, dipoleOfZ,alloutside)! absolute deltat and signed deltat ! z removed
  ! else
    ! betaz=0
  ! end if
  
  ! betaz = dipoleOfZ*0.5d0

end function betaz

subroutine vdi_kSZCheckSizes(arr,sz,arrname,callername)
  implicit none
  real(dl), intent(in) :: arr(:)
  integer, intent(in) :: sz
  character(len=*), intent(in) :: arrname, callername
  
  if(size(arr)/=sz)then
    write(0,*)'Array '//trim(adjustl(arrname))//' has wrong size in '//trim(adjustl(callername))//'.'
    stop
  end if
  
end subroutine vdi_kSZCheckSizes

function vdi_kSZdtaudr(z,t,r,YHe,fb)
  implicit none
  real(dl) :: vdi_kSZdtaudr
  real(dL), intent(in) :: z,t,r, YHe,fb
  
  real(dl) :: rhom, a,Yp, h2, rhom_sigt_over_mp,Xltb, Eofr, sq1p2e
  
  
  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
      r=r,t=t,&
      a=a,Rp=Yp,S=Xltb)
      
!  rhom = VP%Mtilde_paper**2 / (a**2 * Yp)
  ! Mtilde^2 = rho_crit = 11.2345 h^2 proton mass / m^3
  ! sigma_thomson = 6.65246 x 10^(-29) m^2
  h2 = VP%H0**2 * 1.d-4
  ! Mtilde to Mpc = c [km/s]
  ! -> dim[r] in code is Mpc s / km (i.e., code-r * c[kms] = r[Mpc])
  ! ideally we want dtaudr in code dimensions, 1/r = km / s / Mpc.

  ! rho_crit * sigma_thomson / proton mass = 7.47371 x 10^(-28) h^2 m^(-1)
  ! rho_crit * sigma_thomson / proton mass * c[m/s] = 2.24056 x 10^(-19) h^2 s^(-1)    (Mpc/Mpc)
  ! rho_crit * sigma_thomson / proton mass * c[m/s] * Mpc[km] = 6.91365 h^2 km s^(-1) Mpc^(-1)
  ! rhom = VP%Mtilde**2 / (a**2 * Yp) = rho_crit / (a**2 * Yp)
  ! rhom * sigma_thomson / proton mass * c[m/s] * Mpc[km] = 6.91365/(a**2 * Yp) h^2 km s^(-1) Mpc^(-1)

  ! Alternatively:
  ! c sigma_thomson / mprot / GNewton * km[Mpc] = 5.79142d-3 Mpc s / km

! Next lines give same normalization:
!  rhom_sigt_over_mp = 5.79142d-3 / (a**2 * Yp)  * VP%Mtilde**2
  rhom_sigt_over_mp = 6.91365 / (a**2 * Yp) * h2

  Eofr = r**2 * Voidkofr(r) * VP%Mtilde**2

  sq1p2e =  (1.d0 + 2.d0 * Eofr)**0.5d0
  vdi_kSZdtaudr = 0.5d0 * rhom_sigt_over_mp * fb * (2.d0 - YHe) / sq1p2e

end function vdi_kSZdtaudr

! Disable
subroutine vdi_SetComptonY(numz,zre,dtaudr,betaOfZ,rOfZ)
  implicit none
  integer, intent(in) :: numz
  real(dl), intent(in) :: dtaudr(:),betaOfZ(:),rOfZ(:), zre
  !
  real(dl) :: integrand(numz), rreionization, y2(numz), yint(numz)
  real(dl) :: simplint
  
  call vdi_kSZCheckSizes(dtaudr,numz,'dtaudr','vdi_SetComptonY')
  call vdi_kSZCheckSizes(betaOfZ,numz,'betaOfZ','vdi_SetComptonY')
  call vdi_kSZCheckSizes(rOfZ,numz,'rOfZ','vdi_SetComptonY')


  integrand = dtaudr * betaOfZ**2
  rreionization = NR(z=zre)
  where (rOfZ > rreionization) integrand = 0.d0 ! Integrate up to z_reionization

    ! call spline(rOfZ(numz:1:-1),integrand(numz:1:-1),numz,0.d0,0.d0,y2) ! growing r
    ! call spline_integrate(rOfZ(numz:1:-1),integrand(numz:1:-1),y2,yint,numz) ! growing r
  ! Sometimes no splines: sudden zero gives oscillations which make
  ! the second derivative crazy -> wrong result.
  simplint = sum((integrand(2:numz)+integrand(1:numz-1))*0.5*(rOfZ(1:numz-1)-rOfZ(2:numz)))
  if(abs(yint(numz)/simplint-1.d0)>1.d-1)then 
    if(VoidFeedback>1)write(*,'(A,ES14.5,A,ES14.5)')'oscillatory behaviour in compton-y: taking simple integral ',simplint,'instead of',yint(numz)
    yint = simplint
  end if
  
  VP%ComptonY = 0.7d0 * yint(numz)
  VP%ComptonYIsSet = .true.
  if(VoidFeedback>1)write(*,*)'Compton-y is set:',VP%ComptonY

  ! open(12340,file="ComptonY.txt",status="REPLACE")
  ! ! write(12340,*) 'Compton Y:'
  ! write(12340,'(100E20.11)') VP%ComptonY
  ! close(12340)

end subroutine vdi_SetComptonY

function vdi_ComptonY()
  implicit none
  real(dl) :: vdi_ComptonY
  
  if(VP%ComptonYIsSet .or. VoidError)then
    if(VoidError)then
      vdi_ComptonY = 1.d30
    else
      vdi_ComptonY = VP%ComptonY
    end if
  else
    call VoidErrorReport()
    write(0,*)'Asking for Compton-y while it has not been set.'
    stop
  end if

end function vdi_ComptonY


function vdi_kSZ3000()
  implicit none
  real(dl) :: vdi_kSZ3000
  real(dl) :: nm
  if(VP%kSZIsSet .or. VoidError)then
    if(VoidError)then
      vdi_kSZ3000 = 1.d30
    else
    ! see CMB_Cls_Simple.f90 : Theory%cl(l,1:num_clsS) =  nm*Cl_scalar(l,1, scalClOrder(1:num_clsS))
    ! nm = cons/(l*(l+1))
    ! cons = (COBE_CMBTemp*1e6)**2*2*pi
    ! COBE_CMBTemp = 2.726
      nm = 2.726d6**2*2*pi/(3000.*3001.)
      vdi_kSZ3000 = VP%kSZ3000 * nm
    end if
  else
    call VoidErrorReport()
    write(0,*)'Asking for Compton-y while it has not been set.'
    stop
  end if

end function vdi_kSZ3000
