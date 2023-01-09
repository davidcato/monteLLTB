

subroutine vdi_GetIgnObsDEPars()
  implicit none
!  integer, parameter :: npar=2
!  real(dl) :: parameterranges(2*npar) ! min/max
!  real(dl) :: parameterresults(npar)
  integer :: npar=2
  real(dl), allocatable :: parameterranges(:) ! min/max
  real(dl), allocatable :: parametertrial(:) ! first trial
  real(dl), allocatable :: parameterresults(:)
  real(dl), allocatable :: parameterproposals(:) ! proposal densities.
!  real(dl) :: xranges(2) ! xminmax
  integer :: myj
    ! NOTE: take info from third column of table 14 in
    ! http://arxiv.org/abs/1111.1969 , which has combined sigma and number distribution
    ! in one new sigma.
    ! Implement equal spacing in redshift space in dvode integrator, then use redshift dependent
    ! sigmas.
    ! 2-form in order to be able to use asplint.
    real(dl), parameter :: DESredshifts(14) = (/0.0,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.2/)
    real(dl), parameter :: DESrms_over_sqN(14) = (/0.0350,0.0350,0.0140,0.0082,0.0073,0.0068,0.0075,0.0086,0.0120,0.0150,0.0180,0.0200,0.0240,0.0240/)
    ! splined once and for all:
    real(dl), parameter :: DESddrms_over_sqN(14) = (/   2.98812E+00,  -5.97623E+00,   3.83463E+00,  -2.42307E-01,   7.45959E-02,   1.83923E-01,  -9.02901E-02,   4.17238E-01,  -1.98662E-01,   1.37412E-01,  -3.50983E-01,   6.66520E-01,  -1.11510E+00,   5.57549E-01/)
    real(dl), parameter :: DESbinwidth = 1.d-1
    real(dl), allocatable :: sig(:,:)
    integer :: myi, ninside, thisbin, bincount(12)

    if(VoidIncludeFLRWDEFit<1)then
      write(0,*)'you called vdi_GetIgnObsDEPars without specifying VoidIncludeFLRWDEFit:',VoidIncludeFLRWDEFit
      write(0,*)'skipping..'
      return
    end if

  ! change of plan, each value of VoidIncludeFLRWDEFit
  ! represents a pre-defined combination of parameters.
  ! 1: omega_m + w (omega_k=0,wa=0)
  ! 2: omega_m + w + wa (omega_k=0)
  ! 3: omega_m + omega_k (w=-1,wa=0)
  ! 4: omega_m + omega_k + w (wa=0)
  
!  if(VoidIncludeFLRWDEFit>2)then
!    npar=VoidIncludeFLRWDEFit
!  else
!    npar=2
!  end if
  
  if(any(VoidIncludeFLRWDEFit==(/1,3/)))then
    npar=2
  else if(any(VoidIncludeFLRWDEFit==(/2,4/)))then
    npar=3
  end if

  allocate(parameterranges(2*npar),parameterresults(npar),parametertrial(npar),parameterproposals(npar))
 
  do myj=1,size(VPNR%ztr(1,:))
  ! NOTE: make this redshift upper bound an input parameter in params.ini
  ! WV Beijing 13/09/2012
!    if(VPNR%ztr(1,myj)>2.d0)exit
    if(VPNR%ztr(1,myj)>VoidIncludeFLRWDEFitzmax)exit
  end do
  myj=min(myj,size(VPNR%ztr(1,:))) ! limit from above.

!  xranges(1)=0.d0
!  xranges(2)=2.d0

  ! all cases:
!  parameterranges(1) = 0.d0 !omm min!
!  parameterranges(2) = 1.d0 !omm max!
  parametertrial(1) = 1.d0 - VP%omk - VP%omv
  parameterranges(1) = parametertrial(1)*0.9
  parameterranges(2) = parametertrial(1)*1.1
  parameterproposals(1) = 0.2
  if(any(VoidIncludeFLRWDEFit==(/1,2,4/)))then

    parameterranges(3) = VP%w*1.1! -2.d0 !w min!
    parameterranges(4) = VP%w*0.9!-1.d0/3.d0 !w max!
    parametertrial(2) = VP%w
    parameterproposals(2) = 0.2
  end if
  
  if(any(VoidIncludeFLRWDEFit==(/2,4/)))then
    parameterranges(5) = -1.d0 !wa or omk min!
    parameterranges(6) = 1.d0 !wa or omk max!
    if(VoidIncludeFLRWDEFit==2)then
      parametertrial(3) = 0.d0
      parameterproposals(3) = 0.2
    end if
    if(VoidIncludeFLRWDEFit==4)then
      parametertrial(3) = VP%omk
      parameterproposals(3) = 0.005
    end if
  end if

  if(VoidIncludeFLRWDEFit==3)then
    parameterranges(3) = -1.d-1!0 !omk min!
    parameterranges(4) = 1.d-1!0 !omk max!
    parametertrial(2) = 0.d0
    parameterproposals(2) = 0.005
  end if

if(any(parametertrial<parameterranges(1:2*npar-1:2)).or. any(parametertrial>parameterranges(2:2*npar:2)))then
  write(0,*)'ptrial:',parametertrial
  write(0,*)' you have fiducial parameters that are outside of the hardcoded parameter bounds in voidminchi2.f90'
  call VoidErrorReport()
  write(*,*)sig(1,1)
stop
end if
  ! prepare sigmas:
  allocate(sig(1,myj))
ninside=0
! normalized: no more spline.
  if(VoidIncludeFLRWDEFitWeighting==1)then
    bincount = 0
    ! first count bins
    do myi = 2, myj
      if(VPNR%ztr(1,myi)<DESredshifts(14))then
        thisbin = min(size(bincount),max(1,int(VPNR%ztr(1,myi)/DESbinwidth)+1))
      bincount(thisbin)=bincount(thisbin)+ 1
      end if
    end do
    do myi = 2, myj
      if(VPNR%ztr(1,myi)<DESredshifts(14))then
!        call asplint(DESredshifts,DESrms_over_sqN,VPNR%ztr(1,myi),sig(:,myi),DESddrms_over_sqN)
        thisbin = min(size(bincount),max(1,int(VPNR%ztr(1,myi)/DESbinwidth)+1))
        sig(:,myi) = DESrms_over_sqN(thisbin+1) * sqrt(real(bincount(thisbin)))
!        write(*,*)DESrms_over_sqN(1,thisbin+1),sig(:,myi),thisbin,bincount(thisbin)
      else
        if(ninside==0)ninside=myi-1
        sig(:,myi) = 1.d30
      end if
    end do

    if(VPNR%ztr(1,myj)>1.d3)then
      sig(:,myj)=3.2d-3 * VoidDA(VPNR%ztr(1,myj)) ! 0.32% error from Table 9 in http://arxiv.org/abs/astro-ph/0603451
    end if
  else if(VoidIncludeFLRWDEFitWeighting==0)then
    sig=1.d0
  end if
! DES sigmas are only normalized to 12 datapoints. So, renormalize:
!  if(ninside/=0)then
!    sig(1,1:ninside)=sig(1,1:ninside) * sqrt(dble(ninside)/12.d0) ! number must grow for ninside > 12.
!  end if

  ! need to use exact integration steps ztr, because otherwise
  ! the interpolation between the integration steps
  ! leads to a bias in w, always O(2%) too negative.

!  call vdi_GetParsMinimizeChi2(VoidDA,flrwda_z,VPNR%ztr(1,1:myj),sig(1,:),parameterranges, parametertrial,parameterresults,parameterproposals)

  call vdi_GetParsMinimizeChi2(VoidDAModulus,flrwda_z_Modulus,VPNR%ztr(1,2:myj),sig(1,2:myj),parameterranges, parametertrial,parameterresults,parameterproposals)
  if(VoidFeedback>0)then
    write(*,*)'omega_m:',parameterresults(1)
    if(any(VoidIncludeFLRWDEFit==(/1,2,4/)))write(*,*)'w:',parameterresults(2)
    if(VoidIncludeFLRWDEFit==2)write(*,*)'wa:',parameterresults(3)
    if(VoidIncludeFLRWDEFit==4)write(*,*)'omk:',parameterresults(3)
    if(VoidIncludeFLRWDEFit==3)write(*,*)'omk:',parameterresults(2)
  end if

!  stop 'w done'
! pausing:
!read(*,'("Type 1 and hit return to do next line",I2)')myj

!  VP%defit_w=parameterresults(2)
!  if(npar>2)VP%defit_wa=parameterresults(3)
  VP%defit_omm = parameterresults(1)
  VP%defit_w = -1
  VP%defit_wa = 0
  VP%defit_omk = 0
  if(any(VoidIncludeFLRWDEFit==(/1,2,4/)))VP%defit_w=parameterresults(2)
  if(VoidIncludeFLRWDEFit==2)VP%defit_wa=parameterresults(3)
  if(VoidIncludeFLRWDEFit==4)VP%defit_omk=parameterresults(3)
  if(VoidIncludeFLRWDEFit==3)VP%defit_omk=parameterresults(2)
  VP%defit_set=.true.


  deallocate(parameterranges,parameterresults,parametertrial,sig,parameterproposals)

  contains

    function VoidDAModulus(z)
      implicit none
      real(dl), intent(in) :: z
      real(dl) :: VoidDAModulus

      VoidDAModulus = Log( (1+z)**2 * VoidDA(z))
      if(z<1.d-14)VoidDAModulus=0.d0
      if(isnan(VoidDAModulus))VoidDAModulus=0.d0
    end function VoidDAModulus

    function flrwda_z(pars,zarr,sz)
      implicit none
      integer, intent(in) :: sz
      real(dl), intent(in) :: pars(:), zarr(:)
      real(dl) :: flrwda_z(sz), rofz(sz)
      real(dl) :: start, prevz, prevf
      integer :: i, n, imin
!    real(dl), external :: rombint

      n = size(zarr)
      imin=1
      start=zarr(1)
      prevf=0.d0
      if(zarr(1)==0.d0)then
        imin=2
        flrwda_z(1)=0.d0
      else
        imin=1
!        stop 'need to write first integration in vdi_GetIgnObsDEPars'
      end if
      prevz=start
      call storepars(put=pars)
      do i = imin, n
!        rofz(i)=prevf+(zarr(i)-prevz)*MyDrDz(pars,prevz)
if(isnan(f(zarr(i))))then
  rofz = f(zarr(i))
  return
end if
        rofz(i)=prevf+rombint(f,prevz,zarr(i),1.d-12)
        prevf=rofz(i)
        prevz=zarr(i)
        flrwda_z(i)=rofz(i)/(1.d0+zarr(i))
      end do
    
    end function flrwda_z

    function flrwda_z_Modulus(pars,zarr,sz)
      implicit none
      integer, intent(in) :: sz
      real(dl), intent(in) :: pars(:), zarr(:)
      real(dl) :: flrwda_z_Modulus(sz), rofz(sz)
      real(dl) :: start, prevz, prevf
      integer :: i, n, imin
!    real(dl), external :: rombint

      n = size(zarr)
      imin=1
      start=0.d0!zarr(1)
      prevf=0.d0
      flrwda_z_Modulus=0
      if(zarr(1)==0.d0)then
        imin=2
        flrwda_z_Modulus(1)=0.d0
      else
        imin=1
!        stop 'need to write first integration in vdi_GetIgnObsDEPars'
      end if
      prevz=start
      call storepars(put=pars)
      do i = imin, n
!        rofz(i)=prevf+(zarr(i)-prevz)*MyDrDz(pars,prevz)
if(isnan(f(zarr(i))))then
  rofz = f(zarr(i))
  return
end if
        rofz(i)=prevf+rombint(f,prevz,zarr(i),1.d-12)
        prevf=rofz(i)
        prevz=zarr(i)
        flrwda_z_Modulus(i)=Log(rofz(i) * (1.d0+zarr(i))*299792.458/VP%H0)
        if(zarr(i)<1.d-14)flrwda_z_Modulus(i) = 0.d0
        if(isnan(flrwda_z_Modulus(i)))flrwda_z_Modulus(i)=0.d0
      end do

    end function flrwda_z_Modulus

    subroutine storepars(put,get)
      real(dl), intent(in), optional :: put(:)
      real(dl), intent(out), optional :: get(:)
      real(dl), save :: stpars(100)
      integer, save :: parsize=0
      
      if(present(put))then
        stpars=0.d0
        parsize=size(put)
        stpars(1:parsize)=put
      end if
      if(present(get))then
        if(size(get)/=parsize)then
          write(0,*)'Size of get array is not equal to store array, or no array stored in storepars.'
          stop
        end if
        get=stpars(1:parsize)
      end if

    end subroutine storepars
    
    function MyDrDz(mpars,z)
      implicit none
      real(dl), intent(in) :: mpars(:),z
      real(dl) :: MyDrDz
      real(dl) :: omm, omk, oml, w, wa, a

!  if(any(VoidIncludeFLRWDEFit==(/1,2,4/)))write(*,*)'w:',parameterresults(2)
!  if(VoidIncludeFLRWDEFit==2)write(*,*)'wa:',parameterresults(3)
!  if(VoidIncludeFLRWDEFit==4)write(*,*)'omk:',parameterresults(3)
!  if(VoidIncludeFLRWDEFit==3)write(*,*)'omk:',parameterresults(2)

! defaults first:
  w=-1
  wa=0
  omk=0
  omm=mpars(1)
  if(any(VoidIncludeFLRWDEFit==(/1,2,4/)))w=mpars(2)
  if(VoidIncludeFLRWDEFit==2)wa=mpars(3)
  if(VoidIncludeFLRWDEFit==4)omk=mpars(3)
  if(VoidIncludeFLRWDEFit==3)omk=mpars(2)

!      omm=mpars(1)
      oml=1.d0-omm-omk
!      w=mpars(2)
!      if(size(mpars)>2)then
!        wa=mpars(3)
!      else
!        wa=0.d0
!      end if
      a=1.d0/(1.d0+z)
      
      MyDrDz = 1.d0/sqrt(omm*(1+z)**3 + omk*(1+z)**2 + oml*(1+z)**(3*(1.d0+w+(1-a)*wa)))
      
    end function MyDrDz

    function f(z)
      implicit none
      real(dl), intent(in) :: z
      real(dl) :: f
      real(dl) :: p(npar)
      
      call storepars(get=p)
      f=MyDrDz(p,z)
    
    end function f

end subroutine vdi_GetIgnObsDEPars


subroutine vdi_GetParsMinimizeChi2(realfunc,fitfunc,xranges,sigma,parameterranges,parametertrial,parameterresults,parameterproposals)
  implicit none
  interface
    function myfunc1(rr)
      integer, parameter :: dl=kind(1.d0)
      real(dl) :: myfunc1
      real(dl), intent(in) :: rr
    end function myfunc1
  end interface
  interface
    function myfunc2(p,rr,n)
      integer, parameter :: dl=kind(1.d0)
      integer, intent(in) :: n
      real(dl) :: myfunc2(n)
      real(dl), intent(in) :: rr(:)
      real(dl), intent(in) :: p(:)
    end function myfunc2
  end interface
  procedure(myfunc1)  :: realfunc
  procedure(myfunc2)  :: fitfunc
  real(dl), intent(in) :: xranges(:) ! array of points to sample
  real(dl), intent(in) :: sigma(:) ! array of sigmas for each x
  real(dl), intent(in) :: parameterranges(:)
  real(dl), intent(in) :: parametertrial(:),parameterproposals(:)
  real(dl), intent(out) :: parameterresults(:)
  integer, parameter :: nsidepar=16
  integer :: nside
  integer :: steps
  integer :: n, i, j, ipos, iii
  integer, allocatable :: pos(:), cpos(:), corners(:), newcorners(:)
  real(dl), allocatable :: chimatrix(:),parmatrix(:,:) ! fortran layout, arbitrarily dimensional, second dim is par-vector
  real(dl), allocatable :: realfuncarr(:), xarr(:), fitfuncarr(:), xi(:,:), p(:)
  real(dl) :: newranges(size(parameterranges))
!  real(dl) :: newgrid(size(parameterresults),nside), gridsteps(nside)
  real(dl), allocatable :: newgrid(:,:), gridsteps(:)
  real(dl) :: startvol, currvol, fret, temp
  type minsettype
    real(dl) :: chi2
    integer :: ipos
  end type minsettype
  type(minsettype) :: minset
  logical, parameter :: vchdebug=.false.

  n=size(parameterresults)
  nside = nint(nsidepar * (size(parameterresults)+2)/4.) ! because more params are less constrained we need higher reso. :(


  if(size(parameterranges)/=2*n)then
    write(0,*)'size of parameterresults and parameterranges do not match (mod(2)) in vdi_GetParsMinimizeChi.'
    stop
  end if

  steps=size(xranges)
  if(size(sigma)/=steps)then
    write(0,*)'array of sigmas must be same size as array of xranges in voidminch2.f90, vdi_GetParsMinimizeChi.'
    stop
  end if

  allocate(realfuncarr(steps), xarr(steps), fitfuncarr(steps))
  xarr=xranges
  ! generate realfuncarr and zarr:
  do i = 1, steps
!    xarr(i) = (i-1.d0)/(steps-1.d0)*(xranges(2)-xranges(1))+xranges(1)
    realfuncarr(i) = realfunc(xarr(i))
  end do
!write(*,*)'hack realfuncarr.'  
!realfuncarr = fitfunc((/0.3d0,-1.d0/),xarr,steps)

if(.false.)then

  allocate(newgrid(size(parameterresults),nside), gridsteps(nside))
  allocate(chimatrix(nside**n),parmatrix(nside**n,n))
  allocate(pos(n),cpos(n))
  allocate(corners(n**2), newcorners(n**2))
  do i = 1, n**2  ! two points per dimension
      ! 1 is one end, 2 is the other
      pos=serialtomatrix(i,n,2)
      ! get the serial positions of corners of nside**n matrix
      ipos=matrixtoserial(1+(pos(:)-1)*(nside-1),n,nside)
      corners(i)=ipos
  end do
  ! start grid walking
  ! first time get all corners manually
  do i = 1, n**2  ! two points per dimension
    ! 1 is one end, 2 is the other
    pos=serialtomatrix(i,n,2)
    ! get the serial positions of corners of nside**n matrix
    ipos=corners(i)!matrixtoserial(1+(pos(:)-1)*(nside-1),n,nside)
    !write(*,*)pos,matrixtoserial(pos,n,2)
    do j = 1, n
      parmatrix(ipos,j)=parameterranges(2*(j-1)+1)+(pos(j)-1)*(parameterranges(2*j)-parameterranges(2*(j-1)+1))
      !write(*,*)i,j,parmatrix(i,j)
    end do
    chimatrix(ipos)=mychi2(parmatrix(ipos,:))
!    write(*,*)parmatrix(ipos,:),chimatrix(ipos)
  end do
  
  do i = 1, nside
    gridsteps(i)=(i-1.d0)/dble(nside-1.d0)
  end do
  
  ! Next do the loop:
  ! put bounds of region to the corners,
  ! then calculate grid inside this regio
  ! then find smalles n-cube surrounding
  ! bestfit point.
  ! put new bounds in corners etc..
  ! compare starting par-space volume to current.
  startvol = 1
  do j = 1, n
    startvol=startvol*(parameterranges(2*j)-parameterranges(2*(j-1)+1))
  end do
  do iii=1,1000
    ! put corners in newranges
    do i = 1, n
      newranges(2*(i-1)+1) = minval(parmatrix(corners,i))
      newranges(2*i) = maxval(parmatrix(corners,i))
      newgrid(i,:)=gridsteps*(newranges(2*i)-newranges(2*(i-1)+1))+newranges(2*(i-1)+1)
    end do
    currvol = 1
    do j = 1, n
      currvol=currvol*(newranges(2*j)-newranges(2*(j-1)+1))
    end do
    if(vchdebug)write(*,*)'0:',currvol/startvol
    if(abs(currvol/startvol)<1.d-5)exit
    if(vchdebug)write(*,'(2ES14.5)')newranges
    ! walk through matrix
    do i = 1, nside**n
      if(all(i/=corners))then ! we have corners already
        pos=serialtomatrix(i,n,nside)
        do j = 1, n
          parmatrix(i,j)=newgrid(j,pos(j))
        end do     
        chimatrix(i)=mychi2(parmatrix(i,:))
      end if
      ! get smallest value on the fly, start with first point
      if(chimatrix(i)<minset%chi2 .or. i==1 .or. isnan(minset%chi2))then
        minset%chi2=chimatrix(i)
        minset%ipos=i
      end if
      if(vchdebug)write(*,*)'a:',chimatrix(i),minset%chi2,parmatrix(i,:)
    end do
    ! get ipos of neighbouring points around minval(chimatrix)
    pos=serialtomatrix(minset%ipos,n,nside)
    if(vchdebug)write(*,*)'b:',pos
    do i = 1, n**2
      cpos=serialtomatrix(i,n,2) ! 1,2
      cpos=2*cpos-3 ! -1,1
      cpos=cpos*2 ! -2,2, taking it slower
      cpos=pos+cpos ! now its one of the corners, module border lines
      where(cpos<1)cpos=1
      where(cpos>nside)cpos=nside 
      if(vchdebug)write(*,*)'b1:',cpos
      ! now fixed to borders
      newcorners(i) = matrixtoserial(cpos,n,nside)
    end do
    if(vchdebug)write(*,*)'c:',newcorners
    ! now we got the corners. Move new corners to corners of matrix.
    chimatrix(corners)=chimatrix(newcorners)
    parmatrix(corners,:)=parmatrix(newcorners,:) 
    ! done
    if(vchdebug)write(*,'(4ES14.5)')newranges
  
  end do
  
!  allocate(chimatrix(nside**n),parmatrix(nside**n,n))
!  allocate(pos(n),cpos(n))
!  allocate(corners(n**2), newcorners(n**2))
  parameterresults = 0.5d0*(newranges(1:2*n-1:2)+newranges(2:2*n:2))
  
  deallocate(chimatrix,parmatrix,pos,cpos,corners,newcorners)
!  deallocate(realfuncarr,xarr,fitfuncarr)
  deallocate(newgrid,gridsteps)
  if(vchdebug)then
    write(*,*)'old method:'
    write(*,'(10ES14.5)')parameterresults
    write(*,*)mychi2(parameterresults)
  end if
end if


  allocate(xi(n,n),p(n))

  p=parametertrial !(parameterranges(2:2*n:2)+parameterranges(1:2*n-1:2))*0.5d0
  xi=0.d0
  do i = 1, n
    xi(i,i) = (parameterranges(2*i)-parameterranges(2*i-1))*0.3
  end do
  call void_nr_powell(p,xi,1.d-4,i,fret,mychi2)
  parameterresults = p
  if(vchdebug)then
    write(*,*)'direct method:'
    write(*,'(10ES14.5)')parameterresults
    write(*,*)fret
    write(*,*)'chi2:',mychi2(parameterresults)
    write(*,*)'trial chi2:',mychi2(parametertrial)
  end if

!  p=parametertrial !(parameterranges(2:2*n:2)+parameterranges(1:2*n-1:2))*0.5d0
!  p=(parameterranges(2:2*n:2)+parameterranges(1:2*n-1:2))*0.5d0
!  xi=0.d0
!  do i = 1, n
!    xi(i,i) = (parameterranges(2*i)-parameterranges(2*i-1))*0.3
!  end do
!  ! find best starting point:
!  call void_nr_powell(p,xi,1.d-1,i,fret,trialchi2)
!  ! repeat last starting point:
!  call void_nr_powell(p,xi,1.d-4,i,fret,mychi2)
!  parameterresults = p
!  if(vchdebug)then
!    write(*,*)'recursive method:'
!    write(*,'(10ES14.5)')parameterresults
!    write(*,*)fret
!    write(*,*)'chi2:',mychi2(parameterresults)
!  end if

!  p = parameterresults
!  temp = mychi2(p)/steps ! force chi2/datapoint ~ 1
!  call vdi_minimc(p,parameterranges,parameterproposals,temp,mychi2)
!  parameterresults = p
!  if(vchdebug)then
!    write(*,*)'mc method:'
!    write(*,'(10ES14.5)')parameterresults
!    write(*,*)'chi2:',mychi2(parameterresults)
!    write(*,*)'trial chi2:',mychi2(parametertrial)
!  end if

  deallocate(realfuncarr,xarr,fitfuncarr,xi,p)

!  stop
contains

  function trialchi2(tpars) result(myfret)
    implicit none
    real(dl), intent(in) :: tpars(:)
    real(dl) :: myxi(size(tpars),size(tpars)), myp(size(tpars)), myfret
    integer :: i

    myxi=0.d0
    do i = 1, n
      myxi(i,i) = (parameterranges(2*i)-parameterranges(2*i-1))*0.3
    end do
    myp = tpars
    if(mychi2(myp)<1.d29)then
      call void_nr_powell(myp,myxi,1.d-4,i,myfret,mychi2)
    else
      myfret = 1.d30
    end if

!write(*,*)'trialchi2:',myfret
  end function trialchi2

  function mychi2(tpars)
USE iso_c_binding
    implicit none
interface
      function fork() bind(c, name='fork')
        integer :: fork
      end function fork
      subroutine exec(comm) bind(c, name='execlp')
        import
        character(kind=c_char), intent(in) :: comm(*)
      end subroutine exec
end interface
    real(dl), intent(in) :: tpars(:)
    real(Dl) :: mychi2, scale
    integer :: myj, pid
    real(dl) :: S0, S1
    logical, parameter :: liveplot=.false.
    logical, save :: firstplot=.true.
    
    if(any(tpars*(1.d0-(sign(1.d-5,tpars)))>parameterranges(2:2*n:2)).or.any(tpars*(1.d0+(sign(1.d-5,tpars)))<parameterranges(1:2*n-1:2)))then
!      write(*,*)'========='
!      write(*,*)tpars*(1.d0-(sign(1.d-5,tpars)))>parameterranges(2:2*n:2)
!      write(*,*)tpars*(1.d0+(sign(1.d-5,tpars)))<parameterranges(2-1:2*n-1:2)
!      write(*,*)tpars
!      write(*,*)parameterranges(2:2*n:2)
!      write(*,*)parameterranges(2-1:2*n-1:2)
!      write(*,*)'========='
!      stop
      mychi2=1.d30
      return
    end if

    fitfuncarr = fitfunc(tpars,xarr,steps)

    if(.not.any(isnan(fitfuncarr)))then

! Fitting DA, marginalizing over x in DA -> x DA
!      scale=sum(fitfuncarr*realfuncarr)/sum(fitfuncarr**2)
!      mychi2=sum((realfuncarr-scale*fitfuncarr)**2/sigma**2)

! Fitting mu = log((1+z)^2 DA), marginalizing over y in mu -> mu + y
      S1 = sum((realfuncarr-fitfuncarr)/sigma**2)
      S0 = sum(sigma**(-2))
      scale = S1**2/S0
      mychi2= sum ((( realfuncarr - fitfuncarr)/sigma)**2) - scale

      if(liveplot)then
        if(firstplot)then
        write(*,*)'calling system'
  !        write(*,'(A)')'In a separate terminal in the same path, paste the following to see the fitting:'
  !        write(*,'(A)')'if [[ 1 -ne 2 ]]; then echo "set xrange [0:2]"; echo "set yrange [0:2000]"; echo "set y2range [0:50]" ; while [[ 1 -ne 2 ]]; do sleep 0.05; echo "plot ''fort.1234'' u 1:(\$2-\$3)**2 axis x1y2 w boxes fs solid 0.25 ls 3 t ''chi^2'', '''' u 1:(\$2-\$3)**2 sm cs axis x1y2 ls 4 t '''', '''' u 1:2 ls 1 t ''LTB theory'', '''' u 1:3 w l ls 2 t ''FLRW fit''"; done; fi | gnuplot &> /dev/null'
          ! fortran 2008
  !        call EXECUTE_COMMAND_LINE('if [[ 1 -ne 2 ]]; then echo "set xrange [0:2]"; echo "set yrange [0:2000]"; echo "set y2range [0:50]" ; while [[ 1 -ne 2 ]]; do sleep 0.05; echo "plot ''fort.1234'' u 1:(\$2-\$3)**2 axis x1y2 w boxes fs solid 0.25 ls 3 t ''chi^2'', '''' u 1:(\$2-\$3)**2 sm cs axis x1y2 ls 4 t '''', '''' u 1:2 ls 1 t ''LTB theory'', '''' u 1:3 w l ls 2 t ''FLRW fit''"; done; fi | gnuplot &> /dev/null ', &
  !              WAIT=.false.)
          firstplot=.false.
          ! no fortran 2008 yet, so do c-binding, works perfectly on os x
          ! because child process dies with parent.
          pid=fork()
          write(*,*)pid
          if(pid==0)then ! child process
!            call system('if  gnuplot -e q > /dev/null ; then echo "set xrange [0:2]"; echo "set yrange [0:2000]"; echo "set y2range [0:50]" ; while [[ 1 -ne 2 ]]; do sleep 0.05; echo "plot ''fort.1234'' u 1:(\$2-\$3)**2 axis x1y2 w boxes fs solid 0.25 ls 3 t ''chi^2'', '''' u 1:(\$2-\$3)**2 sm cs axis x1y2 ls 4 t '''', '''' u 1:2 ls 1 t ''LTB theory'', '''' u 1:3 w l ls 2 t ''FLRW fit''"; done; fi | gnuplot &> /dev/null ')
            call system('if  gnuplot -e q > /dev/null ; then echo "set xrange [0:2]";  echo "set y2range [0:50]" ; while [[ 1 -ne 2 ]]; do sleep 0.05; echo "plot ''fort.1234'' u 1:(\$2-\$3)**2 axis x1y2 w boxes fs solid 0.25 ls 3 t ''chi^2'', '''' u 1:(\$2-\$3)**2 sm cs axis x1y2 ls 4 t '''', '''' u 1:2 ls 1 t ''LTB theory'', '''' u 1:3 w l ls 2 t ''FLRW fit''"; done; fi | gnuplot &> /dev/null ')
            stop
          end if
        write(*,*)'calling done'
        end if
        open(1234)
        do myj = 1, steps
!            write(1234,'(20ES14.5)')xarr(myj),realfuncarr(myj),scale*fitfuncarr(myj)
            write(1234,'(20ES14.5)')xarr(myj),realfuncarr(myj),fitfuncarr(myj)
        end do
        close(1234)
      end if
  !    read(*,'("Type 1 and hit return to do next line",I2)')myj

    else
      mychi2 = 1.d30
    end if
!write(*,'(A,20ES14.5)')'mychi2:',mychi2,tpars
  end function mychi2
  
  ! map from serialpresentation to matrix position, given npar^nn-matrix
  function serialtomatrix(ii,nnpar,nn)
    integer, intent(in) :: ii, nnpar, nn
    integer :: serialtomatrix(nnpar)
    integer :: it, remain, num
    
    num=ii-1
    remain=num
    do it = 0, nnpar-1
      ! from decimal to base-nn, where 0 -> 1
      serialtomatrix(it+1)=modulo(remain,nn)
      remain=num/nn
    end do
    serialtomatrix=serialtomatrix+1
  end function serialtomatrix

  ! inverse of serialtomatrix
  function matrixtoserial(ii,nnpar,nn)
    integer, intent(in) :: nnpar, ii(nnpar), nn
    integer :: matrixtoserial
    integer :: it, remain, num
    
    matrixtoserial=1
    do it=1,nnpar
      matrixtoserial=matrixtoserial+nn**(it-1)*(ii(it)-1)
    end do
    
  end function matrixtoserial

end subroutine vdi_GetParsMinimizeChi2
