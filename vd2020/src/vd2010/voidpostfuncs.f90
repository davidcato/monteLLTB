
  function z_FLRW(r)
    implicit none
    real(dl) :: r, z_FLRW, z
    real(dl) :: tr(1:2), ztop, zbot!, therealzB
    real(dl), parameter :: z_flrw_prec=1.d-8
    integer :: i
    ! Stupidity
    z=NZ_FLRW(r)
    z_FLRW=z
    return
    ! End Stupidity
    zbot=0.d0

    ztop=VPNR_FLRW%ztr(1,size(VPNR_FLRW%ztr(1,:)))
    if(isnan(ztop))then
       if(voidfeedback.gt.0)write(*,*)'z_FLRW rejects this model.'
       VoidError=.true.
       return
    end if
    z=(ztop+zbot)*0.5d0
    i=0
    if(Lbound(tr,1)>1)then
       call VoidErrorReport()
       stop 'It happened: lbound(tr:,1)>1. In z_flrw(r).'
    end if
    tr(2)=0.d0
    do while (abs(tr(2)/r-1.d0).gt. z_flrw_prec)
       call asplint(VPNR_FLRW%ztr(1,:),VPNR_FLRW%ztr(2:3,:),z,tr,VPNR_FLRW%ddztr(2:3,:))
       if(tr(2).gt.r)then
          ztop=z
       else if(tr(2).lt.r)then
          zbot=z
       else
          if(isnan(z).or.any(isnan(tr)))then
             if(voidfeedback.gt.0)write(*,*)'z_FLRW rejects this model.'
             VoidError=.true.
             return
          end if

          if(voidfeedback.gt.0)write(*,'(A,5E17.8)')' #kmax, zB, omk, H0_out, omdmh2',VP%kmax,VP%zB,VP%omk,VP%H0, VP%omdmh2

!          write(*,'(A,3E14.5)')' this r, r_goal, this z:',tr(2),r,z
          call VoidErrorReport()
          stop 'z_flrw fatal error.'
       end if
       z=(ztop+zbot)*0.5d0

       if(ztop.eq.zbot)then
          if(ztop.eq.VPNR_FLRW%ztr(1,size(VPNR_FLRW%ztr(1,:))))then
             write(*,*) 'z_FLRW hanging.'
             VoidError=.true.
             write(*,*)'r(asked),L', r, VP%L
!             stop ' Increase void_zmax for this model. Cannot find z_flrw.'
             z_FLRW = 1.d30 ! perhaps it is only rtopini who is asking.
             return
          end if
          write(*,*) 'z_FLRW hanging.'
          VoidError=.true.
          return
       end if
          
          
       if(i.gt.20000)then
       if(ztop.eq.VPNR_FLRW%ztr(1,size(VPNR_FLRW%ztr(1,:))))then
          write(*,*) 'z_FLRW hanging. ztop==max(VPNR_FLRW).'
          VoidError=.true.
!          stop ' Increase void_zmax for this model. Cannot find z_flrw.'
          return
       end if
          write(*,*)VPNR_FLRW%ztr(1,:),'\n'
          write(*,*)VPNR_FLRW%ztr(2,:),'\n'
          write(*,*)VPNR_FLRW%ztr(3,:),'\n'
          write(*,*)'Distance to goal (relative):',abs(tr(2)/r-1.d0)
          write(*,'(A,4E30.21)')'z, zbot, ztop, zmax:',z,zbot,ztop,VPNR_FLRW%ztr(1,size(VPNR_FLRW%ztr(1,:)))
          write(*,'(A,3E30.21)')'r(goal), t(now), r(now):',r,tr
          call VoidErrorReport()
          stop
       end if
       i=i+1
    end do
    
    z_FLRW=z


  end function z_FLRW


  function VoidDA(z)
!    use camb, only : AngularDiameterDistance
    real(dl), intent(in) :: z
    real(dl) :: voidda
    
!    if(dabs(VP%zB)>1.d-20)then
   VoidDA = DA(z)

  end function VoidDA

  function DA(z,trin)
    use ode_path
    real(dl) :: DA
    real(dL) :: z
    real(dl) :: tr(1:2)
    real(dl), optional :: trin(2)
    real(dl) :: thisda(2)

    if(VoidError)then
       DA=1.d30
       return
    end if
    
    DA=DA_unnorm(z,trin)  * (VPc/1000.d0) / VP%Mtilde_paper * VP%Mtilde 
    return
!!$
!!$    if(present(trin))then
!!$       if(Lbound(tr,1)>1)then
!!$          call VoidErrorReport()
!!$          stop 'It happened: lbound(tr:,1)>1. In DA.'
!!$       end if
!!$       tr=trin
!!$       if(Lbound(tr,1)>1)then
!!$          call VoidErrorReport()
!!$          stop 'It happened: lbound(tr:,1)>1. In DA.'
!!$       end if
!!$    else
!!$       call asplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),z,tr,VPNR%ddztr(2:3,:))
!!$    end if
!!$    if(Lbound(tr,1)>1)then
!!$       call VoidErrorReport()
!!$       stop 'It happened: lbound(tr:,1)>1. In DA.'
!!$    end if
!!$
!!$!    DA = voidR(tr(2),tr(1)) *(VPc/1000.d0) * VP%Hlocal / VP%H0
!!$    ! in Mpc:
!!$!    DA = voidR(tr(2),tr(1)) * (VPc/1000.d0) / VP%Mtilde_paper * VP%Mtilde / voida(VP%LTB_FLRW_r,vp%t0) * voida(0.d0,vp%t0)
!!$    DA = voidR(tr(2),tr(1)) * (VPc/1000.d0) / VP%Mtilde_paper * VP%Mtilde / avoidS(0.d0,vp%t0)
!!$
!!$    ! Wv new: with collapse we can have Hlocal < 0 (Which is not fitted to data anyway)
!!$    DA=dabs(da)
!!$    ! end wv new
!!$
!!$    ! DA in units of Mpc 
!!$    ! the units are equal to those of camb.
!!$
!!$    if(vp%r0/=0.d0)then
!!$       call asplint(VPNR%rtz(1,:),VPNR%rtz(4:5,:),tr(2),thisda,VPNR%ddrtz(4:5,:))
!!$!       DA = thisda(1) * (VPc/1000.d0) / VP%Mtilde_paper * VP%Mtilde / voida(VP%LTB_FLRW_r,vp%t0) * voida(0.d0,vp%t0)
!!$       DA = thisda(1) * (VPc/1000.d0) / VP%Mtilde_paper * VP%Mtilde / avoidS(VP%r0,vp%t0)
!!$    end if

  end function DA

 
  function DA_unnorm(z,trin)
    use ode_path
    real(dl) :: DA_unnorm
    real(dL) :: z
    real(dl) :: tr(1:2)
    real(dl), optional :: trin(2)
    real(dl) :: thisda(2)
    if(VoidError)then
       DA_unnorm=1.d30
       return
    end if
    
    if(present(trin))then
       if(Lbound(tr,1)>1)then
          call VoidErrorReport()
          stop 'It happened: lbound(tr:,1)>1. In DA.'
       end if
       tr=trin
       if(Lbound(tr,1)>1)then
          call VoidErrorReport()
          stop 'It happened: lbound(tr:,1)>1. In DA.'
       end if
    else
       call asplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),z,tr,VPNR%ddztr(2:3,:))
    end if
    if(Lbound(tr,1)>1)then
       call VoidErrorReport()
       stop 'It happened: lbound(tr:,1)>1. In DA.'
    end if

    ! When taking R(r,t) as distance measure, it already is in orthonormal
    ! coordinates of the observer at r0.
    DA_unnorm = voidR(tr(2),tr(1))

    ! Wv new: with collapse we can have Hlocal < 0 (Which is not fitted to data anyway)
    DA_unnorm=dabs(da_unnorm)
    ! end wv new

    ! DA in units of Mpc 
    ! the units are equal to those of camb.

    if(vp%r0/=0.d0)then
       call asplint(VPNR%rtz(1,:),VPNR%rtz(4:5,:),tr(2),thisda,VPNR%ddrtz(4:5,:))
       DA_unnorm = thisda(1)
       ! Like the above, also this term 'thisda' is already in orthonormal
       ! coordinates of the observer: the solid angle Omega is defined
       ! at the observer, and the surface area of the beam is defined to be
       ! zero at the observer. Hence everything is set.
    end if

  end function DA_unnorm

 

  ! This function, NR, is only called in voidsetL
  ! and we really only need r(zB). However, in the
  ! case of several distances corresponding to one
  ! redshift, we need the largest r (blueshift only 
  ! occurs for z < zB).
  function NR_FLRW(z,t)
    use ode_path
    implicit none
    real(dl) :: NR_FLRW
    real(dl), optional :: z, t
    real(dl), pointer :: tr(:,:)
    real(dl) :: zr(2)
    integer :: n
    real(dl), pointer :: ta(:), zra(:,:), ddzra(:,:)


    if(present(z))then
       
       allocate(tr(2,1))
       
!       call masplint(VPNR_FLRW%ztr(1,:),VPNR_FLRW%ztr(2:3,:),z,tr,VPNR_FLRW%ddztr(2:3,:),VPNR_FLRW%kmaxz,VPNR_FLRW%kminz)
    ! new 23/02/2011
    call MultiValuedInvert(VPNR_FLRW,z,tr)

    if(Lbound(tr,1)>1)then
       call VoidErrorReport()
       stop 'It happened: lbound(tr:,1)>1. In NR_FLRW.'
    end if
       
       NR_FLRW = maxval(tr(2,:))
       
       if(NR_FLRW.gt.1.d28)then
!          write(*,*)tr(1,:)
!          write(*,*)tr(2,:)
          call VoidErrorReport()
          stop 'function NR_FLRW(z) has error.'
       end if
       deallocate(tr)
    else if(present(t)) then

       n = size(VPNR_FLRW%ztr(1,:))
       
       if(VPNR_FLRW%ztr(2,n)<VPNR_FLRW%ztr(2,1))then

          ta => VPNR_FLRW%ztr(2,n:1:-1)
          zra => VPNR_FLRW%ztr(1:3:2,n:1:-1)
       else

          ta => VPNR_FLRW%ztr(2,1:n)
          zra => VPNR_FLRW%ztr(1:3:2,1:n)
       end if

       allocate(ddzra(2,size(VPNR_FLRW%ztr(1,:))))
       call aspline(ta, zra, ddzra) 

!!$       if((t>maxval(ta)).or.(t<minval(ta)))then
!!$          write(*,*)'Fuck fuck fuck fuck fuck.'
!!$       end if

       call asplint(ta,zra,t,zr,ddzra)
!       write(*,'(6E15.4)')zr(:), z_FLRW(zr(2)), t, nt_FLRW(zr(1))
       NR_FLRW = zr(2)

       deallocate(ddzra)
       nullify(zra,ta)

       if(NR_FLRW.gt.1.d28)then
!          write(*,*)zr(1)
!          write(*,*)zr(2)
          call VoidErrorReport()
          stop 'function NR_FLRW(t) has error.'
       end if

    end if
    
  end function NR_FLRW

  

  ! This function, NR, is only called in voidsetL
  ! and we really only need r(zB). However, in the
  ! case of several distances corresponding to one
  ! redshift, we need the largest r (blueshift only 
  ! occurs for z < zB).
  function NR(z)
    use ode_path
    real(dl) :: NR, z
    real(dl), pointer :: tr(:,:)
    
    allocate(tr(2,1))

!    call masplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),z,tr,VPNR%ddztr(2:3,:),VPNR%kmaxz,VPNR%kminz)
    ! new 23/02/2011
    call MultiValuedInvert(VPNR,z,tr)
    if(Lbound(tr,1)>1)then
       call VoidErrorReport()
       stop 'It happened: lbound(tr:,1)>1. In NR(z).'
    end if

    NR = maxval(tr(2,:))
    
    if(NR.gt.1.d28)then
!       write(*,*)tr(1,:)
!       write(*,*)tr(2,:)
       call VoidErrorReport()
!       stop 'function NR has error.'
       write(*,*)'function NR has error.'
       VoidError = .true.
       return
    end if
    deallocate(tr)

    

  end function NR

  function vdi_NTGyr(z,r)
    implicit none
    real(dl) :: vdi_NTGyr
    real(dl), intent(in), optional :: z,r

    vdi_NTGyr = NT(z=z,r=r) * 977.8d0

  end function vdi_NTGyr

  ! This function, NR, is only called in voidsetL
  ! and we really only need r(zB). However, in the
  ! case of several distances corresponding to one
  ! redshift, we need the largest r (blueshift only 
  ! occurs for z < zB).
  function NT(z,r)
!    use ode_path
    real(dl) :: NT
    real(dl), optional :: r, z
    real(dl), pointer :: tr(:,:)
    real(dl) :: tz(2)
    integer  :: n
    real(dl), allocatable :: x(:),y(:,:),ddy(:,:) 

    if(present(z))then

       allocate(tr(2,1))

!       call masplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),z,tr,VPNR%ddztr(2:3,:),VPNR%kmaxz,VPNR%kminz)
    ! new 23/02/2011

    call MultiValuedInvert(VPNR,z,tr)

    if(Lbound(tr,1)>1)then
       call VoidErrorReport()
       stop 'It happened: lbound(tr:,1)>1. In NT(z,r).'
    end if
       
       NT = maxval(tr(1,:))
       deallocate(tr)

    else if(present(r)) then
  
       n=size(VPNR%ztr(1,:))
       allocate(x(n))
       allocate(y(2,n))
       allocate(ddy(2,n))
       
       x=VPNR%ztr(3,:)
       y(1,:)=VPNR%ztr(2,:)
       y(2,:)=VPNR%ztr(1,:)
       
       call aspline(x,y,ddy)
       
       call asplint(x,y,r,tz,ddy)
       
       NT = tz(1)

       deallocate(x,y,ddy)
       
    else
       stop 'function NT called without arguments.'
    end if
    
    if(NT.gt.1.d28)then
!       write(*,*)tr(1,:)
!       write(*,*)tr(2,:)
       call VoidErrorReport()
       stop 'NT fucked up.'
    end if
    
    

  end function NT


  function Nz(r,t)
    real(dl) :: NZ
    real(dl), optional :: r, t

    if(present(r))then
       NZ=NZr(r)
    else if (present(t)) then
       NZ=NZt(t)
    else
       stop 'function nz called with wrong arguments? voidpostfuncs.f90'
    end if

  end function NZ

  function Nzr(r)
    use ode_path
    real(dl) :: NZr, r!,z
    real(dl) :: tz(2)
    integer  :: n
    real(dl), allocatable :: x(:),y(:,:),ddy(:,:) 
    n=size(VPNR%ztr(1,:))
    allocate(x(n))
    allocate(y(2,n))
    allocate(ddy(2,n))

!!$    x=VPNR%ztr(3,:)
!!$    y(1,:)=VPNR%ztr(2,:)
!!$    y(2,:)=VPNR%ztr(1,:)
!!$
!!$    call aspline(x,y,ddy)
!!$
!!$    call asplint(x,y,r,tz,ddy)

    x=VPNR%rtz(1,:)
    y=VPNR%rtz(2:3,:)
    call asplint(x,y,r,tz,VPNR%ddrtz)

    Nzr = tz(2)

    deallocate(x,y,ddy)

  end function NZr

  function Nzt(t)
    use ode_path
    real(dl) :: NZt, t!,z
    real(dl) :: rz(2)
    integer  :: n
    real(dl), pointer :: x(:),y(:,:),ddy(:,:) 
!    real(dl), allocatable :: x(:),y(:,:),ddy(:,:) 
    n=size(VPNR%ztr(1,:))
!    allocate(x(n))
!    allocate(y(2,n))
    allocate(ddy(2,n))

    if ( (t<minval(VPNR%ztr(2,:))) .or. (t>maxval(VPNR%ztr(2,:))) )then
       write(0,*) 't outside of results in nzt (voidpostfuncs.f90)'
       deallocate(ddy)
       write(0,*)ddy(1,1)
       stop 't outside of results in nzt (voidpostfuncs.f90)'       
    end if

    x=>VPNR%ztr(2,n:1:-1)
    y=>VPNR%ztr(3:1:-2,n:1:-1)


    call aspline(x,y,ddy)

    call asplint(x,y,t,rz,ddy)

    Nzt = rz(2)
    deallocate(ddy)
    nullify(x,y)

  end function NZt

  function nz_flrw(r,t)
    implicit none
    real(dl), intent(in), optional :: r,t
    real(dl) :: nz_FLRW
    
    if(present(r))then
       nz_flrw = nz_flrwr(r)
    else if (present(t))then
       nz_flrw = nz_flrwt(t)
    else
       stop 'NZ_Flrw called without arguments?'
    end if

  end function nz_flrw

  function Nz_FLRWr(r)
    implicit none
!    use ode_path
    real(dl) :: NZ_FLRWr, r!,z
    real(dl) :: tz(2)
    integer  :: n
    real(dl), allocatable :: x(:),y(:,:),ddy(:,:) 
    n=size(VPNR_FLRW%ztr(1,:))
    allocate(x(n))
    allocate(y(2,n))
    allocate(ddy(2,n))

    x=VPNR_FLRW%ztr(3,:)
    y(1,:)=VPNR_FLRW%ztr(2,:)
    y(2,:)=VPNR_FLRW%ztr(1,:)

    call aspline(x,y,ddy)

    call asplint(x,y,r,tz,ddy)

    Nz_FLRWr = tz(2)

    deallocate(x,y,ddy)

  end function NZ_FLRWr

  function Nz_FLRWt(t)
    implicit none
!    use ode_path
    real(dl) :: NZ_FLRWt, t!,z
    real(dl) :: zr(2)
    integer  :: n
    real(dl), pointer :: x(:),y(:,:)
    real(dl), allocatable :: ddy(:,:) 
    n=size(VPNR_FLRW%ztr(1,:))
!    allocate(x(n))
!    allocate(y(2,n))
    allocate(ddy(2,n))

    if(VPNR_FLRW%ztr(2,n)<VPNR_FLRW%ztr(2,1))then
       x => VPNR_FLRW%ztr(2,n:1:-1)
       y => VPNR_FLRW%ztr(1:3:2,n:1:-1)
    else
       x => VPNR_FLRW%ztr(2,1:n)
       y => VPNR_FLRW%ztr(1:3:2,1:n)
    end if
    
    call aspline(x,y,ddy)

    call asplint(x,y,t,zr,ddy)

    Nz_FLRWt = zr(1)

    deallocate(ddy)
    nullify(x,y)

  end function NZ_FLRWt

!!$  function Nz_FLRWt(t)
!!$!    use ode_path
!!$    real(dl) :: NZ_FLRWt,t!,z
!!$    real(dl) :: rz(2)
!!$    integer  :: n
!!$    real(dl), allocatable :: x(:),y(:,:),ddy(:,:) 
!!$    n=size(VPNR_FLRW%ztr(1,:))
!!$    allocate(x(n))
!!$    allocate(y(2,n))
!!$    allocate(ddy(2,n))
!!$
!!$    x=-VPNR_FLRW%ztr(2,:)
!!$    y(1,:)=VPNR_FLRW%ztr(3,:)
!!$    y(2,:)=VPNR_FLRW%ztr(1,:)
!!$
!!$    call aspline(x,y,ddy)
!!$
!!$    call asplint(x,y,-t,rz,ddy)
!!$
!!$    Nz_FLRWt = rz(2)
!!$
!!$    deallocate(x,y,ddy)
!!$
!!$  end function NZ_FLRWt

  function NT_FLRW(z, r)
    use ode_path
    real(dl) :: NT_FLRW
    real(dl), optional :: z, r
    real(dl), pointer :: tr(:,:), ddy(:,:)
    real(dl) :: zt(2)
    integer :: n, i


    if(present(z))then
       allocate(tr(2,1))
       
!       call masplint(VPNR_FLRW%ztr(1,:),VPNR_FLRW%ztr(2:3,:),z,tr,VPNR_FLRW%ddztr(2:3,:),VPNR_FLRW%kmaxz,VPNR_FLRW%kminz)
     ! new 23/02/2011
    call MultiValuedInvert(VPNR_FLRW,z,tr)
      
    if(Lbound(tr,1)>1)then
       call VoidErrorReport()
       stop 'It happened: lbound(tr:,1)>1. In NT_FLRW(z,r).'
    end if
       NT_FLRW = maxval(tr(1,:))
       
       if(NT_FLRW.gt.1.d28)then
!          write(*,*)tr(1,:)
!          write(*,*)tr(2,:)
          call VoidErrorReport()
          stop 'NT_FLRW(z) fucked up.'
       end if
       deallocate(tr)
    else if(present(r)) then
       
       n = size(VPNR_FLRW%ztr(3,:))

!       allocate(ddy(2,n))
!       call aspline(VPNR_FLRW%ztr(3,:),VPNR_FLRW%ztr(1:2,:),ddy)
       
!       call asplint(VPNR_FLRW%ztr(3,:),VPNR_FLRW%ztr(1:2,:),r,zt,ddy)
       call asplint(VPNR_FLRW%rtz(1,:),VPNR_FLRW%rtz(3:2:-1,:),r,zt,VPNR_FLRW%ddrtz(3:2:-1,:))
       
       NT_FLRW = zt(2)
       
       if(NT_FLRW.gt.1.d28)then
!          write(*,*)tr(1,:)
!          write(*,*)tr(2,:)
          call VoidErrorReport()
          write(*,*)n, r, maxval(VPNR%ztr(3,:)), VoidError
          do i = 1, n
             write(123,*)VPNR%ztr(:,i)
          end do
          stop 'NT_FLRW(z) fucked up.'
       end if


    end if
    

  end function NT_FLRW

  function VoidComovLTBSize()
    implicit none
    real(dl) :: VoidComovLTBSize

    VoidComovLTBSize = VoidComovDistMPC(0._dl,VP%L,VP%t0,VP%r0)
  end function VoidComovLTBSize
  
  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        function voidrombint(f,a,b,tol)
  !      use Precision
!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, parameter :: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
!        real(dl) f
!        external f
	interface
		function f(x)
			use precision
			real(dl), intent(in) :: x
			real(dl) :: f
		end function f
	end interface
        real(dl) :: voidrombint
        real(dl), intent(in) :: a,b,tol
        integer :: nint, i, k, jmax, j
        real(dl) :: h, gmax, error, g, g0, g1, fourj
!

        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
            go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0._dl
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1._dl
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4._dl*fourj
            g1=g0+(g0-g(j))/(fourj-1._dl)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1._dl-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      voidrombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)  then
          write(*,*) 'Warning: voidrombint failed to converge; '
          write (*,*)'integral, error, tol:', voidrombint,error, tol
        end if
        
        end function voidrombint

  function VoidComovDistMPC(rmin,rmax,t,robs)
    implicit none
    real(dl), intent(in) :: rmin,rmax,t
    real(dl), intent(in), optional :: robs
    real(dl) :: VoidComovDistMPC
!    real(dl), external :: voidrombint
    REAL(DL) :: testl, myrobs	
!	external :: MyIntegrand

    if(present(robs))then
       myrobs=abs(robs)
    else
       myrobs = abs(VP%r0)
    end if
	
	if(VP%is_pure_flrw)then
		VoidComovDistMPC=0.d0
		return
	end	if

    testl=voidrombint(rmin,rmax,VoidDistPrec)
    
    ! At this point, testl is scaled in the orthonormal 
    ! coordinates of an observer who is at r* for which
    ! sqrt(g_{rr}(r*)) = 1.

!    testl=testl/ avoidS(myrobs,t)* VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde 
  !this is wrong. s=int g_rr dr = a scalar. In the orthonormal frame of the
  !observer, s = s. y^{r hat} has norm s.

    ! And now testl is in the coordinates of observer at r_obs
    ! since we rescaled the metric to his orthonormal coordinates.
    testl = testl * VPc*1.d-3 / VP%Mtilde_paper * VP%Mtilde 
    VoidComovDistMPC = testl 

  contains 

    function f(r)
      real(dl) :: r
      real(dl) :: f
	       
      f = aVoidS(r,t) 
      if(f<0.d0)then
        VoidError = .true.
        if(voidfeedback>2)write(0,*)'Shell crossing in VoidComovDistMPC'
      end if
        
    end function f
	
	include 'vd2010/voidrombint.f90'
	
  end function VoidComovDistMPC

  
  function VoidzB()
    implicit none
    real(dl) VoidzB
    VoidzB = VP%zB
  end function VoidzB

 

  function vdii_H0avv(zmin,zmax)
    implicit none
    real(dl), intent(in) :: zmin,zmax
    real(dl) :: vdii_H0avv
    REAL(DL) :: testl

    if(VoidError)then
      vdii_H0avv=0.d0
      return
    end if

    testl=voidrombint(zmin,zmax,VoidDistPrec)
    
    vdii_H0avv = testl / (zmax-zmin)

  contains 

    function f(z)
      real(dl) :: z
      real(dl) :: f
	    real(dl) :: tr(2)   
         
      call asplint(VPNR_FLRW%ztr(1,:),VPNR_FLRW%ztr(2:3,:),z,tr,VPNR_FLRW%ddztr(2:3,:))
      call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=tr(2),t=tr(1),H=f)

      if(isnan(f))then
        writE(0,*)'isnan f in vdii_H0avv'
        call VoidErrorReport()
        stop
      end if
        
    end function f
	
	include 'vd2010/voidrombint.f90'
	
  end function vdii_H0avv

  
  

  
  function vdi_HomogeneityLnLike()
    implicit none
    real(dl) :: vdi_HomogeneityLnLike

    if (VP%is_pure_flrw) then
      vdi_HomogeneityLnLike = 0.d0
    else if (VoidError) then
      vdi_HomogeneityLnLike = 1.d30
    else
      if(VoidDensSigmaIsSet) then
        vdi_HomogeneityLnLike = VP%dens_nsigma_lnlike
      else
      vdi_HomogeneityLnLike = 0.d0
        write(0,*)'trying to get vdi_HomogeneityLnLike, but is not set yet.'
      call VoidErrorReport()
!        stop
      end if
    end if
  end function vdi_HomogeneityLnLike


