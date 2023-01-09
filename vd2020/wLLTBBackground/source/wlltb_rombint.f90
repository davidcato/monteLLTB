module wlltb_numtools
  use wlltb_constants
  implicit none

contains

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function wlltb_rombint(f,a,b,tol)
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
			use wlltb_constants
			real(dl), intent(in) :: x
			real(dl) :: f
		end function f
	end interface
  real(dl) :: wlltb_rombint
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
  do
!10        i=i+1
  i=i+1
!  write(*,*)i
  if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
    exit !    go to 40
!  Calculate next trapezoidal rule approximation to integral.
    g0=0._dl
    do k=1,nint
        g0=g0+f(a+(k+k-1)*h)
    end do
    g0=0.5d0*g(1)+h*g0
    h=0.5d0*h
    nint=nint+nint
    jmax=min(i,MAXJ)
    fourj=1._dl
    do j=1,jmax
!  Use Richardson extrapolation.
      fourj=4._dl*fourj
      g1=g0+(g0-g(j))/(fourj-1._dl)
      g(j)=g0
      g0=g1
    end do
    if (abs(g0).gt.tol) then
      error=1._dl-gmax/g0
    else
      error=gmax
    end if
    gmax=g0
    g(jmax+1)=g0
!      go to 10
    end do
!40      voidrombint=g0
    wlltb_rombint=g0
    if (i.gt.MAXITER.and.abs(error).gt.tol)  then
      write(*,*) 'Warning: voidrombint failed to converge; '
      write (*,*)'integral, error, tol:', wlltb_rombint,error, tol
    end if
      
  end function wlltb_rombint

end module wlltb_numtools
