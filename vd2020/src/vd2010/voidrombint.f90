!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Rombint without the dummy interface: include this fortran file when function f(r) actually exists.
        function voidrombint(a,b,tol)
        use Precision
!  voidrombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be real(dl) and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        implicit none
        integer, parameter :: MAXITER=20
        integer, parameter :: MAXJ=5
        dimension g(MAXJ+1)
!!$	interface
!!$		function f(x)	
!!$			use precision
!!$			real(dl) :: f, x
!!$		end function f
!!$	end interface
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
            if(VoidError)i=i+MAXITER ! NEW VOIDERROR
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
          ! write(*,*)'David suspects that it is a numerical error.'
          write(*,*)'Rejecting model due to numerical trouble.'
          writE(*,'(A,5E17.8)')' #kmax, zB, omk, H0_out, omdmh2'
          call EXIT(voidfeedback - voidfeedback)
        end if
        
        end function voidrombint
