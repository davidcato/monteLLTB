module lltb_complex_tools
use lltb_precision
implicit none

	interface maimag ! my interface to make aimag compatible with reals and doubles
		module procedure maimag_csp,maimag_c, maimag_r, maimag_rsp
	end interface maimag

	interface mreal ! my interface to make aimag compatible with reals and doubles
		module procedure mreal_csp,mreal_c, mreal_rsp, mreal_r
	end interface mreal

	interface killreal
		module procedure killreal_csp,killreal_c, killreal_rsp, killreal_r
	end interface killreal

	interface killimag
		module procedure killimag_csp,killimag_c, killimag_r, killimag_rsp
	end interface killimag


contains

	function maimag_rsp(xin)
		implicit none
		real(sp) :: xin, maimag_rsp

		maimag_rsp = 0._sp

	end function maimag_rsp

	function maimag_r(xin)
		implicit none
		real(dl) :: xin, maimag_r

		maimag_r = 0._dl

	end function maimag_r

	function maimag_csp(xin)
		implicit none
		complex(sp) :: xin
		real(sp) :: maimag_csp

		maimag_csp = aimag(xin)

	end function maimag_csp

	function maimag_c(xin)
		implicit none
		complex(dl) :: xin
		real(dl) :: maimag_c
		maimag_c = real(cmplx(0.d0,-0.5d0)*(xin-conjg(xin)))
	end function maimag_c

	function mreal_rsp(xin)
		implicit none
		real(sp) :: xin, mreal_rsp

		mreal_rsp = xin

	end function mreal_rsp

	function mreal_r(xin)
		implicit none
		real(dl) :: xin, mreal_r

		mreal_r = xin

	end function mreal_r

	function mreal_csp(xin)
		implicit none
		complex(sp) :: xin
		real(sp) :: mreal_csp

		mreal_csp = xin-cmplx(0.e0,maimag(xin))

	end function mreal_csp

	function mreal_c(xin)
		implicit none
		complex(dl) :: xin
		real(dl) :: mreal_c

		mreal_c = xin-dcmplx(0.d0,maimag(xin))

	end function mreal_c

	function killreal_r(xin)
		implicit none
		real(dl) :: killreal_r, xin
		killreal_r = 0.d0
	end function killreal_r

	function killreal_rsp(xin)
		real(sp) :: killreal_rsp, xin
		killreal_rsp = 0.d0
	end function killreal_rsp

	function killreal_c(xin)
		complex(dl) :: killreal_c, xin
		killreal_c = xin - mreal(xin)
	end function killreal_c

	function killreal_csp(xin)
		complex(sp) :: killreal_csp, xin
		killreal_csp = xin - mreal(xin)
	end function killreal_csp

	function killimag_r(xin)
		real(dl) :: killimag_r, xin
		killimag_r = xin
	end function killimag_r

	function killimag_rsp(xin)
		real(sp) :: killimag_rsp, xin
		killimag_rsp = xin
	end function killimag_rsp

	function killimag_c(xin)
		complex(dl) :: killimag_c, xin
		killimag_c = xin - dcmplx(0._dl,maimag(xin))
	end function killimag_c

	function killimag_csp(xin)
		complex(sp) :: killimag_csp, xin
		killimag_csp = xin - cmplx(0._sp,maimag(xin))
	end function killimag_csp



	function msqrt(xx)
		complex(dl) :: msqrt, xx
		real(dl) :: tanin, tanout, root, newnorm, newmod, atanin, atanout

		tanin = aimag(xx)/real(xx)

		msqrt=sqrt(xx)
		tanout = aimag(msqrt)/real(msqrt)

		! which quadrant?
		atanout=atan(tanout)+ (0.5d0 - real(msqrt)/abs(real(msqrt)*2.d0))*pi
		if(atanout>pi)atanout=atanout-2.d0*pi
		if(atanout<-pi)atanout=atanout+2.d0*pi

		atanin=atan(tanin)+ (0.5d0 - real(xx)/abs(real(xx)*2.d0))*pi
		if(atanin>pi)atanin=atanin-2.d0*pi
		if(atanin<-pi)atanin=atanin+2.d0*pi

		root=(atanout)/(atanin)

		if(abs(root/0.5d0-1)>1.d-4)then
			write(*,*)'root:',root
			write(*,*)'atanin:',atanin
			write(*,*)'atanout:',atanout
		   !			stop 'wrong root'
		end if


	end function msqrt

	subroutine SortComplexes(xin)
		implicit none
		complex(dl) :: xin(:), myx1(size(xin)), myx2(size(xin)), myx3(size(xin)), dummy
		integer :: i, n, nreals, m, p, foundpairs, mem1,mem2
		logical :: isreal(size(xin)), notconjugate(size(xin))

! first bluntly sort on Real values!
		n = size(xin)

		mem2=n-1
		mem1=0
		do i = 1,n-1
			do m = 1, n-1 ! slow bubble sort
				if(real(xin(m))<real(xin(m+1)))then
					dummy = xin(m)
					xin(m)=xin(m+1)
					xin(m+1)=dummy
					mem1=m
				end if
				if(m>mem2)exit
			end do
			mem2=mem1
		end do
		isreal=.false.
		nreals = 0
		do i = 1, n
!			if( abs( abs(maimag(xin(i)))/abs((xin(i)))) < accuracy(dl) ) then
!			if(  abs(maimag(xin(i))) < accuracy(dl) ) then
			if(  abs(maimag(xin(i)))/abs(mreal(xin(i))) < accuracy(dl) *acceptance_pure) then
				 isreal(i) = .true.
				 xin(i) =  killimag(xin(i))
				 nreals = nreals + 1
			end if
		end	do

		myx1 = xin
		myx2 = xin

		m = 0
		p = 0
		do i = 1, n
			if(isreal(i))then
				m=m+1
				myx1(m)=xin(i)
			else
				p=p+1
				myx2(p)=xin(i)
			end if
		end do

		if(p/=n-nreals)stop 'bug in SortComplex.'

		m=1
		if(nreals>0)then
			xin=myx1(m:nreals)
			m=nreals+1
		end if
		if(p>0)then
				xin(m:n)=myx2(1:p)
		end if

	end subroutine SortComplexes

	function masinh(x)
		implicit none
		complex(dl), intent(in) :: x
		complex(dl) :: masinh

		masinh = log( x + sqrt(dcmplx(1.,0.) + x**2 ) )

	end function masinh


	function get_nearest_cr(arr_c,x_r)
		complex(dl) get_nearest_cr, my_var
		complex(dl), optional :: arr_c(:)
		real(dl), optional :: x_r
		integer :: n,i

		n = size(arr_c)
		my_var = arr_c(1)
		do i = 2, n
			if ( abs(mreal(my_var) - x_r) > abs(mreal(arr_c(i)) - x_r) ) my_var = arr_c(i)
		end do

		get_nearest_cr = my_var

	end function get_nearest_cr


	function Get_Smallest_posrealval(arr_c,index)
		implicit none
		complex(dl) :: arr_c(:)
		real(dl) :: Get_Smallest_posrealval
		real(dl) :: smallestposrealroot
		integer :: i, n, posrealroots, thisindex
                integer, optional :: index


		n=size(arr_c)
                
                thisindex=-1
		smallestposrealroot = 1.e30_dl
			posrealroots=0
			do i = 1, n
				if (abs( abs(aimag(arr_c(i)))/abs((arr_c(i)))) < accuracy(dl) ) then ! is real, and imaginary part should be killed in sortcomplexes already.
					if(mreal(arr_c(i))>0._dl)then
						posrealroots = posrealroots+1
						if(mreal(arr_c(i))< smallestposrealroot)then
                                                   smallestposrealroot=mreal(arr_c(i))
                                                   thisindex=i
                                                end if
					end if
				end if
			end do
		Get_Smallest_posrealval = smallestposrealroot
                if(present(index))index=thisindex
	end function Get_Smallest_posrealval

end module lltb_complex_tools

