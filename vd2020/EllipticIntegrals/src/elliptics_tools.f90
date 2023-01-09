	function maimag_rsp(xin)
		real(sp) :: xin, maimag_rsp

		maimag_rsp = 0._sp

	end function maimag_rsp

	function maimag_r(xin)
		real(dl) :: xin, maimag_r

		maimag_r = 0._dl

	end function maimag_r

	function maimag_csp(xin)
		complex(sp) :: xin
		real(sp) :: maimag_csp

		maimag_csp = aimag(xin)

	end function maimag_csp

	function maimag_c(xin)
		complex(dl) :: xin
		real(dl) :: maimag_c
		maimag_c = real(cmplx(0.d0,-0.5d0)*(xin-conjg(xin)))
	end function maimag_c

	function mreal_rsp(xin)
		real(sp) :: xin, mreal_rsp

		mreal_rsp = xin

	end function mreal_rsp

	function mreal_r(xin)
		real(dl) :: xin, mreal_r

		mreal_r = xin

	end function mreal_r

	function mreal_csp(xin)
		complex(sp) :: xin
		real(sp) :: mreal_csp

		mreal_csp = xin-cmplx(0.e0,maimag(xin))

	end function mreal_csp

	function mreal_c(xin)
		complex(dl) :: xin
		real(dl) :: mreal_c

		mreal_c = xin-dcmplx(0.d0,maimag(xin))

	end function mreal_c

	function killreal_r(xin)
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

        subroutine elliptics_precisiongoal_sl(goal)
          real, intent(in) :: goal

          call elliptics_precisiongoal_dl(dble(goal))

        end subroutine elliptics_precisiongoal_sl

        subroutine elliptics_precisiongoal_dl(goal)
          real(dl), intent(in) :: goal
          real(dl) :: orig_goal
          real(dl) :: power
          ! We initialize with precision goal = 1.d-16
          
          orig_goal = precision_p(dl) ! parameter accuracy
          
          if(goal < orig_goal)then
             write(*,'(A,ES14.5,A,ES14.5,A)')'You demand precisiongoal ',goal,' which is smaller than the hardcoded lowerlimit of ',orig_goal,'. That will not work.'
             stop 'Too high precisiongoal demanded. Elliptics.'
          end if

          power = log(goal)/log(orig_goal)

          precision=precision_p**power

        end subroutine elliptics_precisiongoal_dl
