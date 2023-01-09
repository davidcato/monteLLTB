!! Carlson symmetric elliptic integral R_D
!! Using algorithm from arXiv:math/9409227 (Carlson)
!!
!! Written by Wessel Valkenburg, 2010
!!
!! License...


	function rd_c(x,y,z)
		implicit none
		integer, parameter :: pr=dl
		complex(pr) :: rd_c

		include 'rd_declarations_c.f90'
		include 'rd_algorithm.f90'

		rd_c = rrdd
	end function rd_c

	function rd_csp(x,y,z)
		implicit none
		integer, parameter :: pr=sp
		complex(pr) :: rd_csp

		include 'rd_declarations_c.f90'
		include 'rd_algorithm.f90'

		rd_csp = rrdd
	end function rd_csp

	function rd_r(x,y,z)
		implicit none
		integer, parameter :: pr=dl
		real(pr) :: rd_r

		include 'rd_declarations_r.f90'
		include 'rd_algorithm.f90'

		rd_r = rrdd
	end function rd_r

	function rd_rsp(x,y,z)
		implicit none
		integer, parameter :: pr=sp
		real(pr) :: rd_rsp

		include 'rd_declarations_r.f90'
		include 'rd_algorithm.f90'

		rd_rsp = rrdd
	end function rd_rsp
