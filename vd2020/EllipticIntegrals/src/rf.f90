!! Carlson symmetric elliptic integral R_F
!! Using algorithm from arXiv:math/9409227 (Carlson)
!!
!! Written by Wessel Valkenburg, 2010
!!
!! License...


	function rf_c(x,y,z)
		implicit none
		integer, parameter :: pr=dl
		complex(pr) :: rf_c

		include 'rf_declarations_c.f90'
		include 'rf_algorithm.f90'

		rf_c = rrff
	end function rf_c

	function rf_csp(x,y,z)
		implicit none
		integer, parameter :: pr=sp
		complex(pr) :: rf_csp

		include 'rf_declarations_c.f90'
		include 'rf_algorithm.f90'

		rf_csp = rrff
	end function rf_csp

	function rf_r(x,y,z)
		implicit none
		integer, parameter :: pr=dl
		real(pr) :: rf_r

		include 'rf_declarations_r.f90'
		include 'rf_algorithm.f90'

		rf_r = rrff
	end function rf_r

	function rf_rsp(x,y,z)
		implicit none
		integer, parameter :: pr=sp
		real(pr) :: rf_rsp

		include 'rf_declarations_r.f90'
		include 'rf_algorithm.f90'

		rf_rsp = rrff
	end function rf_rsp
