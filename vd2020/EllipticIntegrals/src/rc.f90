!! Carlson symmetric elliptic integral R_J
!! Using algorithm from arXiv:math/9409227 (Carlson)
!!
!! Written by Wessel Valkenburg, 2010
!!
!! License...


function rc_c(x,y)
	implicit none
	integer, parameter :: pr=dl
	complex(pr) :: rc_c

!----
	include 'rc_declarations_c.f90'
!---- algorithm also contains declarations: no code in between!
	include 'rc_algorithm.f90'
!----

	rc_c = rrcc

end function rc_c

function rc_csp(x,y)
	implicit none
	integer, parameter :: pr=sp
	complex(pr) :: rc_csp

!----
	include 'rc_declarations_c.f90'
!---- algorithm also contains declarations: no code in between!
	include 'rc_algorithm.f90'
!----

	rc_csp = rrcc

end function rc_csp

function rc_r(x,y)
	implicit none
	integer, parameter :: pr=dl
	real(pr) :: rc_r

!----
	include 'rc_declarations_r.f90'
!---- algorithm also contains declarations: no code in between!
	include 'rc_algorithm.f90'
!----

	rc_r = rrcc

end function rc_r

function rc_rsp(x,y)
	implicit none
	integer, parameter :: pr=sp
	real(pr) :: rc_rsp

!----
	include 'rc_declarations_r.f90'
!---- algorithm also contains declarations: no code in between!
	include 'rc_algorithm.f90'
!----

	rc_rsp = rrcc

end function rc_rsp
