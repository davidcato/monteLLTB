!! Carlson symmetric elliptic integral R_J
!! Using algorithm from arXiv:math/9409227 (Carlson)
!!
!! Written by Wessel Valkenburg, 2010
!!
!! License...


function rj_c(x,y,z,p)
	implicit none
	integer, parameter :: pr=dl
	complex(pr) :: rj_c

!----
	include 'rj_declarations_c.f90'
!---- algorithm also contains declarations: no code in between!
	include 'rj_algorithm.f90'
!----

	rj_c = rrjj

end function rj_c

function rj_csp(x,y,z,p)
	implicit none
	integer, parameter :: pr=sp
	complex(pr) :: rj_csp

!----
	include 'rj_declarations_c.f90'
!---- algorithm also contains declarations: no code in between!
	include 'rj_algorithm.f90'
!----

	rj_csp = rrjj

end function rj_csp

function rj_r(x,y,z,p)
	implicit none
	integer, parameter :: pr=dl
	real(pr) :: rj_r

!----
	include 'rj_declarations_r.f90'
!---- algorithm also contains declarations: no code in between!
	include 'rj_algorithm.f90'
!----

	rj_r = rrjj

end function rj_r

function rj_rsp(x,y,z,p)
	implicit none
	integer, parameter :: pr=sp
	real(pr) :: rj_rsp

!----
	include 'rj_declarations_r.f90'
!---- algorithm also contains declarations: no code in between!
	include 'rj_algorithm.f90'
!----

	rj_rsp = rrjj

end function rj_rsp
