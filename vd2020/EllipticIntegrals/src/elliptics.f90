module ELLIPTICS
	implicit none

	! Everything double
!	integer, parameter :: ddl=16
	integer, parameter :: dl=kind(1.d0)
	integer, parameter :: sp=4
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	real(dl), parameter :: precision_p(16)=(/0.d0,1.d-4,0.d0,1.d-8,0.d0,0.d0,0.d0,1.d-16,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d-32/)
	real(dl) :: precision(16)=precision_p

	real(dl), parameter :: pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068_dl

	INTERFACE rf
		MODULE PROCEDURE rf_c, rf_csp, rf_r, rf_rsp
	END INTERFACE

	INTERFACE rc
		MODULE PROCEDURE rc_c, rc_csp, rc_r, rc_rsp
	END INTERFACE

	INTERFACE rj
		MODULE PROCEDURE rj_r, rj_rsp, rj_c, rj_csp
	END INTERFACE

	INTERFACE rd
		MODULE PROCEDURE rd_r, rd_rsp, rd_c, rd_csp
	END INTERFACE

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

        Interface elliptics_precisiongoal
           module procedure elliptics_precisiongoal_sl, elliptics_precisiongoal_dl
        end interface

private


public rf, rc, rj, rd, elliptics_precisiongoal

contains

include 'elliptics_tools.f90'

include 'rf.f90'
include 'rc.f90'
include 'rj.f90'
include 'rd.f90'

end module ELLIPTICS
