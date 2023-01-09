Module lltb_background_aofht_funcs
use lltb_params
use lltb_complex_tools
use lltb_prefuncs
use lltb_precision
use lltb_background_eqs
use lltb_background_ht
use lltb_background_errors
implicit none


contains

	function aoft_ddeltda0_mk_l(vpi,a)
		implicit none
		type(voidparams) :: vpi
		real(dl) :: a
		complex(dl) :: aoft_ddeltda0_mk_l
		!---------!
		complex(dl) :: om, ol, ok
		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)

		aoft_ddeltda0_mk_l = - 0.5_dl * ol * (om + ok * dcmplx(a,0._dl)) **(-1.5_dl) * dcmplx(a,0._dl)**(3.5_dl)

	end function aoft_ddeltda0_mk_l

	function aoft_ddeltda0_m_l(vpi,a)
		implicit none
		type(voidparams) :: vpi
		real(dl) :: a
		complex(dl) :: aoft_ddeltda0_m_l
		!---------!
		complex(dl) :: om, ol, ok
		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)

		aoft_ddeltda0_m_l = - 0.5_dl * ol * om**(-1.5_dl) * dcmplx(a,0._dl)**(3.5_dl)

	end function aoft_ddeltda0_m_l

	function aoft_ddeltda0_m_k(vpi,a)
		implicit none
		type(voidparams) :: vpi
		real(dl) :: a
		complex(dl) :: aoft_ddeltda0_m_k
		!---------!
		complex(dl) :: om, ol, ok
		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)

		aoft_ddeltda0_m_k = - 0.5_dl * ok * om**(-1.5_dl) * dcmplx(a,0._dl)**(1.5_dl)

	end function aoft_ddeltda0_m_k

	subroutine aoft_m(vpi,a, t)
		implicit none
		type(voidparams) :: vpi
		complex(dl) :: a
		real(dl) :: t
		!-------!
		real(dL) :: mya_pre
		real(dl), parameter :: twothird = 2._dl/3._dl

		mya_pre = 1.5_dl * (vpi%xyz(1))**0.5_dl * t

		a=dcmplx(mya_pre**twothird,0._dl)

	end subroutine aoft_m

	subroutine aoft_mk(vpi, a, t)
		implicit none
		type(voidparams) :: vpi
		complex(dl) :: a
		real(dl) :: t
		!-------!
		real(dL) :: mya_pre
		real(dl), parameter :: twothird = 2._dl/3._dl

!		Oops! Doesn't exist..

	end subroutine aoft_mk

end Module lltb_background_aofht_funcs
