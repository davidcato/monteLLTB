module lltb_params
use lltb_precision

! The main information of this routine:
	type voidparams
		real(dl) :: k =0._dl
		real(dl) :: mtilde
		real(dl) :: lambda
		real(dl) :: xyz(3), workxyz(3), workxyznonzero(3) ! xyz(1) = Omega_m, xyz(2) = Omega_k, xyz(3) = Omega_Lambda
		real(dl) :: workrescalefac, workxyzfac
		complex(dl) :: uvw(3) ! uvw = roots of [xyz(1) + xyz(2) a + xyz(3) a^3 == 0] in a.
		complex(dl) :: uvw_linear(3) ! roots used for linear expansion
		integer :: posrealroots ! number of positive real roots
		real(dl) :: smallestrealroot, smallestlinroot
                integer :: smallestrootindex, smallestlinrootindex
                real(dl) :: tturn
		logical :: crossedbranchcut
		character(len=3) :: ExpCase, ExpCaseInv ! (L)arge, (S)mall, (Z)ero: "LSZ","LLL","LSS" etc.
	end type voidparams
        
! The information exchanged between cosmomc and lltb_background
        type LB_MYMEM
           real(dl) :: H0_inf, Lambda, Mt, r, t, k
           real(dl) :: t0, Mt2, Lambdaover3, H0_inf2
           real(dl) :: Rltb, Rp, Rpd, S, Sd, a, H, ap, tturn, Hp, apd
           real(dl) :: Rpdd, Sdd, apdd, add
        end type LB_MYMEM
        
	integer :: lfeedback = 0

 
	integer, parameter :: xyz_apower(3) = (/0,1,3/)


end module lltb_params
