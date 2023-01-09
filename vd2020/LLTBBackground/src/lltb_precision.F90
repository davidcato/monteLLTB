module lltb_precision

    integer, parameter :: dl=kind(1.d0)
#ifdef NOQUAD
    integer, parameter :: ql=dl
#else
    integer, parameter :: ql=2*dl
#endif
	integer, parameter :: sp=kind(1.e0)
 	real(dl), parameter :: pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068_dl
!	real(dl), parameter :: accuracy_p(16)=(/0.e0_dl,1.e-4_dl,0.e0_dl,1.e-8_dl,0.e0_dl,0.e0_dl,0.e0_dl,1.e-16_dl,0.e0_dl,0.e0_dl,0.e0_dl,0.e0_dl,0.e0_dl,0.e0_dl,0.e0_dl,1.e-32_dl/)
	real(dl), parameter :: accuracy_p(16)=(/1.e30_dl,1.e-4_dl,1.e30_dl,1.e-8_dl,1.e30_dl,1.e30_dl,1.e30_dl,1.e-16_dl,1.e30_dl,1.e30_dl,1.e30_dl,1.e30_dl,1.e30_dl,1.e30_dl,1.e30_dl,1.e-32_dl/)
	real(dl), parameter :: lin_expansion_p(16)=1._dl/sqrt(accuracy_p) ! limit for linear expansions. epsilon =< Sqrt(acc) -> error in expansion is epsilon^2 = accuracy

	real(dl), parameter :: acceptance_pure_p = 1.e2_dl ! ! amount of imaginary numerical residue that we allow for in normal arithmetic, over accuracy.

	real(dl), parameter :: acceptance_p=1.e6_dl ! amount of imaginary numerical residue that we allow for in special functions, over accuracy.
	real(dl), parameter :: acceptance_pole_p=1.e10_dl ! amount of imaginary numerical residue allowed close to a = {u|v|w}.

	real(dl), parameter :: tiny = 1.d-300
        real(dl), parameter :: infinity=1.d300

        real(dl), parameter :: accuracy_downscale_p = 6._dl ! downscale is amount of loss we allow for if i > imax in backgroundaofHt

	real(dl), parameter :: rel_acceptance_p(3) = lin_expansion_p(dl) ! accuracy*rel_acceptance gives limit under which we expand in one of the parameters

!	real(dl), parameter :: omk_limit = 1.e-8 ! lower limit of curvature k: below this assuming k=0 is more accurate.
!	real(dl), parameter :: omm_limit = 1.e0 ! lower limit of omega_matter: below this assuming omega_m=0 is more accurate.
!	real(dl), parameter :: oml_limit = 1.e-2l ! lower limit of omega_lambda: below this assuming lambda=0 is more accurate.
!	real(dl), parameter :: a_limit = 1.e-15_dl ! lower limit for a: below this limit we stop rescaling the parameters, since a vanishes otherwise.

        ! Make all these adjustable parameters for speed:

	real(dl) :: accuracy(16)=accuracy_p
	real(dl) :: lin_expansion(16)=lin_expansion_p

	real(dl) :: acceptance_pure = acceptance_pure_p

	real(dl) :: acceptance=acceptance_p
	real(dl) :: acceptance_pole=acceptance_pole_p

        real(dl) :: accuracy_downscale = accuracy_downscale_p


	real(dl) :: rel_acceptance(3) = rel_acceptance_p ! accuracy*rel_acceptance gives limit under which we expand in one of the parameters
	! or set individually:
!	real(dl), parameter :: rel_acceptance(3) = (/1.d8,1.d8,1.d8/) ! accuracy*rel_acceptance gives limit under which we expand in one of the parameters



end module lltb_precision
