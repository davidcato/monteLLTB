module mod_vd2020
    implicit none

    integer, parameter :: dl = kind(1.d0)
    real(dl), parameter :: pi = 4.d0 * datan(1.d0)

    real(dl), parameter :: VoidDistPrec=1.d-5!  3
    real(dl), parameter :: VoidMinPrec = 1.d-12 ! minimum for voiddistprec: boosted in case of bad VoidSetL
    real(dl), parameter :: VPc = 2.99792458d8 ! voiddistprec => nreps now.
    real(dl), parameter :: u_precision=1.d-9!12!7
    real(dl), parameter :: numerical_derivatives_precision=0.d0!1.d-7!30
    real(dl) :: void_zmax=50.d0
    real(dl), parameter :: void_zmax_fixed=2100.d0 ! 50.d0 ! Changed to 1100 to autimatically check that we can reach last scattering without the metric breaking up (.e.g S going imaginary, curvature to large).
    ! changed zmax to 2100 to be tolerant for extreme z_drag from bao. If zdrag > zmax, reject model.
    real(dl), parameter :: VoidZBPrecision = 1.d-5
    real(dl), parameter :: kmax_max = 1.d3
    real(dl), parameter :: max_double = 1.d307
    real(dl), target :: VoidInfinity = 1.d30, VoidZero = 0.d0

    type, public :: FeedBackSettings
        logical :: parameters = .true.
        logical :: wlltb = .true.
        logical :: lltb = .true.
    end type FeedBackSettings
    type(FeedBackSettings), parameter :: FeedBackDefaults = FeedBackSettings()

    ! the input
    type, public :: VoidParameters
        integer :: profile_function
        real(dl) :: L ! the coordinate distance to the boundary. Mpc.
        real(dl) :: kmax ! The curvature radius at r=0, the maximum value of k(r).
        real(dl) :: w ! the (homogeneous!) Dark Energy equation of state.
        real(dl) :: H0_outside ! the outside FLRW expansion rate. Effectively the age of the universe. km/s/Mpc. Or 1/Mpc. I forgot. Probably 1/Mpc.
        real(dl) :: omm ! Omega_matter: the FLRW dust content
        real(dl) :: omk ! Omega_k: the FLRW curvature content. Mind the sign: omk = 1 - omm - omv
    end type VoidParameters

    type, public :: VoidDistances
        type(VoidParameters) :: params
        type(FeedBackSettings) :: feedback

        ! own parameters, derived from params
        real(dl) :: Mtilde ! mass scale
        real(dl) :: kb_om ! FLRW curvature radius
        real(dl) :: Lambda ! Cosmological constant == 3 H_0^2 Omega_Lambda
        real(dl) :: Hlocal ! \dot a(r, t) / a(r, t) |_{r = 0}. Local expansion rate. Not really related to the observed expansion rate.
        real(dl) :: delta0 ! the energy density at the center, relative to the outer FLRW. Dimensionless.
        real(dl) :: t0 ! the age of the FLRW universe
        real(dl) :: tbar ! forgot what this is.

    contains
        procedure :: init => VDInit
        procedure :: set_derived => VDSetDerived

    end type VoidDistances

contains


    ! Constructor: Set the parameters and compute the derived parameters.
    subroutine VDInit(this, VP, in_feedback)
        implicit none
        class(VoidDistances) :: this
        type(VoidParameters) :: VP
        type(FeedBackSettings), optional :: in_feedback

        this%params = VP

        if ( present(in_feedback) ) then
            this%feedback = in_feedback
        end if

        call this%set_derived()

    end subroutine VDInit

    ! Called by the constructor: compute the derived parameters.
    subroutine VDSetDerived(this)
        use wlltb, only : lltb_normalization
        implicit none

        class(VoidDistances) :: this

        this%kb_om = 4.d0 * this%params%omk / (3.d0 * pi * this%params%omm )

        this%Lambda = 3.d0 * this%params%H0_outside**2 * ( 1 - this%params%omk - this%params%omm)

        call lltb_normalization(this%params%H0_outside, this%kb_om, this%params%kmax, this%Lambda, this%params%w, & ! output:
            this%Mtilde, this%Hlocal, this%delta0, this%t0, this%tbar)

        if ( this%feedback%parameters ) then
            write(*,*) "**************************************"
            write(*,*) "VoidParameters, input and inferred:"
            write(*,'(A20, ES15.6)') "L: ", this%params%L
            write(*,'(A20, ES15.6)') "k_max: ", this%params%kmax
            write(*,'(A20, ES15.6)') "H0(outside): ", this%params%H0_outside
            write(*,'(A20, ES15.6)') "Omega_k: ", this%params%omk
            write(*,'(A20, ES15.6)') "Omega_matter: ", this%params%omm
            write(*,'(A20, ES15.6)') "Omega_DE: ", (1 - this%params%omm - this%params%omk)
            write(*,'(A20, ES15.6)') "w_DE: ", this%params%w
            write(*,'(A20, ES15.6)') "--- derived: "
            write(*,'(A20, ES15.6)') "M_tilde: ", this%Mtilde
            write(*,'(A20, ES15.6)') "Hlocal: ", this%Hlocal
            write(*,'(A20, ES15.6)') "delta0: ", this%delta0
            write(*,'(A20, ES15.6)') "t0: ", this%t0
            write(*,'(A20, ES15.6)') "tbar: ", this%tbar
            write(*,*) "**************************************"
        end if

    end subroutine VDSetDerived

end module mod_vd2020
