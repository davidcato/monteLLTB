module vdii_vdipatch_FakeCMBParams
!  use vdii_newcambparams
  implicit none


  type voidcmbparams
     integer profile
     real zb, Hb_rat, delta0, Hlocal, TVoverTF, Theta, t0
!     real b2, Qb2, nstwo, param5, param2, param6, ToverT
     real alpha, beta
     real L2_1, d0_2
     logical :: invoid=.false., frommcmc=.true., effobs=.false., outside=.false.
     logical :: novoid=.false., void = .true.
     real MCMComk, MCMComdmh2
     ! camb new:
     logical :: splitISW, withISW
  end type voidcmbparams
! end void new

  Type FakeCMBParams
     real nuisance(1:5)
      !unit Gaussians for experimental parameters
     real norm(1:5)
      !These are fast parameters controling amplitudes, calibrations, etc.
     real InitPower(1:5)
      !These are fast paramters for the initial power spectrum
     !Now remaining (non-independent) parameters
     real omb, omc, omv, omnu, omk, omdm
     real ombh2, omch2, omnuh2, omdmh2
     real omegabh2, omegadmh2
     real zre, nufrac
     real h, H0
     real w
     real YHe, nnu, zreio
     real reserved(5)
! void new
     type(VoidcmbParams) :: vcp
! end void new
  end Type FakeCMBParams




end module vdii_vdipatch_FakeCMBParams
