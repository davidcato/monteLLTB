module wlltb_constants

  integer, parameter :: stdout = 6, stderr = 0


  integer, parameter :: dl=kind(1.d0), sp=kind(1.e0)
  real(dl), parameter :: pi=3.14159265358979324_dl
  
  type normalizations
    real(dl) :: aBB, abar, a0
  end type normalizations
  type(normalizations), parameter :: norms = normalizations(1.d-10,1.d0/3.001d3, 1.d0)
  
  type indices
    integer :: altb, adot, raprime, For3, rhow, aout
    integer :: max
  end type indices
  
  type(indices), parameter :: index = indices(1,2,3,4,5,6,6)

  integer :: lfeedback = 0
  
end module wlltb_constants