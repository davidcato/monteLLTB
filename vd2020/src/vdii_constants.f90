module vdii_constants

  integer, parameter :: stdout = 6, stderr = 0


  integer, parameter :: dl=kind(1.d0), sp=kind(1.e0)
  real(dl), parameter :: pi=3.14159265358979324_dl

  public isnan
interface isnan
    module procedure isnan_sc, isnan_vec, isnan_2d
end interface isnan

contains


    function isnan_sc(var)
        implicit none
        real(kind(1.d0)) var
        logical isnan_sc

        isnan_sc = var /= var
    end function isnan_sc

    function isnan_vec(var)
        implicit none
        real(kind(1.d0)) var(:)
        logical isnan_vec(size(var))

        isnan_vec = var /= var
    end function isnan_vec

    function isnan_2d(var)
        implicit none
        real(kind(1.d0)) var(:,:)
        logical isnan_2d(size(var(:,1)),size(var(1,:)))

        isnan_2d = var /= var
    end function isnan_2d


end module vdii_constants