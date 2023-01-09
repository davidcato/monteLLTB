module wlltb_types
  use wlltb_constants
  use kofrattributes
  implicit none
  
  type test
    real(dl) :: x
  end type test
  
  ! Awesome fortran 2003 feature:
  pointer :: kofratt_point
  interface
     function kofratt_point(r,att)
        use wlltb_constants
        use kofrattributes
        implicit none
        real(dl) :: kofratt_point
        real(dl), intent(in) :: r
        type(kofr_attributes), intent(in) :: att
     end function kofratt_point
  end interface

!  pointer :: kofr_point
  interface
     function kofr_point_p(r) result(result)
        use wlltb_constants
        implicit none
        real(dl) :: result
        real(dl), intent(in) :: r
     end function kofr_point_p
  end interface
  procedure(kofr_point_p), pointer :: kofr_point
  type AllCosmology
     procedure(kofratt_point), pointer, nopass :: kofr, dkdr, tbbofr, dtbbdr
     type(kofr_attributes) :: kofr_att
  end type AllCosmology

  type funcscontainer
     procedure(kofr_point), pointer, nopass :: kofr, dkdr, tbbofr, dtbbdr
  end type funcscontainer

  type wlltbparams
    real(dl) :: dummy
  end type wlltbparams


!                      write(thisstream,'(A20,10ES31.21)')' H0_inf: ', mylbmem%H0_inf
!                      write(thisstream,'(A20,10ES31.21)')' Lambda: ', mylbmem%Lambda
!                      write(thisstream,'(A20,10ES31.21)')' Mt: ', mylbmem%Mt
!                      write(thisstream,'(A20,10ES31.21)')' r: ', mylbmem%r
!                      write(thisstream,'(A20,10ES31.21)')' t: ', mylbmem%t
!                      write(thisstream,'(A20,10ES31.21)')' k: ', mylbmem%k
!                      write(thisstream,'(A20,10ES31.21)')' t0: ', mylbmem%t0

  type lb_mymem
    real(Dl) :: H0_inf, Lambda, Mt, r, t, k, t0, w, dk, tbb, dtbb
           real(dl) :: Rltb, Rp, Rpd, S, Sd, a, H, ap, tturn, Hp, apd
           real(dl) :: Rpdd, Sdd, apdd, add, F, rhow, For3, deltarhow, deltam
  end type lb_mymem

end module wlltb_types