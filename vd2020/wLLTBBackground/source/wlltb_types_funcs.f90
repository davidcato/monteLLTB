!*
! *  wlltb_types_funcs.f90
! *  wLLTBBackground
! *
! *  Created by Wessel Valkenburg on 29/11/2011.
! *  Copyright 2011 Wessel Valkenburg. All rights reserved.
! *
! */

! 
module wlltb_types_funcs
  use wlltb_types
  implicit none

  interface operator (.eq.)
    module procedure IsEqual
  end interface 

  private
  
  public :: operator(.eq.), lbmemInitZero

contains 

  function IsEqual(lb1,lb2)
    implicit none
    type(lb_mymem), intent(in) :: lb1, lb2
    logical IsEqual
    !    real(Dl) :: H0_inf, Lambda, Mt, r, t, k, t0, w
    
    IsEqual = (lb1%H0_inf == lb2%H0_inf) .and. (lb1%Lambda == lb2%Lambda) .and. AlmostEq(lb1%Mt,lb2%Mt) &
              .and. (lb1%r == lb2%r) .and. (lb1%t == lb2%t) .and. AlmostEq(lb1%k, lb2%k) &
              .and. (lb1%w == lb2%w) & !.and. AlmostEq(lb1%t0, lb2%t0)
              .and. AlmostEq(lb1%tbb, lb2%tbb) &
              .and. AlmostEq(lb1%dtbb, lb2%dtbb) &
              .and. AlmostEq(lb1%dk, lb2%dk) 
              
  end function IsEqual


  function AlmostEq(x,y)
    implicit none
    real(dl), intent(in) :: x,y
    logical AlmostEq

    if(x==0.d0 .or. y==0.d0)then
      AlmostEq=abs(x-y)<1.d-15
    else
      AlmostEq = abs(x/y - 1.d0)<1.d-15
    end if
    
  end function AlmostEq

  subroutine lbmemInitZero(lb1)
    implicit none
    type(lb_mymem), intent(inout) :: lb1
    
    lb1= lb_mymem(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
        0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
    
  end subroutine lbmemInitZero

end module wlltb_types_funcs
