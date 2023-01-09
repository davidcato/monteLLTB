module lltb_background_eqs
use lltb_params
use lltb_complex_tools
!use lltb_prefuncs
use lltb_precision
use lltb_background_errors
implicit none

contains


	function HoverH0_squared(omokol,aa,nocheck)
		real(dl), intent(in) :: omokol(:),aa
		real(dl) :: HoverH0_squared
    logical, optional :: nocheck

		HoverH0_squared = omokol(1) / aa**3 + omokol(2) / aa**2 + omokol(3)
    
    if(present(nocheck))then
      if(nocheck)return
    end if
		if(HoverH0_squared<0._dl-accuracy(dl)*acceptance_pole)then
                   write(*,*)'For a:',aa
                   write(*,*)'H^2:',HoverH0_squared
                   write(*,*)'Distance from root, relative:',aa
                   call lltb_error('H^2 < 0.')
                else if(HoverH0_squared<0._dl)then
                   HoverH0_squared=0._dl
                end if

	end function HoverH0_squared

	function HoverH0(omokol,aa)
		real(dl), intent(in) :: omokol(:),aa
		real(dl) :: HoverH0

		HoverH0 = (HoverH0_squared(omokol,aa))**(0.5_dl)
!		write(*,*)'sq(h(a)), h(a)', dadt_oversqH0, HoverH0_squared(omokol,aa)

	end function HoverH0

	function dadt_overH0(omokol,aa)
		real(dl), intent(in) :: omokol(:),aa
		real(dl) :: dadt_overH0

		dadt_overH0 = (HoverH0(omokol,aa)) * aa
!		write(*,*)'sq(h(a)), h(a)', dadt_oversqH0, HoverH0_squared(omokol,aa)

	end function dadt_overH0

	function dH2daoverH0(omokol,aa)
		real(dl), intent(in) :: omokol(:),aa
		real(dl) :: dH2daoverH0

		dH2daoverH0 = -3._dl*omokol(1) / aa**4 - 2._dl*omokol(2) / aa**3
!		write(*,*)'sq(h(a)), h(a)', dadt_oversqH0, HoverH0_squared(omokol,aa)

	end function dH2daoverH0

	function adotdotoverH0(omokol,aa)
          real(dl), intent(in) :: omokol(:),aa
          real(dl) :: adotdotoverH0
          
          adotdotoverH0 = -0.5_dl*omokol(1) / aa**2 + omokol(3) * aa
          !		write(*,*)'sq(h(a)), h(a)', dadt_oversqH0, HoverH0_squared(omokol,aa)
          
        end function adotdotoverH0

	function dHdaoverH0(omokol,aa)
		real(dl), intent(in) :: omokol(:),aa
		real(dl) :: dHdaoverH0

		dHdaoverH0 = (dH2daoverH0(omokol,aa)) / HoverH0(omokol,aa) * 0.5_dl
		!		write(*,*)'sq(h(a)), h(a)', dadt_oversqH0, HoverH0_squared(omokol,aa)

	end function dHdaoverH0


end module lltb_background_eqs
