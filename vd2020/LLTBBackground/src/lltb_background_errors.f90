module lltb_background_errors
use lltb_params
use lltb_precision

contains


	subroutine lltb_set_error(vpi,lbmem,attachement)
		type(voidparams), optional :: vpi
                type(lb_mymem), optional :: lbmem
                character(len=*), optional :: attachement

		if(present(vpi))call lltb_error_routine(vpi=vpi,display=.false.)
		if(present(lbmem))then
                   call lltb_error_routine(lbmem=lbmem,display=.false.)
                end if
		if(present(attachement))call lltb_error_routine(attachement=attachement,display=.false.)

	end subroutine lltb_set_error

	subroutine lltb_error(msg)
		character(len=*), optional :: msg

		if(present(msg))then
			call lltb_error_routine(msg=msg,display=.true.)
		else
			call lltb_error_routine(display=.true.)
		end if

	end subroutine lltb_error

	subroutine lltb_error_routine(vpi,lbmem,attachement,msg,display)
		implicit none
			real, allocatable :: CrashMe(:)
		type(voidparams), optional :: vpi
		type(lb_mymem), optional :: lbmem
		character(len=*), optional :: msg
		character(len=*), optional :: attachement
		type(voidparams), save :: myvpi
		type(lb_mymem), save :: mylbmem
                character(len=10000), save :: myattachement
		logical :: display
                logical, save :: haslbmem=.false., hasattachement=.false.
                integer :: i, thisstream

		if(.not.display)then
			if(present(vpi))myvpi=vpi
			if(present(lbmem))then
                           haslbmem=.true.
                           mylbmem=lbmem
                        end if
			if(present(attachement))then
                           hasattachement=.true.
                           myattachement=trim(adjustl(attachement))
                        end if
      return
		end if
                
                do i = 1, 1 ! Skip the one to standard output
                   if(i==1)thisstream=0
                   if(i==2)thisstream=6

                   write(thisstream,*)
                   write(thisstream,*)
                   write(thisstream,'(A)')'		---------------------------------------'
                   write(thisstream,'(A)')'		---------------------------------------'
                   write(thisstream,'(A)')'		EEEEE  RRRRR   RRRRR     OOOO    RRRRR'
                   write(thisstream,'(A)')'		E      R   RR  R   RR   OO  OO   R   RR'
                   write(thisstream,'(A)')'		EEE    RRRRR   RRRRR   OO    OO  RRRRR'
                   write(thisstream,'(A)')'		E      R RR    R RR     OO  OO   R RR'
                   write(thisstream,'(A)')'		EEEEE  R  RRR  R  RRR    OOOO    R  RRR'
                   write(thisstream,'(A)')'		---------------------------------------'
                   write(thisstream,'(A)')'		---------------------------------------'
                   write(thisstream,*)
                   write(thisstream,*)


                   if(present(msg))write(thisstream,'(A20,A)')'Error message:','     '//trim(msg)

                   ! voidparams:
                   !		real(dl) :: k
                   !		real(dl) :: mtilde
                   !		real(dl) :: lambda
                   !		real(dl) :: xyz(3), workxyz(3), workxyznonzero(3) ! xyz(1) = Omega_m, xyz(2) = Omega_k, xyz(3) = Omega_Lambda
                   !		real(dl) :: workrescalefac, workxyzfac
                   !		complex(dl) :: uvw(3) ! uvw = roots of [xyz(1) + xyz(2) a + xyz(3) a^3 == 0] in a.
                   !		integer :: posrealroots ! number of positive real roots
                   !		logical :: crossedbranchcut

                   write(thisstream,'(A20,10ES31.21)')' k: ', myvpi%k
                   write(thisstream,'(A20,10ES31.21)')' mtilde: ', myvpi%mtilde
                   write(thisstream,'(A20,10ES31.21)')' lambda: ', myvpi%lambda
                   write(thisstream,'(A20,10ES31.21)')' xyz: ', myvpi%xyz
                   write(thisstream,'(A20,10ES31.21)')' workxyz: ', myvpi%workxyz
                   write(thisstream,'(A20,10ES31.21)')' workxyznonzero: ', myvpi%workxyznonzero
                   write(thisstream,'(A20,10ES31.21)')' workrescalefac: ', myvpi%workrescalefac
                   write(thisstream,'(A20,10ES31.21)')' workxyzfac: ', myvpi%workxyzfac
                   write(thisstream,'(A20,3("   (",ES13.4,",",ES13.4,")"))')'uvw: ', myvpi%uvw
                   write(thisstream,'(A20,3("   (",ES13.4,",",ES13.4,")"))')'uvw_linear: ', myvpi%uvw_linear
                   write(thisstream,'(A20,10I9)')' posrealroots: ', myvpi%posrealroots
                   write(thisstream,'(A20,10ES31.21)')' smallestrealroot: ', myvpi%smallestrealroot
                   write(thisstream,'(A20,10ES31.21)')' smallestlinroot: ', myvpi%smallestlinroot
                   write(thisstream,'(A20,10I9)')' smallestrootindex: ', myvpi%smallestrootindex
                   write(thisstream,'(A20,10I9)')' smallestlinrootindex: ', myvpi%smallestlinrootindex
                   write(thisstream,'(A20,10ES31.21)')' tturn: ', myvpi%tturn
                   write(thisstream,'(A20,3L5)')'crossedbranchcut: ', myvpi%crossedbranchcut
                   write(thisstream,'(A20,A5,A)')'ExpCase: ', myvpi%ExpCase, ' (expansion for t(a))'
                   write(thisstream,'(A20,A5,A)')'ExpCaseInv: ', myvpi%ExpCaseInv,' (expansion for inversion a(t))'

                   if(haslbmem)then
                      write(thisstream,*)
                      write(thisstream,'(A)')'		---------------------------------------'
                      write(thisstream,'(A)')'                From lltb_background:'
                      write(thisstream,'(A)')'		---------------------------------------'
                      write(thisstream,*)
                      write(thisstream,'(A20,10ES31.21)')' H0_inf: ', mylbmem%H0_inf
                      write(thisstream,'(A20,10ES31.21)')' Lambda: ', mylbmem%Lambda
                      write(thisstream,'(A20,10ES31.21)')' Mt: ', mylbmem%Mt
                      write(thisstream,'(A20,10ES31.21)')' r: ', mylbmem%r
                      write(thisstream,'(A20,10ES31.21)')' t: ', mylbmem%t
                      write(thisstream,'(A20,10ES31.21)')' k: ', mylbmem%k
                      write(thisstream,'(A20,10ES31.21)')' t0: ', mylbmem%t0
                      write(thisstream,'(A20,10ES31.21)')' Mt2: ', mylbmem%Mt2
                      write(thisstream,'(A20,10ES31.21)')' Lambdaover3: ', mylbmem%Lambdaover3
                      write(thisstream,'(A20,10ES31.21)')' H0_inf2: ', mylbmem%H0_inf2

                   end if

                   if(hasattachement)then
                      write(thisstream,*)
                      write(thisstream,'(A)')'		---------------------------------------'
                      write(thisstream,'(A)')'                Attachement:'
                      write(thisstream,'(A)')'		---------------------------------------'
                      write(thisstream,*)
                      write(thisstream,'(A)')trim(adjustl(myattachement))
                   end if

                end do

		write(0,*)'lltb_error! Provoking a crash so you can backtrace..'
		! write(*,*)CrashMe(1)
    ! call exit(-1)
		! stop 'lltb_error'
    call exit(0)

	end subroutine

end module lltb_background_errors
