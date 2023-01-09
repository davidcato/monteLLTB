module lltb_background_ht_ddok_funcs
use lltb_background_errors
use lltb_params
use lltb_complex_tools
use lltb_prefuncs
use lltb_precision
implicit none


contains

! Deprecated: see lltb_background_limits now.
!!$	subroutine ddok_tofamax(vpi, tc, t)
!!$		use elliptics
!!$		implicit none
!!$		type(voidparams) :: vpi
!!$!		real(dl), intent(in) :: xyz(3)
!!$!		complex(dl), intent(in) :: uvw(3)
!!$!		real(dl), intent(in) :: a
!!$		real(dl), intent(out), optional :: t
!!$		complex(dl), intent(out), optional :: tc
!!$
!!$
!!$		complex(dl) :: umv, umw, vmu, vmw, wmu, wmv, uvw(3)
!!$		integer :: i, thisroot, j, indices(3)
!!$!		real(dl) :: am1
!!$		complex(dl) :: sqzp3, rd_args(3), rd_res(2), squvwm1, f1, tcmplx, um1c
!!$
!!$		sqzp3 = (dcmplx(vpi%workxyz(3),0._dl))**(1.5_dl)
!!$
!!$		um1c=cmplx(1._dl,0._dl)/vpi%smallestlinroot
!!$		squvwm1 = dcmplx(1._dl,0._dl) / (vpi%uvw_linear(1)*vpi%uvw_linear(2)*vpi%uvw_linear(3))**(0.5_dl)
!!$
!!$
!!$                ! which of the three is the actual smallest positive real root?
!!$                ! vpi%smallestlinrootindex
!!$                thisroot=vpi%smallestlinrootindex
!!$                if(thisroot<1)thisroot=vpi%smallestrootindex
!!$
!!$                indices(1)=thisroot
!!$
!!$                rd_args(1)=cmplx(0._dl,0._dl)
!!$               j=0
!!$		do i = 1, 3
!!$                   if(i==thisroot)cycle
!!$                   j=j+1
!!$                   indices(j+1)=i
!!$                   rd_args(j+1) = um1c - 1._dl / vpi%uvw_linear(i)
!!$                   if(j==2)exit
!!$		end do
!!$
!!$		uvw=vpi%uvw_linear
!!$		umv=1._dl/(uvw(indices(1))-uvw(indices(2)))
!!$		umw=1._dl/(uvw(indices(1))-uvw(indices(3)))
!!$		vmu=1._dl/(uvw(indices(2))-uvw(indices(1)))
!!$		vmw=1._dl/(uvw(indices(2))-uvw(indices(3)))
!!$		wmu=1._dl/(uvw(indices(3))-uvw(indices(1)))
!!$		wmv=1._dl/(uvw(indices(3))-uvw(indices(2)))
!!$
!!$
!!$		f1 = dcmplx(0.,1._dl)/3._dl/sqzp3 * squvwm1
!!$
!!$
!!$		rd_res(1)=f1*(vmu*vmw-umv*umw)*rd(rd_args(3),rd_args(1),rd_args(2))
!!$		rd_res(2)=f1*(wmu*wmv-umv*umw)*rd(rd_args(1),rd_args(2),rd_args(3))
!!$
!!$		tcmplx = (rd_res(1)+rd_res(2))
!!$
!!$
!!$		if(vpi%crossedbranchcut)then
!!$!			write(*,*)'Crossed branchcut, taking account in ddok_tofa_mkl.'
!!$!			call lltb_error('Stopping, this scenario has not been tested yet. Test it!')
!!$			tcmplx=-tcmplx
!!$		end if
!!$
!!$		if(present(tc))tc = tcmplx
!!$		if(present(t))t=killimag(tcmplx)
!!$
!!$
!!$        end subroutine ddok_tofamax
!!$

	subroutine ddok_tofa_mkl(vpi,a, tc, t)
		use elliptics
		implicit none
		type(voidparams) :: vpi
!		real(dl), intent(in) :: xyz(3)
!		complex(dl), intent(in) :: uvw(3)
		real(dl), intent(in) :: a
		real(dl), intent(out), optional :: t
		complex(dl), intent(out), optional :: tc
		complex(dl) :: umv, umw, vmu, vmw, wmu, wmv, uvw(3)
		integer :: i
		real(dl) :: am1
		complex(dl) :: sqzp3, rd_args(3), rd_res(3), squvwm1, f1, tcmplx, am1c, myuvw

		sqzp3 = (dcmplx(vpi%workxyz(3),0._dl))**(1.5_dl)
		am1 = 1._dl / a
		am1c=dcmplx(am1,0._dl)
                
                myuvw = product(vpi%uvw_linear)
                myuvw=killimag(myuvw) ! uvw is real by definition, small residue may fuck it up.

!		squvwm1 = dcmplx(1._dl,0._dl) / (vpi%uvw_linear(1)*vpi%uvw_linear(2)*vpi%uvw_linear(3))**(0.5_dl)
		squvwm1 = dcmplx(1._dl,0._dl) / (myuvw)**(0.5_dl)



		do i = 1, 3
			rd_args(i) = am1c - 1._dl / vpi%uvw_linear(i)
!			rj_args(i) = (- a + vpi%uvw(i))/(a*vpi%uvw(i))
		end do

		uvw=vpi%uvw_linear
		umv=1._dl/(uvw(1)-uvw(2))
		umw=1._dl/(uvw(1)-uvw(3))
		vmu=1._dl/(uvw(2)-uvw(1))
		vmw=1._dl/(uvw(2)-uvw(3))
		wmu=1._dl/(uvw(3)-uvw(1))
		wmv=1._dl/(uvw(3)-uvw(2))


		f1 = dcmplx(0.,-1._dl)/3._dl/sqzp3 * squvwm1

		rd_res(1)=f1*umv*umw*rd(rd_args(2),rd_args(3),rd_args(1))
		rd_res(2)=f1*vmu*vmw*rd(rd_args(3),rd_args(1),rd_args(2))
		rd_res(3)=f1*wmu*wmv*rd(rd_args(1),rd_args(2),rd_args(3))

		tcmplx = - (rd_res(1)+rd_res(2)+rd_res(3))

		if(vpi%crossedbranchcut)then
!			write(*,*)'Crossed branchcut, taking account in ddok_tofa_mkl.'
!			call lltb_error('Stopping, this scenario has not been tested yet. Test it!')
			tcmplx=-tcmplx
		end if

		if(present(tc))tc = tcmplx
		if(present(t))t=killimag(tcmplx)

	end subroutine ddok_tofa_mkl

	subroutine ddok_tofa_mk(vpi,a, tc, t, adotterm)
		use elliptics
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(in) :: a
		real(dl), intent(out), optional :: t
		complex(dl), intent(out), optional :: tc
		real(dl), intent(out), optional :: adotterm
		integer :: i
		complex(dl) :: om, ol, ok
		complex(dl) :: tcmplx, tctest

		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)

                if(present(adotterm))adotterm = -0.5_Dl*(3_dl*om**1.5_dl*Log(Sqrt(ok)))/ok**2.5_dl


		tcmplx = -0.5_Dl* (   ((3_Dl*om**2.5_Dl*Log(Sqrt(ok)*Sqrt(om)))/ok**2.5_Dl + &
         (om*Sqrt((A*om)/(A*ok + om))*  &
            (Sqrt(A)*Sqrt(ok)*(A*ok + 3_Dl*om) - &
              3_Dl*om*Sqrt(A*ok + om)*Log(Sqrt(A)*ok + Sqrt(ok)*Sqrt(A*ok + om))))/ &
          (Sqrt(A)*ok**2.5_Dl))/om**1.5_Dl  )

!		if(vpi%crossedbranchcut)then
!			write(*,*)'Crossed branchcut, taking account in ddok_tofa_mk.'
!			call lltb_error('Stopping, this scenario has not been tested yet. Test it!')
! !! Do not do this	tcmplx=-tcmplx
!		end if

		if(present(tc))tc = tcmplx
		if(present(t))t=killimag(tcmplx)
	end subroutine ddok_tofa_mk


	subroutine ddok_tofa_m(vpi,a, tc, t)
		use elliptics
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(in) :: a
		real(dl), intent(out), optional :: t
		complex(dl), intent(out), optional :: tc
		integer :: i
		complex(dl) :: om, ol, ok
		complex(dl) :: tcmplx

		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)

		tcmplx = - 0.2_dl * om**(-1.5_dl) * a**(2.5_dl)

!		if(vpi%crossedbranchcut)then
!			write(*,*)'Crossed branchcut, taking account in ddok_tofa_m.'
!			call lltb_error('Stopping, this scenario has not been tested yet. Test it!')
! !! Do not do this	tcmplx=-tcmplx
!		end if

		if(present(tc))tc = tcmplx
		if(present(t))t=killimag(tcmplx)



	end subroutine ddok_tofa_m


	subroutine ddok_pert_mk_l(vpi,a, tc, t)
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(in) :: a
		real(dl), intent(inout), optional :: t
		real(dl) :: ayx
		complex(dl) :: om, ol, ok
		complex(dl), intent(inout), optional :: tc
		complex(dl) :: tcmplx, sqayx, xy32, tcmplx_2nd

		ayx = a * vpi%workxyz(2) / vpi%workxyz(1)
		sqayx = sqrt(dcmplx(ayx,0.))
		xy32 = dcmplx(vpi%workxyz(1) / vpi%workxyz(2),0.)**(1.5_dl)

		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)


		tcmplx =  ol*   (3_dl*((105_dl*om**3*Log(Sqrt(ok)*Sqrt(om)))/(8._dl*ok**5.5_dl) + &
           (Sqrt(A)*Sqrt(ok)*(8_dl*A**4*ok**4 - 18_dl*A**3*ok**3*om + 63_dl*A**2*ok**2*om**2 + &
                 420_dl*A*ok*om**3 + 315_dl*om**4) - &
              315_dl*om**3*(A*ok + om)**1.5_dl*Log(Sqrt(A)*ok + Sqrt(ok)*Sqrt(A*ok + om)))/ &
            (24._dl*ok**5.5_dl*(A*ok + om)**1.5_dl)))/4._dl

		tcmplx_2nd =    ol**2 *     (-15*((9009*om**5*Log(Sqrt(ok)*Sqrt(om)))/(128.*ok**8.5) + &
           (Sqrt(A)*Sqrt(ok)*(128*A**7*ok**7 - 240*A**6*ok**6*om + &
                 520*A**5*ok**5*om**2 - 1430*A**4*ok**4*om**3 + &
                 6435*A**3*ok**3*om**4 + 69069*A**2*ok**2*om**5 + 105105*A*ok*om**6 + &
                 45045*om**7) - 45045*om**5*(A*ok + om)**2.5* &
               Log(Sqrt(A)*ok + Sqrt(ok)*Sqrt(A*ok + om)))/ &
            (640.*ok**8.5*(A*ok + om)**2.5)))/16.


		if(present(tc))then
!      if(abs(tcmplx_2nd)>abs(tcmplx))write(*,*)'This may be a thing here.'
!			if (abs(tcmplx_2nd)/abs(tc)>accuracy(dl) * rel_ACCEPTANCE(3) ) then
			if (abs(tcmplx_2nd)/abs(tc)>accuracy(dl) * rel_ACCEPTANCE(3) .or. abs(tcmplx)>abs(tc)) then
			! check for blowing up expansion: in that case the full elliptic integral is more accurate, so send a message t = -1
				if(lfeedback>1)then
					writE(*,*)'2nd order in ddok_pert_mk_l', tc, tc+tcmplx_2nd
					write(*,*)vpi%workxyz(1)/a/vpi%workxyz(2),">0 or <-1"
				end if
				tc=-1._dl
      else
        tc = tc + tcmplx
			end if
		end if
		if(present(t).and.present(tc))then
			t=killimag(tc)
		else if (present(t)) then
			if (abs(tcmplx)>abs(t)) then
        t=-1._dl
      else
        t=t+killimag(tcmplx)
      end if
		end if


	end subroutine ddok_pert_mk_l

	subroutine ddok_pert_m_k(vpi,a, tc, t)
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(in) :: a
		real(dl), intent(inout), optional :: t
		complex(dl) :: ok,om, ol
		complex(dl), intent(inout), optional :: tc
		complex(dl) :: tcmplx!, sqayx, xy32
! For testing:
		complex(dl) :: tctest

		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)
!
!
		tcmplx = ok * (3._dl*A**3.5_dl)/(14._dl*om**2.5_dl)


		if(present(tc))tc = tc + tcmplx
		if(present(t).and.present(tc))then
			t=killimag(tc)
		else if (present(t)) then
			t=t+killimag(tcmplx)
		end if


	end subroutine ddok_pert_m_k

	subroutine ddok_pert_m_l(vpi,a, tc, t)
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(in) :: a
		real(dl), intent(inout), optional :: t
		complex(dl) :: om, ol, ok
		complex(dl), intent(inout), optional :: tc
		complex(dl) :: tcmplx!, sqayx, xy32


		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)


		tcmplx = ol * (3_dl*A**5.5_dl)/(22._dl*om**2.5_dl)


		if(present(tc))tc = tc + tcmplx
		if(present(t).and.present(tc))then
			t=killimag(tc)
		else if (present(t)) then
			t=t+killimag(tcmplx)
		end if


	end subroutine ddok_pert_m_l




end module lltb_background_ht_ddok_funcs
