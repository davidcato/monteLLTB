module lltb_background_ht_funcs
use lltb_background_errors
use lltb_params
use lltb_complex_tools
use lltb_prefuncs
use lltb_precision
use lltb_background_ht_ddok_funcs
implicit none


contains



	subroutine tofa_mkl(vpi,a, tc, t)
		use elliptics
		implicit none
		type(voidparams) :: vpi
!		real(dl), intent(in) :: xyz(3)
!		complex(dl), intent(in) :: uvw(3)
		real(dl), intent(in) :: a
		real(dl), intent(out), optional :: t
		complex(dl), intent(out), optional :: tc

		complex(dl) :: myuvw

		integer :: i
		real(dl) :: am1
		complex(dl) :: sqz, rj_args(3), rj_res, squvwm1, f1, tcmplx, am1c

		sqz = sqrt(dcmplx(vpi%workxyz(3),0.))
		am1 = 1._dl / a
		am1c=dcmplx(am1,0._dl)

!		squvwm1 = dcmplx(1._dl,0._dl) / sqrt(vpi%uvw_linear(1)*vpi%uvw_linear(2)*vpi%uvw_linear(3))
                myuvw=vpi%uvw_linear(1)*vpi%uvw_linear(2)*vpi%uvw_linear(3)
                ! the product uvw is real by definition. Negative imaginary residue leads to wrong sign of sqrt.
                myuvw=killimag(myuvw)
                squvwm1=dcmplx(1._dl,0._dl)/sqrt(myuvw)
		rj_args = am1c - dcmplx(1._dl,0._dl) / vpi%uvw_linear
!write(*,'(3(" (",ES17.10,",",ES17.10,") "))')rj_args
		rj_res = rj(rj_args(1),rj_args(2),rj_args(3),am1c)



		f1 = dcmplx(0.,2)/dcmplx(3._dl,0._dl)/sqz * squvwm1
		tcmplx = f1*rj_res


		if(mreal(tcmplx)<0.)then
			write(*,*)'BRANCH CUT!'
      write(*,*)'for parameters xyz:',vpi%xyz
      write(*,*)'for a:',a
      write(*,*)'for realroot:',vpi%smallestrealroot
if(a<1.d200)then
  lfeedback=-1
else
  tcmplx=-2.d0
end if
if(lfeedback<0)write(*,'("(",ES12.5,",",ES13.5,")")')f1,sqrt(vpi%uvw_linear(1))*sqrt(vpi%uvw_linear(2))*sqrt(vpi%uvw_linear(3))/sqrt(vpi%uvw_linear(1)*vpi%uvw_linear(2)*vpi%uvw_linear(3)),sqz,squvwm1,sqrt(vpi%uvw_linear(1))*sqrt(vpi%uvw_linear(2))*sqrt(vpi%uvw_linear(3)),sqrt(vpi%uvw_linear(1)*vpi%uvw_linear(2)*vpi%uvw_linear(3)),vpi%uvw_linear,product(vpi%uvw_linear)
			vpi%crossedbranchcut = .true.
			tcmplx=-tcmplx ! Sometimes we end up with the wrong square root..
if(lfeedback<0)write(*,*)vpi%crossedbranchcut
if(lfeedback<0)then
   if(vpi%crossedbranchcut)stop
end if
		end if

		if(present(tc))tc = tcmplx
		if(present(t))t=killimag(tcmplx)


	end subroutine tofa_mkl


	subroutine tofa_kl(vpi,a, tc, t)
		use elliptics
		implicit none
		type(voidparams) :: vpi
!		real(dl), intent(in) :: xyz(3)
!		complex(dl), intent(in) :: uvw(3)
		real(dl), intent(in) :: a
		real(dl), intent(out), optional :: t
		complex(dl), intent(out), optional :: tc
		complex(dl) :: as_arg, tcmplx

		as_arg = dcmplx(a,0.)*sqrt(dcmplx(vpi%workxyz(3)/vpi%workxyz(2),0.))

		tcmplx = masinh(as_arg) / sqrt(dcmplx(vpi%workxyz(3),0.))

		if(vpi%workxyz(2)<0._dl)tcmplx=tcmplx*dcmplx(0,-1.)

		if(mreal(tcmplx)<0.)tcmplx=-tcmplx ! Sometimes we end up with the wrong square root.. Branch cut.

		if(present(tc))tc = tcmplx
		if(present(t))t=killimag(tcmplx)

	end subroutine tofa_kl


	subroutine tofa_km(vpi,a, tc, t)
		use elliptics
		implicit none
		type(voidparams) :: vpi
!		real(dl), intent(in) :: xyz(3)
!		complex(dl), intent(in) :: uvw(3)
		real(dl), intent(in) :: a
		real(dl), intent(out), optional :: t
		real(dl) :: ayx
		complex(dl), intent(out), optional :: tc
		complex(dl) :: tcmplx, sqayx, xy32

		ayx = a * vpi%workxyz(2) / vpi%workxyz(1)
		sqayx = sqrt(dcmplx(ayx,0.))
                ! Goes wrong for negative k:
		! xy32 = dcmplx(vpi%workxyz(1) / vpi%workxyz(2),0.)**(1.5_dl)
                ! Goes right for negative k: (see appendix b)
                xy32 = dcmplx(1._dl/vpi%workxyz(1) * vpi%workxyz(2),0.)**(-1.5_dl)

		tcmplx = (sqayx * sqrt(dcmplx(1._dl+ayx,0.)) - masinh(sqayx)) * xy32 / sqrt(dcmplx(vpi%workxyz(1),0.))

!		if(vpi%workxyz(2)<0._dl)tcmplx=tcmplx*dcmplx(0,-1.)
		if(mreal(tcmplx)<0.)then
                   !write(*,*)'km branch cut.'
                   vpi%crossedbranchcut=.true.
                   tcmplx=-tcmplx ! Sometimes we end up with the wrong square root.. Branch cut.
                end if

		if(present(tc))tc = tcmplx
		if(present(t))t=killimag(tcmplx)


	end subroutine tofa_km



	subroutine tofa_pert_ml_k(vpi,a, tc, t)
		implicit none
		type(voidparams) :: vpi!, myvpi
		real(dl), intent(in) :: a
		real(dl), intent(inout), optional :: t
		real(dl) :: ayx, mya, myt
		complex(dl) :: om, ol, ok
		complex(dl), intent(inout), optional :: tc
		complex(dl) :: tcmplx, sqayx, xy32, tcmplx_2nd, mytc

		if(present(tc))mytc=tc
		if(present(t))myt=t
		mya=a

!		myvpi=vpi
!		myvpi%xyz(2)=0._dl
!		call void_makeworkxyz(myvpi,a)
!		call void_uvw(myvpi%workxyz,myvpi%uvw, myvpi%posrealroots, myvpi%smallestrealroot)
		call ddok_tofa_mkl(vpi,mya, tc=tcmplx, t=myt)
!


		tcmplx= tcmplx*vpi%xyz(2)

		if(present(tc))then
			tc = tc + tcmplx
		end if
		if(present(t).and.present(tc))then
			t=killimag(tc)
		else if (present(t)) then
			t=t+killimag(tcmplx)
		end if

	end subroutine tofa_pert_ml_k

	subroutine tofa_pert_mk_l(vpi,a, tc, t)
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(in) :: a
		real(dl), intent(inout), optional :: t
		real(dl) :: ayx, estd_2nd(2)
		complex(dl) :: om, ol, ok
		complex(dl), intent(inout), optional :: tc
		complex(dl) :: tcmplx, sqayx, xy32, tcmplx_2nd
		ayx = a * vpi%workxyz(2) / vpi%workxyz(1)
		sqayx = sqrt(dcmplx(ayx,0.))
		xy32 = dcmplx(vpi%workxyz(1) / vpi%workxyz(2),0.)**(1.5_dl)

		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
		ol=dcmplx(vpi%workxyz(3),0._dl)


		tcmplx =     -(ol*((35_dl*Sqrt(1_dl/om)*om**3.5_dl*Log(Sqrt(ok)*Sqrt(om)))/(8._dl*ok**4.5_dl) + &
            (Sqrt(A/(A*ok + om))*(Sqrt(A)*Sqrt(ok)*(8_dl*A**3*ok**3 - 14_dl*A**2*ok**2*om + 35_dl*A*ok*om**2 + 105_dl*om**3) - &
                 105_dl*om**3*Sqrt(A*ok + om)*Log(Sqrt(A)*ok + Sqrt(ok)*Sqrt(A*ok + om))))/(24._dl*Sqrt(A)*ok**4.5_dl)))/2._dl

		tcmplx_2nd =    (3_dl*ol**2*((3003_dl*Sqrt(1_dl/om)*om**5.5_dl*Log(Sqrt(ok)*Sqrt(om)))/(128._dl*ok**7.5_dl) + &
           (Sqrt(A/(A*ok + om))*(Sqrt(A)*Sqrt(ok)*Sqrt(A*ok + om)*(384_dl*A**6*ok**6 - 624_dl*A**5*ok**5*om + 1144_dl*A**4*ok**4*om**2 - 2574_dl*A**3*ok**3*om**3 + &
                   9009_dl*A**2*ok**2*om**4 + 60060_dl*A*ok*om**5 + 45045_dl*om**6) - 45045_dl*om**5*(A*ok + om)**2*Log(Sqrt(A)*ok + Sqrt(ok)*Sqrt(A*ok + om))))/ &
            (1920._dl*Sqrt(A)*ok**7.5_dl*(A*ok + om)**1.5_dl)))/8._dl

                ! This expansion seems unstable. Perform crosscheck, estimate error:
                estd_2nd(1)=(abs(tcmplx)/abs(tc))**2
                estd_2nd(2)=abs(tcmplx_2nd)/abs(tc)

                if(abs(maimag(tcmplx))/abs(mreal(tcmplx))>accuracy(dl)*acceptance)then
                   estd_2nd=1.d3
                end if

		if(present(tc))then
			tc = tc + tcmplx
			if (any(estd_2nd>accuracy(dl) * rel_ACCEPTANCE(3)) ) then
!			if (abs(tcmplx_2nd)/abs(tc)>accuracy(dl) * rel_ACCEPTANCE(3) ) then
!			if (abs(tcmplx_2nd)/abs(tc)>accuracy(dl) * acceptance_pure ) then
			! check for blowing up expansion: in that case the full elliptic integral is more accurate, so send a message t = -1
				if(lfeedback>1)then
					writE(*,*)'2nd order in tofa_pert_mk_l', tc, tc+tcmplx_2nd
					write(*,*)vpi%workxyz(1)/a/vpi%workxyz(2),">0 or <-1"
				end if
				tc=-1._dl
			end if
		end if
		if(present(t).and.present(tc))then
			t=killimag(tc)
		else if (present(t)) then
			t=t+killimag(tcmplx)
		end if


	end subroutine tofa_pert_mk_l


	subroutine tofa_pert_m_k(vpi,a, tc, t)
		implicit none
		type(voidparams) :: vpi
		real(dl), intent(in) :: a
		real(dl), intent(inout), optional :: t
		complex(dl) :: ok!,om, ol
		complex(dl), intent(inout), optional :: tc
		complex(dl) :: tcmplx!, sqayx, xy32
! For testing:
		complex(dl) :: tctest

!		om=dcmplx(vpi%workxyz(1),0._dl)
		ok=dcmplx(vpi%workxyz(2),0._dl)
!		ol=dcmplx(vpi%workxyz(3),0._dl)
!
!
!		tcmplx = - 0.2_dl * ok * om**(-1.5_dl) * a**(2.5_dl)
		call ddok_tofa_m(vpi,a, tcmplx)

		! testing:
!		call ddok_tofa_mK(vpi,a, tctest)
!		write(*,*)'Compare mk: Known, new:',tcmplx,tctest
!		call ddok_tofa_mKl(vpi,a, tctest)
!		write(*,*)'Compare mkl: Known, new:',tcmplx,tctest
		tcmplx=tcmplx*ok


		if(present(tc))tc = tc + tcmplx
		if(present(t).and.present(tc))then
			t=killimag(tc)
		else if (present(t)) then
			t=t+killimag(tcmplx)
		end if


	end subroutine tofa_pert_m_k

	subroutine tofa_pert_m_l(vpi,a, tc, t)
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


		tcmplx = - 1._dl/9._dl * ol * om**(-1.5_dl) * a**(4.5_dl)


		if(present(tc))tc = tc + tcmplx
		if(present(t).and.present(tc))then
			t=killimag(tc)
		else if (present(t)) then
			t=t+killimag(tcmplx)
		end if


	end subroutine tofa_pert_m_l




end module lltb_background_ht_funcs
