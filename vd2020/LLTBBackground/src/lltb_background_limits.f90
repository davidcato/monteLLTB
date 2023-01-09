module lltb_background_limits
use elliptics
use lltb_precision
use lltb_params
use lltb_prefuncs
!!$!use lltb_prefuncs
!!$!use lltb_background_ht
!!$use lltb_background_ht_ddok
!!$use lltb_background_ht_ddok_funcs
!!$!use lltb_background_ht_funcs
!!$use lltb_background_aofht
!!$!use lltb_background_aofht_funcs
!!$use lltb_background_errors
use lltb_background_eqs
!use lltb_complex_tools

contains

  function aprime_amax(vpi,a, yprime, dtbbdr, hsign, adotterm)
    use elliptics
    implicit none
    type(voidparams) :: vpi, myvpi
    real(dl), intent(in) :: a, yprime, dtbbdr, hsign
    real(dl), optional, intent(out) :: adotterm
    real(dl) :: aprime_amax
!    complex(dl) :: thislinuvw(3)
    

    if(any(vpi%expcase==(/'LLL','LLS','LSL','SLL','LSS','SLS','SSL'/)))then

       myvpi=vpi
       call CopyUVWtoLinear(myvpi)
       aprime_amax = aprime_amax_mkl(myvpi,a, yprime, dtbbdr, hsign, adotterm)

    else if (any(vpi%expcase==(/'LLZ','LSZ','SLZ'/)))then
       aprime_amax = aprime_amax_mk(vpi,a, yprime, dtbbdr, hsign, adotterm)
    else if (dabs(yprime)<1.d-30) then
       aprime_amax=0.d0
       adotterm=0.d0
    else
       call lltb_error('Unknown case in aprime_amax')
       stop 'case not known in aprime_amax'
    end if


  end function aprime_amax

  
  function aprime_amax_mkl(vpi,a, yprime, dtbbdr, hsign, adotterm)
    use elliptics
    implicit none
    type(voidparams) :: vpi
    real(dl), intent(in) :: a, yprime, dtbbdr, hsign
    real(dl), optional, intent(out) :: adotterm
    real(dl) :: aprime_amax_mkl
    complex(dl) :: umv, umw, vmu, vmw, wmu, wmv, uvw(3)
    integer :: i, thisroot, j, indices(3)
    real(dl) :: am1
    real(dl) :: thisadot, mysign
    complex(dl) :: sqzp3, sqz,rd_args(3), rd_res(3), squvwm1, f1, tcmplx, am1c, myuvw!, um1c

    thisadot = hsign * a * HoverH0(vpi%xyz,a)
    sqz = (dcmplx(vpi%workxyz(3),0._dl))**(.5_dl)
    sqzp3 = (dcmplx(vpi%workxyz(3),0._dl))**(1.5_dl)
    am1 = 1._dl / a
    am1c=dcmplx(am1,0._dl)
!    squvwm1 = dcmplx(1._dl,0._dl) / (vpi%uvw_linear(1)*vpi%uvw_linear(2)*vpi%uvw_linear(3))**(0.5_dl)

    myuvw=product(vpi%uvw)
    myuvw=killimag(myuvw) ! prod uvw = real by definition, small imaginary residue fucks up sqrt.
    
!    squvwm1 = dcmplx(1._dl,0._dl) / (vpi%uvw(1)*vpi%uvw(2)*vpi%uvw(3))**(0.5_dl)
    squvwm1 = dcmplx(1._dl,0._dl) / (myuvw)**(0.5_dl)

    thisroot=vpi%smallestlinrootindex
    if(thisroot<1)thisroot=vpi%smallestrootindex

    indices(1)=thisroot

    rd_args(1)=am1c - cmplx(1._dl,0._dl)/vpi%smallestlinroot
    j=0
    do i = 1, 3
       if(i==thisroot)cycle
       j=j+1
       indices(j+1)=i
       rd_args(j+1) = am1c - 1._dl / vpi%uvw_linear(i)
       if(j==2)exit
    end do


    uvw=vpi%uvw_linear

    umv=1._dl/(uvw(indices(1))-uvw(indices(2)))
    umw=1._dl/(uvw(indices(1))-uvw(indices(3)))
    vmu=1._dl/(uvw(indices(2))-uvw(indices(1)))
    vmw=1._dl/(uvw(indices(2))-uvw(indices(3)))
    wmu=1._dl/(uvw(indices(3))-uvw(indices(1)))
    wmv=1._dl/(uvw(indices(3))-uvw(indices(2)))



    f1 = dcmplx(0.,1._dl)/3._dl/sqzp3 * squvwm1


    rd_res(2)=f1*(vmu*vmw-umv*umw)*rd(rd_args(3),rd_args(1),rd_args(2))
    rd_res(3)=f1*(wmu*wmv-umv*umw)*rd(rd_args(1),rd_args(2),rd_args(3))


    tcmplx = (rd_res(2)+rd_res(3))
!!$writE(*,*)'in aprimeamax:',tcmplx
!!$write(*,*)'tcmplx=',tcmplx
!!$write(*,*)'f1:',f1
!!$write(*,*)'(vmu*vmw-umv*umw)',(vmu*vmw-umv*umw)
!!$write(*,*)'(wmu*wmv-umv*umw)',(wmu*wmv-umv*umw)
!!$write(*,*)'rd_res a:',rd_res(2)
!!$write(*,*)'rd_res b:',rd_res(3)
!!$write(*,*)'thisroot:',thisroot
!!$write(*,*)'smallestlinrootindex:',vpi%smallestlinrootindex
!!$write(*,*)'smallestrootindex:',vpi%smallestrootindex
!!$write(*,*)'indices:',indices
!!$write(*,*)'u:',uvw(1)
!!$write(*,*)'v:',uvw(2)
!!$write(*,*)'w:',uvw(3)
!!$write(*,*)'rd_args(1):',rd_args(1)
!!$write(*,*)'rd_args(2):',rd_args(2)
!!$write(*,*)'rd_args(3):',rd_args(3)

    if(vpi%crossedbranchcut)then
       mysign=-1._dl
    else
       mysign=1._dl
    end if


    tcmplx=( mysign*yprime * tcmplx + dtbbdr)
    if(present(adotterm))adotterm=-mreal(tcmplx)
    aprime_amax_mkl=-thisadot*mreal(tcmplx)-mysign*yprime*mreal(dcmplx(1._dl,0._dl)/vpi%xyz(3) *umv*umw*a)

    if(isnan(real(tcmplx)).or.isnan(aimag(tcmplx)))then
       write(*,*)rd_res(1)
       write(*,*)rd_res(2)
       write(*,*)f1
       write(*,*)vmu*vmw-umv*umw
       write(*,*)wmu*wmv-umv*umw
       write(*,*)rd(rd_args(3),rd_args(1),rd_args(2))
       write(*,*)rd(rd_args(1),rd_args(2),rd_args(3))

       call lltb_error('NaN in aprime_amax_mkl')
       
    end if
  end function aprime_amax_mkl

  function aprime_amax_mk(vpi,a, yprime, dtbbdr, hsign, adotterm)
!    use  lltb_background_ht_ddok_funcs
    implicit none
    type(voidparams) :: vpi
    real(dl), intent(in) :: a, yprime, dtbbdr, hsign
    real(dl), optional, intent(out) :: adotterm
    real(dl) :: aprime_amax_mk

    complex(dl) :: thisadotterm, part2 !thisddok, 
    complex(dl) :: om,ok,ol
    real(dl) :: thisadot

!    call ddok_tofa_mk(vpi, a=a, t=thisddok,adotterm=thisadotterm)
    thisadot = hsign * a * HoverH0(vpi%xyz,a)
!    thisadot = a * HoverH0(vpi%xyz,a)

    om=cmplx(vpi%workxyz(1),0._dl)
    ok=cmplx(vpi%workxyz(2),0._dl)
    ol=cmplx(vpi%workxyz(3),0._dl)

    thisadotterm = 0.5_Dl* ((3_dl*om*Log(Sqrt(ok)*Sqrt(om)))/ok**2.5_dl - &
       (3_dl*om*Log(Sqrt(A)*ok + Sqrt(ok)*Sqrt(A*ok + om)))/ok**2.5_dl )

!-0.5_Dl*(3_dl*Sqrt(1_dl/om)*om**1.5_dl*Log(Sqrt(ok)*Sqrt(om)))/ok**2.5_dl
!                            !(3_dl*om**1.5_dl*Log(Sqrt(ok)))/ok**2.5_dl

    adotterm = mreal(thisadotterm) * yprime - dtbbdr



    part2=(A*ok + 3_dl*om)/ok**2_dl


!    part2 = (Sqrt(A)*Sqrt(ok)*(A*ok + 3*om) - &
!         3*om*Sqrt(A*ok + om)*Log(Sqrt(A)*ok + Sqrt(ok)*Sqrt(A*ok + om)))/ &
!       (Sqrt(A)*ok**2.5) 

    aprime_amax_mk=0.5_dl*part2*yprime + thisadot * adotterm


  end function aprime_amax_mk


  function Hprime_amax(vpi,a,yprime,hsign,adotterm,aprime)
    implicit none
    real(dl), intent(in) :: a,yprime,hsign,adotterm, aprime
    type(voidparams) :: vpi
    real(dl) :: Hprime_amax
    if(any(vpi%expcase==(/'LLL','LLS','LSL','SLL','LSS','SLS','SSL'/)))then
       Hprime_amax = Hprime_amax_mkl(vpi,a, yprime, hsign, adotterm)
    else if (any(vpi%expcase==(/'LLZ','LSZ','SLZ'/)))then
       Hprime_amax = Hprime_amax_mk(vpi,a, yprime, hsign, adotterm,aprime)
    else if (dabs(yprime)<1.d-30) then
       Hprime_amax=0.d0
    else
       write(*,*)yprime,'yprime'
       call lltb_error('Unknown case in Hprime_amax')
       stop 'case not known in Hprime_amax'
    end if

  end function Hprime_amax

  function Hprime_amax_mkl(vpi,a,yprime,hsign,adotterm)
    implicit none
    real(dl), intent(in) :: a,yprime,hsign,adotterm!,aprime
    type(voidparams) :: vpi
    real(dl) :: Hprime_amax_mkl
    real(dl) :: thishp!, thishp_old
    real(dl) :: thisdhdalike
    real(dl) ::  x, y, u, z!, thish
    complex(dl) :: theroot, xyu ! z can be < 0
    integer :: i, thisroot, j, indices(3)
!    real(dl) :: term1, term2, term2a

    ! Next line only valid for adot = 0:
    thisdHdalike = -(3._dl * vpi%xyz(1) / a**3 + 2._dl * vpi%xyz(2)/a**2) * 0.5_dl
    ! Normally always valid: 
!    thisadotdot_overa =  - 0.5_dl * vpi%xyz(1) / a**3 + vpi%xyz(3)
!    thish = hsign * HoverH0(vpi%xyz,a)

    x=vpi%xyz(1)
    y=vpi%xyz(2)
    z=vpi%xyz(3)
    u=vpi%smallestlinroot
    xyu = 1.5_dl * x / sqrt(dcmplx(z*a**3,0.d0)) / u / (3._dl*x/u + 2._dl * y) 

    thisroot=vpi%smallestlinrootindex
    if(thisroot<1)thisroot=vpi%smallestrootindex

    indices(1)=thisroot

    j=0
    do i = 1, 3
       if(i==thisroot)cycle
       j=j+1
       indices(j+1)=i
       if(j==2)exit
    end do
    theroot = - hsign * sqrt((a-u)/((a-vpi%uvw(indices(2)))*(a-vpi%uvw(indices(3)))))

    thishp = thisdhdalike * adotterm + & 
         yprime * mreal(xyu*theroot)

    Hprime_amax_mkl = thishp
  end function Hprime_amax_mkl

  function Hprime_amax_mk(vpi,a,yprime,hsign,adotterm,aprime)
    implicit none
    real(dl), intent(in) :: a,yprime,hsign,adotterm, aprime
    type(voidparams) :: vpi
    real(dl) :: Hprime_amax_mk
    real(dl) :: thishp, thish
    real(dl) :: thisadotdot_overa
    real(dl) :: thisadot
    real(dl) :: om, ok, thiszero


    thisadotdot_overa =  - 0.5_dl * vpi%xyz(1) / a**3 !+ vpi%xyz(3)=0
    
    thish = hsign * HoverH0(vpi%xyz,a)

    thisadot = a * thish

    om = vpi%workxyz(1)
    ok = vpi%workxyz(2)

    thiszero = 0.5_dl*yprime * thisadot * (( -1.5_dl * om / a / ok**2_dl + 1._dl/ok) ) / a
! ==    thiszero = 0.5_dl*yprime * thisadot *  -1.5_dl * om / a**2 / ok**2_dl + 0.5_dl*yprime * thisadot * 1._dl/ok / a

    thishp = thisadotdot_overa * adotterm + thiszero !+ & 
!         hsign*yprime * 0.5_dl / a * mreal(theroot)*xyu

    thishp = thishp - thish / a * aprime

    Hprime_amax_mk = thishp 

    
!    Hprime_amax_mk = 0.

  end function Hprime_amax_mk



end module lltb_background_limits
