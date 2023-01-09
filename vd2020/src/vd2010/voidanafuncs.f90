  subroutine voidtestanafunc()
    implicit none
    integer :: i, k
    real(dl) :: r, t, u(3), S, Sn, Sd, Sdn


    call setrtopini
    VP%L = VP%rtopini / 3.d0

    k = 100

    do i = 1,k
       r = VP%L*dble(i-1) / dble(k-1)
       t = VP%t0*(1.d0-dble(k-i+1) / dble(2*k))
       u = voidu(r,t)
       S = avoidS(r,t,u)
       Sd = avoidSdot(r,t,u)
! include voidnumfuncs.f90 in mainfile for this:
!       Sn = nvoidS(r,t,u)
!       Sdn = nvoidSdot(r,t,u)
       
       write(*,'(10ES15.6)')r,t,S,Sn,Sd,Sdn,u(:)
    end do

    stop 'voidtestanafunc done.'

  end subroutine voidtestanafunc

  function voidtofa(r,a)
    implicit none
    real(dl) :: voidtofa
    real(dl) :: r,a

    ! testing:
    real(dl) :: myval
    call lltb_functions_a(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=voidtofa,a=a)
    return
    
  end function voidtofa
  

  function voida(r,t,u1)
    implicit none
    real(dl) :: voida
    real(dl) :: r,t,u(3)
    real(dl), optional :: u1(3)

! testing:
    real(dl) :: myval
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,a=voida)
return

    if(present(u1))then
       u=u1
    else
       u=voidu(r,t)
    end if

    if(u(1).gt.0.d0)then
       voida= 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcosh(u(1)) - 1.d0  ) 
    else if(u(2).gt.0.d0) then
       voida= 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcos(u(2)) - 1.d0  ) 
    else if(any(u.lt.0.d0)) then
       voida = VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
    end if


    if(abs(myval/voida-1.0_dl)>1.)then
       call VoidErrorReport()
       write(*,'(3(A,ES16.7))')'a oldcode:',voida,' newcode:',myval,' rel diff:',abs(myval/voida-1.0_dl)
    end if

    if(isnan(voida).or.voida.gt.1.d29)then
          voida = VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
          if(voidfeedback.gt.2)then
             write(*,*)'Applying analytical approximation for voida for r close to L.'
             write(*,'(A,5E14.5)')' voida,r,t,L:', voida,r,t,VP%L
          end if
          return
       end if
   end function voida

   function voidaout(t,u1)
    implicit none
    real(dl) :: voidaout
    real(dl) :: r,t,u(3)
    real(dl), optional :: u1(3)

! testing:
    real(dl) :: myval

    r= 1.2*VP%L 

    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,a=voidaout)
return

    if(present(u1))then
       u=u1
    else
       u=voidu(r,t)
    end if

    if(u(1).gt.0.d0)then
       voidaout= 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcosh(u(1)) - 1.d0  ) 
    else if(u(2).gt.0.d0) then
       voidaout= 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcos(u(2)) - 1.d0  ) 
    else if(any(u.lt.0.d0)) then
       voidaout = VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
    end if


    if(abs(myval/voidaout-1.0_dl)>1.)then
       call VoidErrorReport()
       write(*,'(3(A,ES16.7))')'a oldcode:',voidaout,' newcode:',myval,' rel diff:',abs(myval/voidaout-1.0_dl)
    end if

    if(isnan(voidaout).or.voidaout.gt.1.d29)then
          voidaout = VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
          if(voidfeedback.gt.2)then
             write(*,*)'Applying analytical approximation for voidaout for r close to L.'
             write(*,'(A,5E14.5)')' voidaout,r,t,L:', voidaout,r,t,VP%L
          end if
          return
       end if
   end function voidaout

  function voidaH(r,t,u1)
    implicit none
    real(dl) :: voidaH
    real(dl) :: r,t,u(3)
    real(dl), optional :: u1(3)
! testing:
    real(dl) :: myval,myval2
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,a=myval,H=myval2)
    voidaH=myval*myval2
    return

    if(present(u1))then
       u=u1
    else
       u=voidu(r,t)
    end if

    if(u(1).gt.0.d0)then
       voidaH= dsqrt(2.d0 * voidkofr(r)) * VP%Mtilde * dsinh(u(1)) / ( dcosh(u(1)) - 1.d0  ) 
    else if (u(2).gt.0.d0)then
       voidaH= - dsqrt(2.d0 * dabs(voidkofr(r))) * VP%Mtilde * dsin(u(2)) / ( dcos(u(2)) - 1.d0  ) 
    else if (any(u.lt.0.d0)) then
       voidaH =2.d0/3.d0* VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(-1.d0/3.d0)
    end if
 
    if(isnan(voidaH).or.voidaH.gt.1.d29)then
          voidaH =2.d0/3.d0* VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(-1.d0/3.d0)
          if(voidfeedback.gt.2)then
             write(*,*)'Applying analytical approximation for voidaH for r close to L.'
             write(*,'(A,5E14.5)')' voidaH,r,t,L:', voidaH,r,t,VP%L
          end if
          return
       end if

    write(*,'(3(A,ES16.7))')'aH oldcode:',voidaH,' newcode:',myval,' rel diff:',abs(myval/voidaH-1.0_dl)

   end function voidaH




  function voidH(r,t,u_opt)
    implicit none
    real(dl) :: voidH
    real(dl) :: r,t,u(3)
    real(dl) :: k
    real(dl), optional :: u_opt(3)
! testing:
    real(dl) :: myval

    
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,H=voidH)
    
    return
    
! Next is never executed anymore

    if(present(u_opt))then
       u=u_opt
    else
       u=voidu(r,t)
    end if

!    if((r.ge.VP%L .and. VP%kb_om.eq.0.d0) .or. any(u(1:2).lt.0.d0))then
    if((r.ge.VP%LTB_maxr .and. VP%kb_om.eq.0.d0) .or. any(u(1:2).lt.0.d0))then
       ! FLRW
       voidH=2.d0/3.d0/t
       return
    end if

    k = voidkofr(r)
    
    if(u(1).gt.0.d0)then
       voidH = 3.d0 * dabs(k)**(1.5d0) * VP%Mtilde / dsqrt(2.d0) / pi &
            * dsinh(u(1)) &
            / (dcosh(u(1))-1.d0)**2
    else if(u(2).gt.0.d0) then
       voidH = 3.d0 * dabs(k)**(1.5d0) * VP%Mtilde / dsqrt(2.d0) / pi &
            * dsin(u(2)) &
            / (dcos(u(2))-1.d0)**2
       if(voidH.lt.0.d0)then
         ! stop 'error in function voidH: H<0.'
          write(*,*) 'error in function voidH: H<0.'
          VoidError = .true.
          return
       end if
    else
       write(*,*) 'error in u, function voidH.'
       VoidError=.true.
       return
    end if

    write(*,'(3(A,ES16.7))')'H oldcode:',voidH,' newcode:',myval,' rel diff:',abs(myval/voidH-1.0_dl)
  end function voidH


  function voidHaout(t,u_opt)
    implicit none
    real(dl) :: voidHaout
    real(dl) :: r,t,u(3)
    real(dl) :: k
    real(dl), optional :: u_opt(3)
! testing:
    real(dl) :: myval

    r= 1.2*VP%L 

    
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,H=voidHaout)
    
    return
    
! Next is never executed anymore

    if(present(u_opt))then
       u=u_opt
    else
       u=voidu(r,t)
    end if

!    if((r.ge.VP%L .and. VP%kb_om.eq.0.d0) .or. any(u(1:2).lt.0.d0))then
    if((r.ge.VP%LTB_maxr .and. VP%kb_om.eq.0.d0) .or. any(u(1:2).lt.0.d0))then
       ! FLRW
       voidHaout=2.d0/3.d0/t
       return
    end if

    k = voidkofr(r)
    
    if(u(1).gt.0.d0)then
       voidHaout = 3.d0 * dabs(k)**(1.5d0) * VP%Mtilde / dsqrt(2.d0) / pi &
            * dsinh(u(1)) &
            / (dcosh(u(1))-1.d0)**2
    else if(u(2).gt.0.d0) then
       voidHaout = 3.d0 * dabs(k)**(1.5d0) * VP%Mtilde / dsqrt(2.d0) / pi &
            * dsin(u(2)) &
            / (dcos(u(2))-1.d0)**2
       if(voidHaout.lt.0.d0)then
         ! stop 'error in function voidH: H<0.'
          write(*,*) 'error in function voidHaout: H<0.'
          VoidError = .true.
          return
       end if
    else
       write(*,*) 'error in u, function voidHaout.'
       VoidError=.true.
       return
    end if

    write(*,'(3(A,ES16.7))')'H oldcode:',voidHaout,' newcode:',myval,' rel diff:',abs(myval/voidHaout-1.0_dl)
  end function voidHaout


  function voidR(r,t,u1)
    implicit none
    real(dl) :: voidR
    real(dl) :: r,t,u(3)
    real(dl), optional :: u1(3)
    real(dl), parameter :: maxu(2)=(/acosh(max_double),asinh(max_double)/)
! testing:
    real(dl) :: myval
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,Rltb=voidR)
return

    if(present(u1))then
       u=u1
    else
       u=voidu(r,t)
    end if
    
!    if(r.gt.VP%L.and. VP%kb_om.eq.0.d0)then
    if(r.gt.VP%LTB_maxr.and. VP%kb_om.eq.0.d0)then
       voidR= r * VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
       return
    end if
    
    if(any(u(1)>maxu).or.any(u(2)>maxu))then
       voidR = 1.d31
       voiderror = .true.
       return
    end if

    if(u(1).gt.0.d0)then
       voidR= r * 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcosh(u(1)) - 1.d0  ) 
    else if(u(2).gt.0.d0)then
       voidR= r * 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcos(u(2)) - 1.d0  ) 
    else if (any(u.lt.0.d0)) then
       voidR = r * VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
    end if
       
 
    if(isnan(voidR).or.voidR.gt.1.d29)then
          voidR = r * VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
          if(voidfeedback.gt.2)then
             write(*,*)'Applying analytical approximation for voidR for r close to L.'
             write(*,'(A,5E14.5)')' voidR,r,t,L:', voidR,r,t,VP%L
          end if
          return
       end if
       write(*,'(3(A,ES16.7))')'R oldcode:',voidR,' newcode:',myval,' rel diff:',abs(myval/voidR-1.0_dl)
   end function voidR


   function Debug_my_aprime(r,t,u1)
    implicit none
    real(dl) Debug_my_aprime
    real(dl) :: r,t,u(3)
    real(dl) :: k, thisRp, newap, my_Err
    real(dl), optional :: u1(3)

    if(present(u1))then
       u=u1
    else
       u=voidu(r,t)
    end if

    Debug_my_aprime = ( avoidRp(r,t,u) - voida(r,t,u) ) / r
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,ap=newap)

    if(Debug_my_aprime==0.d0 .and. newap==0.d0)then
       my_Err = 0.d0
    else
       my_err = abs(newap/Debug_my_aprime-1.d0)
    end if

    writE(*,'(3(A,ES15.4))')' aprima. old:', Debug_my_aprime, ' new:', newap, ' rel: ', my_err
    

  end function Debug_my_aprime

  function avoidS(r,t,u1)
    implicit none
    real(dl) avoidS
    real(dl) :: r,t,u(3)
    real(dl) :: k, thisRp
    real(dl), optional :: u1(3)
! testing:
    real(dl) :: myval
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,S=avoidS)
return

    if(present(u1))then
       u=u1
    else
       u=voidu(r,t)
    end if

    k=voidkofr(r)

    thisRp= avoidRp(r,t,u)

    if( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2 .le. 0.d0) then
       avoidS = 1.d30
       voiderror = .true.
!    else if( thisRp .le. 0.d0 .or. isnan(thisRp)) then
    else if(isnan(thisRp)) then
       avoidS = 1.d30
       voiderror = .true.
    else

       avoidS = thisRp / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
    end if

!       write(*,'(A,E14.5,A,E14.5,A,E14.5,A,E14.5)')'voidS:',avoidS,' voidRp:',thisRp,' r:',r,' t:',t

    ! Error tracking:
!    DER%S = DER%Rp / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
    if(abs(myval/avoidS-1.0_dl)>1.)write(*,'(3(A,ES16.7))')'S oldcode:',avoidS,' newcode:',myval,' rel diff:',abs(myval/avoidS-1.0_dl)
    if(abs(myval/avoidS-1.0_dl)>1.)then
       call VoidErrorReport()
       write(*,*)'at r:',r
       write(*,*)'at t:', t
       write(*,*)'Calling extra instances of voida and voidH.'
       myval=voida(r,t)
       myval=voidH(r,t)
       write(*,*)'Calling Debug_my_aprime.'
       myval=Debug_my_aprime(r,t,u)
       write(*,*)'Extra calling done.'
    end if


  end function avoidS

  function avoidSout(t,u1)
    implicit none
    real(dl) avoidSout
    real(dl) :: r,t,u(3)
    real(dl) :: k, thisRp
    real(dl), optional :: u1(3)
! testing:
    real(dl) :: myval

    r= 1.2*VP%L

    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,S=avoidSout)
return

    if(present(u1))then
       u=u1
    else
       u=voidu(r,t)
    end if

    k=voidkofr(r)

    thisRp= avoidRp(r,t,u)

    if( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2 .le. 0.d0) then
       avoidSout = 1.d30
       voiderror = .true.
!    else if( thisRp .le. 0.d0 .or. isnan(thisRp)) then
    else if(isnan(thisRp)) then
       avoidSout = 1.d30
       voiderror = .true.
    else

       avoidSout = thisRp / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
    end if

!       write(*,'(A,E14.5,A,E14.5,A,E14.5,A,E14.5)')'voidS:',avoidSout,' voidRp:',thisRp,' r:',r,' t:',t

    ! Error tracking:
!    DER%S = DER%Rp / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
    if(abs(myval/avoidSout-1.0_dl)>1.)write(*,'(3(A,ES16.7))')'S oldcode:',avoidSout,' newcode:',myval,' rel diff:',abs(myval/avoidSout-1.0_dl)
    if(abs(myval/avoidSout-1.0_dl)>1.)then
       call VoidErrorReport()
       write(*,*)'at r:',r
       write(*,*)'at t:', t
       write(*,*)'Calling extra instances of voida and voidH.'
       myval=voida(r,t)
       myval=voidH(r,t)
       write(*,*)'Calling Debug_my_aprime.'
       myval=Debug_my_aprime(r,t,u)
       write(*,*)'Extra calling done.'
    end if


  end function avoidSout

  function avoidSdot(r,t,u_opt)
    implicit none
    real(dl) :: avoidSdot
!    real(dl) :: voidRpdot
    real(dl) :: r,t,u(3)
    real(dl) :: Rpd !, Rp1, Rp2
!    real(dl) :: distance
    real(dl), parameter :: myprecision = numerical_derivatives_precision !myprecision = 1.d-7
    integer, parameter :: mymax = 10000
!    integer :: counter
!    real(dl) :: sinhu_u, k, prev(3)
    real(dl) :: k

    real(dl), optional :: u_opt(3)
! Numderical recipes:
    real(dl), parameter :: con = 1.4, con2=con*con, big=1.d30, safe = 2.d0
    integer, parameter :: Ntab=10
!    integer :: nri, nrj
!    real(dl) :: err, errt, fac, nra(ntab,ntab), thisdfdx
! END NR
! testing:
    real(dl) :: myval
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,Sd=avoidSdot)
return

    if(present(u_opt))then
       u=u_opt
    else
       u=voidu(r,t)
    end if

    ! Sdot / S = R'dot / R'.
    ! a = R / r
    ! R'dot = r a'dot + adot
    !       = r (aH)' + aH
    ! We know a and H analytically as function of t and r.

!    if((r.ge.VP%L .and. VP%kb_om.eq.0.d0) .or.  any(u(1:2).lt.0.d0))then
    if((r.ge.VP%LTB_maxr .and. VP%kb_om.eq.0.d0) .or.  any(u(1:2).lt.0.d0))then
       ! FLRW: Sdot / S = H(t)
       avoidSdot=aVoidSDotoverS(r,t,u) * avoidS(r,t,u) !2.d0/3.d0/t
       return
    end if

    ! And the next line fucks every thing down (i.e. compensating fucking up):
    k = voidkofr(r)
    ! do not forget next time!

    Rpd = aVoidRpdot(r,t,u)

    if( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2 <= 0.d0) then
       avoidSdot = 1.d30
       voiderror = .true.
    else
       avoidSdot=Rpd / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
    end if
    !write(*,*)voidSdotoverS
    !stop'Sdot'

    ! Error tracking:
!    DER%Sdot = DER%Rpdot / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
    if(abs(myval/avoidSdot-1.0_dl)>1.)write(*,'(3(A,ES16.7))')'Sdot oldcode:',avoidSdot,' newcode:',myval,' rel diff:',abs(myval/avoidSdot-1.0_dl)
  end function avoidSdot

  function avoidSdotoverS(r,t,u_opt)
    implicit none
    real(dl) :: avoidSdotoverS
!    real(dl) :: voidRpdot
    real(dl) :: r,t,u(3)
    real(dl) :: Rpd !, Rp1, Rp2
!    real(dl) :: distance
    real(dl), parameter :: myprecision = numerical_derivatives_precision !myprecision = 1.d-7
    integer, parameter :: mymax = 10000
!    integer :: counter
!    real(dl) :: sinhu_u, prev(3)

    real(dl), optional :: u_opt(3)
! Numderical recipes:
    real(dl), parameter :: con = 1.4, con2=con*con, big=1.d30, safe = 2.d0
    integer, parameter :: Ntab=10
!    integer :: nri, nrj
!    real(dl) :: err, errt, fac, nra(ntab,ntab), thisdfdx
! END NR
! testing:
    real(dl) :: myval,myval2
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,Sd=myval,S=myval2)
    avoidSdotoverS=myval/myval2
    if(.not.isnan(avoidSdotoverS))then
       return
    else
       write(*,*)'nan in sdotovers.'
       write(*,*)'S:',myval2
       write(*,*)'Sdot:',myval

    end if


    if(present(u_opt))then
       u=u_opt
    else
       u=voidu(r,t)
    end if

    ! Sdot / S = R'dot / R'.
    ! a = R / r
    ! R'dot = r a'dot + adot
    !       = r (aH)' + aH
    ! We know a and H analytically as function of t and r.

!    if((r.ge.VP%L .and. VP%kb_om.eq.0.d0) .or.  any(u(1:2).lt.0.d0))then
    if((r.ge.VP%LTB_maxr .and. VP%kb_om.eq.0.d0) .or.  any(u(1:2).lt.0.d0))then
       ! FLRW: Sdot / S = H(t)
       avoidSdotoverS=2.d0/3.d0/t
       return
    end if

    Rpd = aVoidRpdot(r,t,u)

    avoidSdotoverS=Rpd / avoidRp(r,t,u)
    !write(*,*)voidSdotoverS
    !stop'Sdot'
    write(*,'(3(A,ES16.7))')'Sd/S oldcode:',avoidSdotoverS,' newcode:',myval,' rel diff:',abs(myval/avoidSdotoverS-1.0_dl)
!stop 'fuck it.'
  end function avoidSdotoverS



  function avoidRp(r,t,u_opt)
    implicit none
    real(dl) :: avoidRp
    real(dl) :: r,t,u(3)
    real(dl) :: sinhu_u!, prev(3)
    real(dl), optional :: u_opt(3)
    real(dl) :: chu, shu, kp, k, cq, sq
    real(dl), parameter :: maxu(2)=(/dacosh(max_double),dasinh(max_double)/)
! testing:
    real(dl) :: myval
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,Rp=avoidRp)
return

    if(present(u_opt))then
       u=u_opt
    else
       u=voidu(r,t)
    end if

    ! Analytical approximation if voidu uses it too:
!    sinhu_u = t * 3.d0 * voidkofr(r)**(1.5d0) * VP%Mtilde / dsqrt(2.d0)/pi
!    if((dabs(6.d0*sinhu_u)**(2.d0/3.d0).lt.u_precision).or.r.ge.VP%L)then 
!    if(all(u.lt.0.d0).or. (r.ge.VP%L .and. VP%kb_om.eq.0.d0))then 
    if(all(u.lt.0.d0).or. (r.ge.VP%LTB_maxr .and. VP%kb_om.eq.0.d0))then 
       avoidRp = VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
       !   write(*,*)voidRp,r,t
       !   stop 'r=0'

!       if(voidfeedback.gt.4 .or. (voidfeedback.gt.2 .and. r.lt. VP%L))then
       if(voidfeedback.gt.4 .or. (voidfeedback.gt.2 .and. r.lt. VP%LTB_maxr))then
          write(*,*)'Applying analytical approximation for avoidRp.'
          write(*,'(A,5E14.5)')' avoidRp,r,t,sinhu_u, L:', avoidRp,r,t,sinhu_u, VP%L
       end if
       return
    end if

!!    Reminder of the use of u:
!!    if(u(1).gt.0.d0)then
!!       voidR= r * 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcosh(u(1)) - 1.d0  ) 
!!    else if(u(2).gt.0.d0)then
!!       voidR= r * 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcos(u(2)) - 1.d0  ) 
!!    else if (any(u.lt.0.d0)) then
!!       voidR = r * VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
!!    end if
!!       

    if(u(1)>0.d0 .and. all(u(1)<maxu))then ! u is real, sinh's and cosh's
       
       if(u(3)>0.d0)then ! the full monty
          chu = dcosh(u(1))
          shu = dsinh(u(1))
          k = voidkofr(r)
          kp = voiddkdr(r)
          
          avoidRp = 2.d0 * pi / 3.d0 / k * ( &
               (chu - 1.d0)*(1.d0 - r*kp/k) + &
               3.d0 * r * kp / 2.d0 / k * shu * (shu - u(1)) / (chu - 1.d0))
       else ! small u approximation
 
          k = voidkofr(r)
          kp = voiddkdr(r)
          avoidRp = pi * u(1)**2 / 3.d0 / k + (5.d0 * pi * k + 3.d0 * pi * r * kp) / 180.d0 / k**2 * u(1)**4 ! + O(u(1)^6)
       end if

    else if (u(2)>0.d0 .and. all(u(2)<maxu)) then ! u is imaginary, sin's and cos's

       if(u(3)>0.d0)then ! the full monty

          cq = dcos(u(2))
          sq = dsin(u(2))
          k = voidkofr(r)
          kp = voiddkdr(r)
          
          avoidRp = 2.d0 * pi / 3.d0 / k * ( &
               (cq - 1.d0)*(1.d0 - r * kp / k) + &
               3.d0 * r * kp / 2.d0 / dabs(k) * sq * &
               (sq - u(2)) / (cq - 1.d0))
               

       else ! small u approximation

          k = voidkofr(r)
          kp = voiddkdr(r)

          avoidRp = (pi*(k*kp*r + (-k + kp*r)*dAbs(k))*u(2)**2)/(3.d0*k**2*dAbs(k)) + (pi*(-8.d0*k*kp*r + 5.d0*(k - kp*r)*dAbs(k))*u(2)**4)/(180.d0*k**2*dAbs(k))

       end if
    else! other situations do not occur at this stage.
       if(VP%omk<1.d0)then
          voiderror = .true.
          avoidRp = 1.d30
          return
       else if( any(u(2)>maxu).or.  any(u(1)>maxu) )then
          voiderror = .true.
          avoidRp = 1.d30
          return
       else
          
          writE(*,*)u
          write(*,'(A,4ES18.9)')'zB, kmax, omk, H_out:',VP%zB,VP%kmax, VP%omk, VP%H0
          call VoidErrorReport()
          stop 'This point should not occur in avoidRp'
       end if
    end if

!    DER%Rp = err
    if(abs(myval/avoidRp-1.0_dl)>1.)write(*,'(3(A,ES16.7))')'Rp oldcode:',avoidRp,' newcode:',myval,' rel diff:',abs(myval/avoidRp-1.0_dl)
  end function avoidRp


!!$  function voidudot(r,t,u_opt)
!!$    implicit none
!!$    real(dl) :: voidudot
!!$    real(dl) :: r,t,u(3)
!!$    real(dl) :: sinhu_u, prev(3)
!!$    real(dl), optional :: u_opt(3)
!!$    real(dl) :: chu, k
!!$
!!$    if(present(u_opt))then
!!$       u=u_opt
!!$    else
!!$       u=voidu(r,t)
!!$    end if
!!$
!!$    if(u(1)>0.d0)then ! u is real, sinh's and cosh's
!!$       
!!$
!!$       chu = dcosh(u(1))
!!$       k = voidkofr(r)
!!$
!!$       voidudot = 3.d0 * VP%Mtilde / 2.d0**0.5d0 / pi * k**(1.5d0) / (chu - 1.d0)
!!$
!!$    else if (u(2)>0.d0) then ! u is imaginary, sin's and cos's
!!$       stop 'not written yet: voidudot'
!!$    else! other situations do not occur at this stage.
!!$       stop 'This point should not occur in voidudot'
!!$    end if 
!!$
!!$
!!$  end function voidudot
!!$

  function avoidRpdot(r,t,u_opt)
    implicit none
    real(dl) :: avoidRpdot
    real(dl) :: r,t,u(3)
    real(dl) :: sinhu_u!, prev(3)
    real(dl), optional :: u_opt(3)
!    real(dl) :: chu, shu, kp, k, cq, sq
    real(dl) :: kp, k
    real(dl), parameter :: maxu(2)=(/dacosh(max_double),dasinh(max_double)/)
! testing:
    real(dl) :: myval
    call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,Rpd=avoidRpdot)
return

    if(present(u_opt))then
       u=u_opt
    else
       u=voidu(r,t)
    end if

    ! Analytical approximation if voidu uses it too:
!    sinhu_u = t * 3.d0 * voidkofr(r)**(1.5d0) * VP%Mtilde / dsqrt(2.d0)/pi
!    if((dabs(6.d0*sinhu_u)**(2.d0/3.d0).lt.u_precision).or.r.ge.VP%L)then 
!    if(all(u.lt.0.d0).or. (r.ge.VP%L .and. VP%kb_om.eq.0.d0))then 
    if(all(u.lt.0.d0).or. (r.ge.VP%LTB_maxr .and. VP%kb_om.eq.0.d0))then 
       avoidRpdot = VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
       !   write(*,*)voidRp,r,t
       !   stop 'r=0'
write(*,*)u
write(*,*)r,t
call VoidErrorReport()

       stop 'write analytical approximation for avoidrpdot.'
       if(voidfeedback.gt.2)then
          write(*,*)'Applying analytical approximation for avoidRpdot.'
          write(*,'(A,5E14.5)')' avoidRpdot,r,t,sinhu_u, L:', avoidRpdot,r,t,sinhu_u, VP%L
       end if
       return
    end if

!!    Reminder of the use of u:
!!    if(u(1).gt.0.d0)then
!!       voidR= r * 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcosh(u(1)) - 1.d0  ) 
!!    else if(u(2).gt.0.d0)then
!!       voidR= r * 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcos(u(2)) - 1.d0  ) 
!!    else if (any(u.lt.0.d0)) then
!!       voidR = r * VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
!!    end if
!!       
    
    if(u(1).gt.0d0 .and. all(u(1)<maxu))then
       if(u(3).gt.0.d0 )then ! numerical evaluation worked:
          !       avoidRpdot = avoidRpdotoverudot(r,t,u) * voidudot(r,t,u)
          k = voidkofr(r)
          kp = voiddkdr(r)
          
          avoidRpdot = &
               dSinh(u(1)*0.5d0)**2 / 2.d0**0.5d0 / k**0.5d0 / (dcosh(u(1))-1.d0)**3 * &
               (2.d0*k*(-2.d0 * dsinh(u(1))+ dsinh(2.d0 * u(1))) + &
               r * kp * ( 6.d0 * u(1) - 8.d0 * dsinh(u(1)) + dsinh(2.d0*u(1))))
          avoidRpdot = avoidRpdot  * VP%Mtilde
       else ! small u expansion:
          
          k = voidkofr(r)
          kp = voiddkdr(r)
          avoidRpdot = 2.d0**1.5d0 * k**0.5d0 / u(1) + &
               (5.d0 * k + 6.d0 * r * kp) / 15.d0 / 2.d0**0.5d0 / k**0.5d0 * u(1)
          avoidRpdot = avoidRpdot  * VP%Mtilde
          
       end if
    else if(u(2)>0 .and. all(u(2)<maxu)) then ! imaginary u
       if(u(3).gt.0.d0)then ! numerical evaluation worked:
          !       avoidRpdot = avoidRpdotoverudot(r,t,u) * voidudot(r,t,u)
          k = voidkofr(r)
          kp = voiddkdr(r)
          
          avoidRpdot = (dSqrt(dAbs(k))*dSin(u(2)/2.d0)**2*(16.d0*(k - kp*r)*dAbs(k)*dCos(u(2)/2.d0)*dSin(u(2)/2.d0)**3 + &
                3.d0*k*kp*r*(-4.d0*dSin(u(2)) + dSin(2.d0*u(2)) + 2.d0*u(2))))/(dSqrt(2.d0)*k**2*(-1.d0 + dCos(u(2)))**3)

          avoidRpdot = avoidRpdot  * VP%Mtilde
       else ! small u expansion:
          
          k = voidkofr(r)
          kp = voiddkdr(r)
          avoidRpdot = (-2.d0*dSqrt(2.d0)*dSqrt(dAbs(k))*(-(k*kp*r) + k*dAbs(k) - kp*r*dAbs(k)))/(k**2*u(2)) + &
               (dSqrt(dAbs(k))*(-924.d0*k*kp*r + 420.d0*k*dAbs(k) - 420.d0*kp*r*dAbs(k))*u(2))/(1260.d0*dSqrt(2.d0)*k**2) + &
               (dSqrt(dAbs(k))*(-31.d0*k*kp*r + 7.d0*k*dAbs(k) - 7.d0*kp*r*dAbs(k))*u(2)**3)/(1260.d0*dSqrt(2.d0)*k**2)
          ! long live FortranForm in Mathematica..
          avoidRpdot = avoidRpdot  * VP%Mtilde

       end if
    else
       voiderror = .true.
       avoidRpdot = 1.d30
       return
       
    end if

    !    DER%Rp = err
    if(abs(myval/avoidRpdot-1.0_dl)>1.)write(*,'(3(A,ES16.7))')'Rpd oldcode:',avoidRpdot,' newcode:',myval,' rel diff:',abs(myval/avoidRpdot-1.0_dl)
  end function avoidRpdot
  
!!$
!!$  function avoidRpdotoverudot(r,t,u_opt)
!!$    implicit none
!!$    real(dl) :: avoidRpdotoverudot
!!$    real(dl) :: r,t,u(3)
!!$    real(dl) :: sinhu_u, prev(3)
!!$    real(dl), optional :: u_opt(3)
!!$    real(dl) :: chu, shu, kp, k
!!$
!!$    if(present(u_opt))then
!!$       u=u_opt
!!$    else
!!$       u=voidu(r,t)
!!$    end if
!!$
!!$    ! Analytical approximation if voidu uses it too:
!!$!    sinhu_u = t * 3.d0 * voidkofr(r)**(1.5d0) * VP%Mtilde / dsqrt(2.d0)/pi
!!$!    if((dabs(6.d0*sinhu_u)**(2.d0/3.d0).lt.u_precision).or.r.ge.VP%L)then 
!!$    if(all(u.lt.0.d0).or. (r.ge.VP%L .and. VP%kb_om.eq.0.d0))then 
!!$       stop 'write analytical approximation rpdotoverudot.'
!!$       if(voidfeedback.gt.2)then
!!$          write(*,*)'Applying analytical approximation for avoidRpdotoverudot.'
!!$          write(*,'(A,5E14.5)')' avoidRpdotoverudot,r,t,sinhu_u, L:', avoidRpdotoverudot,r,t,sinhu_u, VP%L
!!$       end if
!!$       return
!!$    end if
!!$
!!$!    if(.not.(any(u(1:2).lt.0.d0)).and.(u(3).lt.0.d0))then
!!$!       write(*,*)' small u approximation.'
!!$!    else
!!$!       write(*,*)' Numerical u.'
!!$!    end if
!!$
!!$!!    Reminder of the use of u:
!!$!!    if(u(1).gt.0.d0)then
!!$!!       voidR= r * 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcosh(u(1)) - 1.d0  ) 
!!$!!    else if(u(2).gt.0.d0)then
!!$!!       voidR= r * 2.d0 * pi / 3.d0 / voidkofr(r) * ( dcos(u(2)) - 1.d0  ) 
!!$!!    else if (any(u.lt.0.d0)) then
!!$!!       voidR = r * VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
!!$!!    end if
!!$!!       
!!$
!!$    if(u(1)>0.d0)then ! u is real, sinh's and cosh's
!!$       
!!$
!!$       chu = dcosh(u(1))
!!$       shu = dsinh(u(1))
!!$       k = voidkofr(r)
!!$       kp = voiddkdr(r)
!!$
!!$       avoidRpdotoverudot = 2.d0 * pi / 3.d0 / k * (&
!!$            (1.d0 + r * kp / 2.d0 / k) * shu + &
!!$            3.d0 / 2.d0 * r * kp / k * (shu - u(1))/(chu - 1.d0) * &
!!$            (chu - (shu - u(1))**2 / (chu - 1.d0)))
!!$
!!$
!!$    else if (u(2)>0.d0) then ! u is imaginary, sin's and cos's
!!$       stop 'not written yet: avoidRpdotoverudot'
!!$    else! other situations do not occur at this stage.
!!$       stop 'This point should not occur in avoidRpdotoverudot'
!!$    end if 
!!$
!!$!    DER%Rp = err
!!$  end function avoidRpdotoverudot
!!$
!!$

  ! y(1) = t
  ! y(2) = dt/du = - (z + 1)
  ! y(3) = k
  ! du = - dz
  subroutine ltb_aderivsdr(v,y,dy)
    implicit none
    real(dl), intent(in) :: v
    real(dl), intent(in) :: y(:)
    real(dl), intent(out) :: dy(:)
    real(dl) :: u(3),t,r, xx, thisS, thisSdot, thisk
    ! new off-center
    real(dl) :: drdl, thisa
    logical :: fatalerror
    ! end

!    logical, save :: warned=.false.
!    if(size(y).ne.3)stop 'Wrong dimension in ltb_derivs.'
    if(abs(size(y)-4).ne.1)stop 'Wrong dimension in ltb_derivs.'
    ! this says y==(3 or 5).
    
    if(voiderror)then
      dy=1.d300
      return ! Don't go on with errors!
    end if
    
    r=v!y(3)
    t=y(1)
    u = voidu(r,t)


    fatalerror=.false.
!    if(t.lt.0.d0 .or. t.gt.1.e300_dl .or. t>1.d2*vp%t0)then ! Why again t<1.d2t0???
    if(t.lt.0.d0 .or. t.gt.1.e200_dl )then
          if(voidfeedback>2)call VoidRedMessage('t = unphysical, <0 or >1.d200.')
          VoidError = .true.
    end if

!  WRITE(0,*)'Hoi23',y,dy,VP%t0,voiderror


    ! Check that we do not exceed the maximum radius of a closed universe:
    thisk=voidkofr(R)
    if(1._dl + 2._dl * r**2 * VP%Mtilde**2 * thisk .le. 0.d0)then
      voiderror=.true.
      fatalerror=.true.
      if(voidfeedback>2)call VoidRedMessage('curvature term in metric < 0')
    end if  

!    write(*,'(2(A,ES15.8))')'t:',y(1),'   z:',-y(2)-1._dl
!    write(*,*) voiderror
    
    if(.not.voiderror)then
       
       call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=r,t=t,Sd=thisSdot,S=thisS, a=thisa)
        if(.not. thisa>0.d0)VoidError=.true. ! New: detecting collapsed solution from wlltb_background.
       if(abs(abs(vp%tdir)-1.d0)>1.d-30)then
        call VoidErrorReport()
        write(0,*)'VP%tdir == 0'
        stop
      end if
             
       dy(1) = - thisS * vp%tdir
       dy(2) = thisSdot * y(2)* vp%tdir


!       write(*,*)-y(2)-1._dl,-dy(2)/ vp%tdir, VP%H0 * ( (1.-VP%omv-VP%omk)*(-y(2))**3 + VP%omk*(-y(2))**2 + VP%omv*(-y(2))**(3*(1+VP%w)))**0.5d0

 ! if(voidtesting)write(*,*)dy(1:2),y(1:2),t,r
!!$    dy(1) = - avoidS(r,t,u)
!!$    dy(2) = avoidSdot(r,t,u) * y(2)

       ! new off-center:
       if(size(dy)>4)then
          drdl = -y(2) / thisS
          dy(4) = y(5) / drdl
          dy(5) = - y(4) * 4._dl * pi * VP%Mtilde**2 / sqrt(1._dl + 2._dl * thisk * r**2 * VP%Mtilde**2) * (- y(2))/thisa**2 
       end if


       if(dy(2)* vp%tdir.gt.0.d0)then
          if(VoidDisallowBlueshift)then
             call VoidBigAlert('You chose not to allow for blueshift, and now we start rejecting models with blueshift.')
             
             VoidError = .true.
          end if
       end if
       
       
       ! Simple check of physicalness of parameters:
       if(VP%tdir*dy(1).gt.0.d0)then
          if(voidfeedback.gt.0)call VoidAlertAtMostOnceEveryModel('Rejecting unphysical model: shell density diverges, clock runs backwards / shell crossing.')
          if(voidfeedback.gt.1)write(*,'(A)')' Rejecting unphysical model: shell density diverges, clock runs backwards / shell crossing.'
!          writE(*,*)'L:',VP%L
!          write(*,*)'S:',thisS
!          write(*,*)'r:',r
          dy=0.d0
          dy=1.d300
          voiderror = .true.
          ! write(*,*)'David suspects that it is a unphysical point.'
          write(*,*)'Rejecting model due to numerical trouble.'
          writE(*,'(A,5E17.8)')' #kmax, zB, omk, H0_out, omdmh2'
          call EXIT(voidfeedback - voidfeedback)
          fatalerror = .true.
          return
       end if
       
    
!-----------------------------------------------------------!
! trick to follow the curvature: steep peak can be overseen:
!       dy(3)= voiddkdr(r)
       dy(3)= voidkofr(r)!*maxval(abs(dy((/1,2,4,5/))))*1.d3
!-----------------------------------------------------------!
       if(any(isnan(dy)).or.any(dy.gt.1.d30))then
          if(voidfeedback.gt.1)then
             write(*,'(A,3E14.5)')'r:',r
             write(*,'(A,3E14.5)')'y(1:3):',y
             write(*,'(A,3E14.5)')'dy(1:3):',dy
             xx=avoidS(r,t,u)
             write(*,*)'voidS:',xx
             xx=avoidSdotoverS(r,t,u)
             write(*,*)'voidSdotoverS:',xx
             xx=voidR(r,t,u)
             write(*,*)'voidR:',xx
             xx=avoidRp(r,t,u)
             write(*,*)'voidRp:',xx
          end if
          if(voidfeedback.gt.1)write(*,*)'error in ltb_aderivsdr.'
          !stop 'error.'
          VoidError=.true.
          !       return
       end if


!!$
       !    xx=avoidS(r,t,u)
       !    write(*,'(7E14.5)')dy,xx,y
!       write(*,*)thisS,thisSdot, vp%tdir
!!$    stop'dy'
       
    end if
    

    if(voiderror .and. .not. fatalerror)then
       dy=1.d300
       voiderror=.false.
    end if


    ! Error tracking:
!    DER%dy(1) = DER%S
!    DER%dy(2) = DER%Sdot * y(2)
 !   DER%dy(3) = 0.d0
 !   DER%dy(4) = 0.d0

end subroutine ltb_aderivsdr


function vdi_voiddtauda(a)
  implicit none
  real(dl) :: vdi_voiddtauda
  real(dl), intent(in) :: a
  
  real(dl) :: z
  real(dl), pointer :: tr(:,:)
  real(dl) :: Ypd,Yp,thisH
  
  if(voiderror)then
    vdi_voiddtauda = 1.d0
    return
  end if
  
  z = 1.d0/a - 1.d0
  
  allocate(tr(2,1))

  call MultiValuedInvert(VPNR,z,tr)
  if(Lbound(tr,1)>1)then
     call VoidErrorReport()
     stop 'It happened: lbound(tr:,1)>1. In vdi_voiddtauda(z).'
  end if

  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
    r=tr(2,1),t=tr(1,1),&
    Rp=Yp,Rpd=Ypd)
    
  deallocate(tr)
  
  thisH = Ypd/Yp
  
  vdi_voiddtauda = 1.d-3 * vpc / (thisH * a**2)
  
end function vdi_voiddtauda

function voidHL(z)
  implicit none
  real(dl) :: voidHL
  real(dl), intent(in) :: z
  
  real(dl), pointer :: tr(:,:)
  real(dl) :: Ypd,Yp,thisH
  
  if(voiderror)then
    voidHL = 1.d0
    return
  end if
  
  allocate(tr(2,1))

  call MultiValuedInvert(VPNR,z,tr)
  if(Lbound(tr,1)>1)then
     call VoidErrorReport()
     stop 'It happened: lbound(tr:,1)>1. In voidHL(z).'
  end if

  call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,&
    r=tr(2,1),t=tr(1,1),&
    Rp=Yp,Rpd=Ypd)
    
  deallocate(tr)
  
  thisH = Ypd/Yp
  
  voidHL = thisH
  
end function voidHL