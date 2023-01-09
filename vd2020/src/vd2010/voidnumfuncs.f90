  function nvoidS(r,t,u1)
    implicit none
    real(dl) nvoidS
    real(dl) :: r,t,u(3)
    real(dl) :: k, thisRp
    real(dl), optional :: u1(3)

    if(present(u1))then
       u=u1
    else
       u=voidu(r,t)
    end if

    k=voidkofr(r)

    thisRp= nvoidRp(r,t,u)

    nvoidS = thisRp / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )

    !    write(*,'(A,E14.5,A,E14.5,A,E14.5,A,E14.5)')'voidS:',voidS,' voidRp:',thisRp,' r:',r,' t:',t

    ! Error tracking:
!    DER%S = DER%Rp / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
  end function nvoidS

  function nvoidRpdot(r,t,u_opt)

!  function voidSdotoverS(r,t,u_opt)
    implicit none
!    real(dl) :: voidSdotoverS
    real(dl) :: nvoidRpdot
    real(dl) :: r,t,u(3)
    real(dl) :: aH1, aH, dr,aHp1,aHp2, r1, r2, u1(3), u2(3), Rpd, Rpd1, Rpd2, aH2!, Rp1, Rp2
    real(dl) :: distance
    real(dl), parameter :: myprecision = numerical_derivatives_precision !myprecision = 1.d-7
    integer, parameter :: mymax = 10000
    integer :: counter
    real(dl) :: sinhu_u, k, prev(3)

    real(dl), optional :: u_opt(3)
! Numderical recipes:
    real(dl), parameter :: con = 1.4, con2=con*con, big=1.d30, safe = 2.d0
    integer, parameter :: Ntab=10
    integer :: nri, nrj
    real(dl) :: err, errt, fac, nra(ntab,ntab), thisdfdx
! END NR

    if(present(u_opt))then
       u=u_opt
    else
       u=voidu(r,t)
    end if

!!$    ! Sdot / S = R'dot / R'.
!!$    ! a = R / r
!!$    ! R'dot = r a'dot + adot
!!$    !       = r (aH)' + aH
!!$    ! We know a and H analytically as function of t and r.
!!$
!!$    if((r.ge.VP%L .and. VP%kb_om.eq.0.d0) .or.  any(u(1:2).lt.0.d0))then
!!$       ! FLRW: Sdot / S = H(t)
!!$       voidSdotoverS=2.d0/3.d0/t
!!$       return
!!$    end if
!!$
!!$    ! Analytical approximation if voidu uses it too:
!!$!    sinhu_u = t * 3.d0 * voidkofr(r)**(1.5d0) * VP%Mtilde / dsqrt(2.d0)/pi
!!$!    if(dabs(6.d0*sinhu_u)**(2.d0/3.d0).lt.u_precision)then 
!!$    if(u(3).lt.0.d0)then 
!!$!       k=voidkofr(r)
!!$       voidSdotoverS = 2.d0 / 3.d0 / t       
!!$       !       voidSdot = voidRpdot(r,t)/ ( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
!!$       if(voidfeedback.gt.2)then
!!$          write(*,*)'Applying analytical approximation for voidSdotoverS for r close to L.'
!!$          write(*,'(A,5E14.5)')' voidSdotoverS,r,t,L:', voidSdotoverS,r,t,VP%L
!!$       end if
!!$       return
!!$    end if

    k = voidkofr(r)

! Numerical derivative here
    if(r.eq.0.d0)then
       dr=1.d-1 * VP%L
    else
       dr=r*1.d-1
    end if
    if(dr.eq.0.d0)dr=1.d-2

    aH=voidaH(r,t,u)! * voidH(r,t,u)
    
    prev(:)=1.d30

    distance=1.d0
    counter = 0

    r1=r+dr*0.5d0
    u1=voidu(r1,t)
    aH1=voidaH(r1,t,u1)! * voidH(r1,t,u1)
    !aHp1=(aH1 - aH)/dr*2.d0

    if(r-0.5d0*dr.gt.0.d0)then
       r2=r-0.5d0*dr
    else
       r2=r+dr
    end if
    u2=voidu(r2,t)
    aH2= voidaH(r2,t,u2)


    nra(1,1)=r * (aH1-aH2)/(r1-r2) + aH
    
    err = BIG

    do nri = 2,ntab
       dr = dr/con

       r1=r+dr*0.5d0
       u1=voidu(r1,t)
       aH1=voidaH(r1,t,u1)! * voidH(r1,t,u1)
       !aHp1=(aH1 - aH)/dr*2.d0
       
       if(r-0.5d0*dr.gt.0.d0)then
          r2=r-0.5d0*dr
       else
          r2=r+dr
       end if
       u2=voidu(r2,t)
       aH2= voidaH(r2,t,u2)
       
       
       nra(1,nri)=r * (aH1-aH2)/(r1-r2) + aH

       fac = con2
       do nrj=2,nri
          nra(nrj,nri)=(nra(nrj-1,nri)*fac-nra(nrj-1,nri-1))/(fac-1.d0)
          fac=con2*fac
          errt=max(abs(nra(nrj,nri)-nra(nrj-1,nri)),abs(nra(nrj,nri)-nra(nrj-1,nri-1)))
          if(errt.le.err)then
             err=errt
             thisdfdx=nra(nrj,nri)
          end if
       end do
       if(abs(nra(nri,nri)-nra(nri-1,nri-1)).ge.SAFE*err)exit
    end do

    Rpd1=thisdfdx

    Rpd = Rpd1
    nvoidRpdot = Rpd


    ! Error tracking:
!    DER%Rpdot = err
  end function nvoidRpdot

  function nvoidSdotoverS(r,t,u_opt)
    implicit none
    real(dl) :: nvoidSdotoverS
!    real(dl) :: voidRpdot
    real(dl) :: r,t,u(3)
    real(dl) :: aH1, aH, dr,aHp1,aHp2, r1, r2, u1(3), u2(3), Rpd, Rpd1, Rpd2, aH2!, Rp1, Rp2
    real(dl) :: distance
    real(dl), parameter :: myprecision = numerical_derivatives_precision !myprecision = 1.d-7
    integer, parameter :: mymax = 10000
    integer :: counter
    real(dl) :: sinhu_u, k, prev(3)

    real(dl), optional :: u_opt(3)
! Numderical recipes:
    real(dl), parameter :: con = 1.4, con2=con*con, big=1.d30, safe = 2.d0
    integer, parameter :: Ntab=10
    integer :: nri, nrj
    real(dl) :: err, errt, fac, nra(ntab,ntab), thisdfdx
! END NR

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
       nvoidSdotoverS=2.d0/3.d0/t
       return
    end if

    Rpd = nVoidRpdot(r,t,u)

    nvoidSdotoverS=Rpd / nvoidRp(r,t,u)
    !write(*,*)voidSdotoverS
    !stop'Sdot'
  end function nvoidSdotoverS

  function nvoidSdot(r,t,u_opt)
    implicit none
    real(dl) :: nvoidSdot
!    real(dl) :: voidRpdot
    real(dl) :: r,t,u(3)
    real(dl) :: aH1, aH, dr,aHp1,aHp2, r1, r2, u1(3), u2(3), Rpd, Rpd1, Rpd2, aH2!, Rp1, Rp2
    real(dl) :: distance
    real(dl), parameter :: myprecision = numerical_derivatives_precision !myprecision = 1.d-7
    integer, parameter :: mymax = 10000
    integer :: counter
    real(dl) :: sinhu_u, k, prev(3)

    real(dl), optional :: u_opt(3)
! Numderical recipes:
    real(dl), parameter :: con = 1.4, con2=con*con, big=1.d30, safe = 2.d0
    integer, parameter :: Ntab=10
    integer :: nri, nrj
    real(dl) :: err, errt, fac, nra(ntab,ntab), thisdfdx
! END NR

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
       nvoidSdot=nVoidSDotoverS(r,t,u) * nvoidS(r,t,u) !2.d0/3.d0/t
       return
    end if

    ! And the next line fucks every thing down (i.e. compensating fucking up):
    k = voidkofr(r)
    ! do not forget next time!

    Rpd = nVoidRpdot(r,t,u)

    nvoidSdot=Rpd / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
    !write(*,*)voidSdotoverS
    !stop'Sdot'

    ! Error tracking:
!    DER%Sdot = DER%Rpdot / dsqrt( 1.d0 + 2.d0 * r**2 * k * VP%Mtilde**2  )
  end function nvoidSdot


  function nvoidRp(r,t,u_opt)
    implicit none
    real(dl) :: nvoidRp
    real(dl) :: r,t,u(3)
    real(dl) :: RR, RR1, RR2, dr,Rp1,Rp2, r1, r2, u1(3), u2(3)
    real(dl) :: distance
    real(dl), parameter :: mymindr=1.d-14
    real(dl), parameter :: myprecision = numerical_derivatives_precision !myprecision = 1.d-7
    integer, parameter :: mymax = 10000
    integer :: counter
    real(dl) :: sinhu_u, prev(3)
    real(dl), optional :: u_opt(3)
! Numderical recipes:
    real(dl), parameter :: con = 1.4, con2=con*con, big=1.d30, safe = 2.d0
    integer, parameter :: Ntab=10
    integer :: nri, nrj
    real(dl) :: err, errt, fac, nra(ntab,ntab), thisdfdx
! END NR

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
       nvoidRp = VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
       !   write(*,*)voidRp,r,t
       !   stop 'r=0'

       if(voidfeedback.gt.2)then
          write(*,*)'Applying analytical approximation for nvoidRp.'
          write(*,'(A,5E14.5)')' nvoidRp,r,t,sinhu_u, L:', nvoidRp,r,t,sinhu_u, VP%L
       end if
       return
    end if

! Numerical derivative here
    if(r.eq.0.d0)then
       dr=1.d-1 * VP%L
    else
       dr=r*1.d-1
       if(r+dr.gt.VP%L)dr=(VP%L-r)*0.5d0
    end if
    if(dr.eq.0.d0)dr=1.d-2

    r1=r+dr*0.5d0
    u1=voidu(r1,t)
    RR1= voidR(r1,t,u1)

    if(r-0.5d0*dr.gt.0.d0)then
       r2=r-dr*0.5d0
    else
       r2=r+dr
    end if
    u2=voidu(r2,t)
    RR2 = voidR(r2,t,u2)
    
    prev(:)=1.d30

    distance=1.d0
    counter = 0



    nra(1,1)= (RR2-RR1)/(r2-r1)
    
    err = BIG

    do nri = 2,ntab
       dr = dr/con
!       if(r+dr.gt.VP%L)dr=(VP%L-r)*0.5d0
       r1=r+dr*0.5d0
       u1=voidu(r1,t)
       RR1= voidR(r1,t,u1)

       if(r-0.5d0*dr.gt.0.d0)then
          r2=r-dr*0.5d0
       else
          r2=r+dr
       end if
       u2=voidu(r2,t)
       RR2 = voidR(r2,t,u2)
  
       nra(1,nri)= (RR2-RR1)/(r2-r1)

       fac = con2
       do nrj=2,nri
          nra(nrj,nri)=(nra(nrj-1,nri)*fac-nra(nrj-1,nri-1))/(fac-1.d0)
          fac=con2*fac
          errt=max(abs(nra(nrj,nri)-nra(nrj-1,nri)),abs(nra(nrj,nri)-nra(nrj-1,nri-1)))
          if(errt.le.err)then
             err=errt
             thisdfdx=nra(nrj,nri)
          end if
       end do
       if(abs(nra(nri,nri)-nra(nri-1,nri-1)).ge.SAFE*err)exit
    end do
    Rp1=thisdfdx
!!write(*,*)Rp1

!--
!!$
!!$
!!$    if(r.eq.0.d0)then
!!$       dr=1.d-5
!!$    else
!!$       dr=r*5.d-2
!!$       if(r+dr.gt.VP%L)dr=(VP%L-r)*0.5d0
!!$    end if
!!$    u = voidu(r,t)
!!$    RR=voidR(r,t,u)
!!$    
!!$    prev(:)=1.d30
!!$
!!$    distance=1.d0
!!$    counter = 0
!!$    do while (distance.gt.myprecision)
!!$       r1=r+dr*0.5d0
!!$       u1=voidu(r1,t)
!!$       RR1= voidR(r1,t,u1)
!!$       Rp1=(RR1 - RR)/dr*2.d0
!!$
!!$       if(r-0.5d0*dr.gt.0.d0)then
!!$          r2=r-dr*0.5d0
!!$       else
!!$          r2=r+dr
!!$       end if
!!$       u2=voidu(r2,t)
!!$       Rp2=(voidR(r2,t,u2) - RR)/(r2-r)
!!$
!!$       dr=dr*0.9d0
!!$
!!$       distance = dabs(Rp1/Rp2-1.d0)
!!$       if(isnan(distance))distance = dabs(Rp2/Rp1-1.d0)
!!$       if(isnan(distance))distance = dabs(Rp2-Rp1)
!!$
!!$!if(abs(r/ 4.194358881034884E-003-1.d0).lt.1.d-5)then
!!$!   write(*,'(4E16.7)')Rp1,Rp2, RR1, RR
!!$!end if
!!$
!!$       if(voidfeedback.gt.3.and.(isnan(Rp1).or.isnan(Rp2)))then
!!$          write(*,'(A,4E14.5)')' r,t,u,L:',r,t,u,VP%L
!!$          write(*,*)'voidR(r,t,u):',RR
!!$          write(*,*)'voidR(r1,t,u):',RR1
!!$          write(*,*)'Rp1:',Rp1
!!$          write(*,*)'Rp2:',Rp2
!!$          write(*,*)'dr:',dr
!!$          write(*,*)'error in voidRp.'
!!$!          stop 'error in voidRp.'
!!$       end if
!!$       if(dr.lt.mymindr)exit
!!$       counter=counter+1
!!$       if(counter.gt.mymax)then
!!$          write(*,'(A,E16.7,A,E16.7,A,E16.7)')' voidRp cannot dR/dr, for r:',r,'t:',t
!!$          !stop 'voidRp error.'
!!$          VoidError=.true.
!!$          return
!!$       end if
!!$       if(prev(3)/distance-1.d0.lt.0.d0)then
!!$          if(voidfeedback.gt.2)then
!!$             write(*,*)'Accepting low precision in voidRp, because precision decreases for even smaller dr.'
!!$             write(*,*)'present precision (distance):',distance
!!$             write(*,*)'previous precision, accepting:',prev(3)
!!$             write(*,*)'r',r
!!$             write(*,*)'t',t
!!$          end if
!!$          Rp1=prev(1)
!!$          Rp2=prev(2)
!!$          if(prev(3).gt.1.d-1)then
!!$!             VoidError=.true.
!!$!             return
!!$             Rp1 = 0.d0
!!$             Rp2 = 0.d0 ! This will lead to exclusion of this timestep.
!!$          end if
!!$          exit
!!$       end if
!!$
!!$       prev(1)=Rp1
!!$       prev(2)=Rp2
!!$       prev(3)=distance
!!$!write(*,*)Rp1, Rp2
!!$    end do
!!$
!!$!write(*,*)Rp1
!!$!stop
    nvoidRp=Rp1
    if(nvoidRp.eq.0.d0)then
       ! Analytical approximation if possible
       if(dabs(r/VP%L-1.d0).lt.1.d-1)then 
          nvoidRp = VP%Mtilde**(2.d0/3.d0) * (6.d0 * pi)**(1.d0/3.d0) * dabs(t)**(2.d0/3.d0)
          !   write(*,*)voidRp,r,t
          !   stop 'r=0'
          
          if(voidfeedback.gt.2)then
             write(*,*)'Applying analytical approximation for voidRp for r close to L.'
             write(*,'(A,5E14.5)')' voidRp,r,t,sinhu_u,L:', nvoidRp,r,t,sinhu_u, VP%L
          end if
          return
       end if
       if(voidfeedback.gt.2)then
          write(*,'(A,4E14.5)')'r,t,u,L:',r,t,u,VP%L
          write(*,*)'voidR(r,t,u):',RR
          write(*,*)'Rp1:',Rp1
          write(*,*)'voidR(r1,t,u):',RR1
          write(*,*)'Rp2:',Rp2
          write(*,*)'dr:',dr
       end if
!!$       write(*,*)'Error in voidRp: inexplicable voidRp=0.'
!!$       VoidError=.true.
!!$       return
       !      stop 'error in voidRp.'
    end if
    !    if(voidfeedback.gt.1)then
    !       write(*,*)'voidRp:', voidRp
    !    end if
    
    
    ! Error tracking:
!    DER%Rp = err
  end function nvoidRp



  ! y(1) = t
  ! y(2) = dt/du = - (z + 1)
  ! y(3) = k
  ! du = - dz
  subroutine ltb_nderivsdr(v,y,dy)
    implicit none
    real(dl), intent(in) :: v
    real(dl), intent(in) :: y(:)
    real(dl), intent(out) :: dy(:)
    real(dl) :: u(3),t,r, xx
    if(size(y).ne.3)stop 'Wrong dimension in ltb_derivs.'


    r=v!y(3)
!!$    if(r.gt.VP%L)then
!!$       dy=1.d30
!!$       write(*,*)'Forcing smaller stepsize, r > L.',v
!!$       return
!!$    end if
    t=y(1)
    u = voidu(r,t)


! r:

    dy(1) = - nvoidS(r,t,u)
!    dy(2) = voidSdotoverS(r,t,u) * voidS(r,t,u) * y(2)
    dy(2) = nvoidSdot(r,t,u) * y(2)
!    dy(3) = 1.d0

! z:
!!$    dy(1)= 1.d0 / (voidSdotoverS(r,t,u) * y(2))
!!$
!!$    dy(3)= - dy(1) / voidS(r,t,u)
!!$    dy(2)= - 1.d0

! affine:
!!$    dy(1)= y(2)
!!$
!!$    dy(3)= - y(2)/ voidS(r,t,u)
!!$    dy(2)= - voidSdotoverS(r,t,u) * dy(1)**2

!-----------------------------------------------------------!
! trick to follow the curvature: steep peak can be overseen:
!    dy(4)=voiddkdr(r)*dy(3)
    dy(3)=voiddkdr(r)
!-----------------------------------------------------------!

    if(any(isnan(dy)).or.any(dy.gt.1.d30))then
       if(voidfeedback.gt.1)then
          write(*,'(A,3E14.5)')'r:',r
          write(*,'(A,3E14.5)')'y(1:3):',y
          write(*,'(A,3E14.5)')'dy(1:3):',dy
          xx=nvoidS(r,t,u)
          write(*,*)'voidS:',xx
          xx=nvoidSdotoverS(r,t,u)
          write(*,*)'voidSdotoverS:',xx
          xx=voidR(r,t,u)
          write(*,*)'voidR:',xx
          xx=nvoidRp(r,t,u)
          write(*,*)'voidRp:',xx
       end if
       if(voidfeedback.gt.1)write(*,*)'error in ltb_derivs.'
       !stop 'error.'
!       VoidError=.true.
!       return
    end if
!!$
!!$    xx=voidS(r,t,u)
!    write(*,'(7E14.5)')dy,xx,y
!!$    stop'dy'
if(voiderror)then
   dy=1.d30
   voiderror=.false.
end if

    ! Error tracking:
!    DER%dy(1) = DER%S
!    DER%dy(2) = DER%Sdot * y(2)
 !   DER%dy(3) = 0.d0
 !   DER%dy(4) = 0.d0

end subroutine ltb_nderivsdr

