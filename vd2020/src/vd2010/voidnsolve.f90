
  function voidu(r,t)
    implicit none
    real(dl) :: voidu(3) ! (real, imaginary, approximation or numerical)
    real(dl), intent(in) :: r,t
    real(dl) :: sinhu_u, sinq_q, y(2), thisk
    real(dl), parameter :: my_incr=1.1d0, my_precision=u_precision!1.d-7
    integer, parameter :: mymax = 10000
!    integer :: counter

    if(t.lt.0.d0)then
       !stop 'voidu(r,t) not written for negative t.'
!       if(.not.VoidError)write(*,*) 'voidu(r,t) not written for negative t.'
       VoidError=.true.
    end if
    if(t.eq.0.d0)then
       voidu=0.d0
       return
    end if

    if(r.eq.0.d0 .and. VP%L.gt.0.d0)then
       thisk= VP%kmax + VP%kb_om
    else if (r.eq.0.d0) then
       thisk = VP%kb_om
    else
       thisk =  voidkofr(r)
    end if

    if(thisk.gt.0.d0)then
       sinhu_u = t * 3.d0 * thisk**(1.5d0) * VP%Mtilde / dsqrt(2.d0)/pi

       y = uofsinhu_u(sinhu_u)
       voidu(1)=dabs(y(1)) !sign(y,t)
       voidu(2)=0.d0
       voidu(3)=y(2)

    else if(thisk.lt.0.d0) then
       sinq_q = t * (-3.d0)  * dabs(thisk)**(1.5d0) * VP%Mtilde  / dsqrt(2.d0) / pi

       y = qofsinq_q(sinq_q)

       voidu(1) = 0.d0
       voidu(2) = dabs(y(1)) ! sinq - q < 0 ---> q positive
       voidu(3)=y(2)
    else if (thisk.eq.0.d0) then
       voidu=-1.d0
    else
       write(*,*) 'voidkofr(r) is not a number in function voidu.'
       VoidError=.true.
       return
    end if

  end function voidu

  function uofsinhu_u(sinhu_uin)
    implicit none
    real(dl) :: sinhu_u, sinhu_uin
    real(dl) ::  uofsinhu_u(2)
    real(dl) :: r,t, f, df, dy
    real(dl) :: sinhy_y, y, ytop, ybot!, thisk
    real(dl), parameter :: my_incr=1.1d0, my_precision=u_precision!1.d-7
    real(dl), parameter :: maxy = asinh(max_double)
    integer, parameter :: mymax = 10000
    integer :: counter
    logical :: voiduerror
    if(voiderror)then
       uofsinhu_u = 1.d30
       return
    end if

    voiduerror = .false.

    sinhu_u = dabs(sinhu_uin)

    ! Analytical approximation if precision allows for it:
    if(dabs(6.d0*sinhu_u)**(2.d0/3.d0).lt.u_precision)then 
       ! 2/3, because u ~ sinhu_u^(1/3), error in sinhu_u is ~u^5, 
       ! relative error is u^5 / sinhu_u simeq sinhu_u^2/3
       sinhy_y = dabs(sinhu_u)
       y=(6.d0*sinhu_u)**(1.d0/3.d0)
       uofsinhu_u(1)=y
       uofsinhu_u(2)=-1.d0
       return
    end if
    uofsinhu_u(2)=1.d0



    ! first guess: u=t
    ! do algorithm for positive values always
    y=(6.d0*sinhu_u)**(1.d0/3.d0)
    if (y>maxy) y=maxy
    

    ytop=y
    ybot=0.d0
    sinhy_y = dsinh(y) - y
    if(sinhy_y.lt.sinhu_u)then
       ytop=ytop*2.d0
    end if
    counter =0
    do while (dabs(dabs(sinhy_y/sinhu_u)-1.d0).gt.my_precision*1.d3)
       sinhy_y = dsinh(y) - y
       if(sinhy_y.gt.dabs(sinhu_u))then
          ytop=y
          y=(y+ybot)*5.d-1
       else
          ybot=y
          if(y.ge.ytop)then
             y=y*my_incr
             ytop=y
          else
             y=(y+ytop)*0.5d0
          end if
       end if
       counter=counter+1
       if(counter.gt.mymax)then
          write(*,'(A,E17.8,A,E17.8)')'void u cannot find u(r,t), for r:',r,'t:',t
          write(*,'(4E17.8)')sinhy_y, sinhu_u, y, dabs(dabs(sinhy_y/sinhu_u)-1.d0)
          !stop 'voidu error.'
          voiduerror = .true.
          exit
!          VoidError=.true.
 !         return
       end if

    end do
    
    counter = 0 
!write(*,*)sinhu_uin/(dsinh(y)-y)-1.d0

    ! Improve accuracy with simple Newton-Raphson method from NR:
    do while (dabs(dabs(sinhy_y/sinhu_u)-1.d0).gt.0.d0)!my_precision/1.d2)
       f = dsinh(y) - y - sinhu_u
       df = dcosh(y) - 1.d0
       dy = f / df
       y = y - dy
!       write(*,'(5E14.5)')dy, f, df, y,sinhu_u/(dsinh(y)-y)-1.d0
       counter = counter + 1
       if(counter .ge. 20) exit
    end do
    sinhy_y = dsinh(y) - y
    if(voiduerror)then
       if( dabs(dabs(sinhy_y/sinhu_u)-1.d0).lt.1.d2*my_precision)then
          write(*,*)'Accepting:',dabs(dabs(sinhy_y/sinhu_u)-1.d0)
          voiduerror = .false.
       else
          voiderror = .true.
       end if
    end if
    !       if(isnan(sinhy_y)
    
    uofsinhu_u(1)=y
    !write(*,*)sinhu_uin/(dsinh(y)-y)-1.d0
    !stop
    
  end function uofsinhu_u

  function qofsinq_q(sinq_qin)
    implicit none
    real(dl) :: sinq_q, sinq_qin
    real(dl) :: qofsinq_q(2)
    real(dl) :: r,t, f, df, dy
    real(dl) :: siny_y, y, ytop, ybot!, thisk
    real(dl), parameter :: my_incr=1.1d0, my_precision=u_precision!1.d-4
    integer, parameter :: mymax = 10000
    integer :: counter
    logical :: voidqerror
    if(voiderror)then
       qofsinq_q = 1.d30
       return
    end if
    voidqerror = .false.

    sinq_q = - dabs(sinq_qin)

    ! Analytical approximation if precision allows for it:
    if(dabs(6.d0*sinq_q)**(2.d0/3.d0)/20.d0.lt.u_precision)then 
       ! 2/3, because u ~ sinhu_u^(1/3), error in sinhu_u is ~u^5, 
       ! relative error is u^5 / sinhu_u simeq sinhu_u^2/3
       y=(6.d0*(-sinq_q))**(1.d0/3.d0)
       qofsinq_q(1)=y
       qofsinq_q(2)=-1.d0
       return
    end if
    qofsinq_q(2)=1.d0



    ! first guess: u=t
    ! do algorithm for positive values always
    y=(-6.d0*sinq_q)**(1.d0/3.d0)

    ytop=y
    ybot=0.d0
    siny_y = dsin(y) - y
    if(siny_y.lt.sinq_q)then
       ytop=ytop*2.d0
    end if
    counter =0
    do while (dabs(dabs(siny_y/sinq_q)-1.d0).gt.my_precision*1.d3)
       siny_y = dsin(y) - y
       if(siny_y.lt.sinq_q)then
          ytop=y
          y=(y+ybot)*5.d-1
       else
          ybot=y
          if(y.ge.ytop)then
             y=y*my_incr
             ytop=y
          else
             y=(y+ytop)*0.5d0
          end if
       end if
       counter=counter+1
       if(counter.gt.mymax)then
          write(*,*)'qofsinq_q cannot find u(r,t), for r:',r,'t:',t
          write(*,*)siny_y, sinq_q, y
          !stop 'voidu error.'
          voidqerror = .true.
          exit
!          VoidError=.true.
       end if

    end do

    qofsinq_q(1)=y
    counter = 0 


    ! Improve accuracy with simple Newton-Raphson method from NR:
    do while (dabs(dabs(siny_y/sinq_q)-1.d0).gt.0.d0)!my_precision/1.d2)
       f = dsin(y) - y - sinq_q
       df = dcos(y) - 1.d0
       dy = f / df
       y = y - dy
!       write(*,'(5E14.5)')dy, f, df, y,sinq_q/(dsin(y)-y)-1.d0
       counter = counter + 1
       if(counter .ge. 10) exit
    end do
    
    siny_y = dsin(y) - y
    
    if(voidqerror)then
       if( dabs(dabs(siny_y/sinq_q)-1.d0).lt.1.d2*my_precision)then
          voidqerror = .false.
       else
          voiderror = .true.
       end if
    end if

    qofsinq_q(1)=y


  end function qofsinq_q

