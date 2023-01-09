

  ! Set your favorite profile here.
  subroutine VoidkofrContainer(k,dk,r)
    implicit none
    real(dl), optional :: k, dk
    real(dL), intent(in) :: r
    real(dl) :: rr
   real(dl) :: kofr, dkdr, rsign
!    real(dl), optional :: k, dk, r
!    real(dl) :: kofr, dkdr
    real(dl) :: a, b, kmax, kb, L!, maxr
    real(dl) :: x,y, mycosterm, xx
    integer :: powern

    rsign = r
    rr = dabs(r)
    


    ! Set variables:
    kmax = VP%kmax
    a = VP%alpha
    b = VP%beta
    L = VP%L
    kb = VP%kb_om

    if(VoidProfile == 1) then ! Our compensated profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          kofr = kmax * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - kmax * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if (VoidProfile == 2) then ! Our asymptotic profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VoidInfinity
       if(rr.lt.VP%LTB_maxr)then
          kofr = kmax /(1+(rr/L)**b) + kb
          dkdr = - kmax /(1+(rr/L)**b)**2 * b * rr**(b-1) / L**b
       else
          kofr = kb
          dkdr = 0.d0
       end if

    else if (VoidProfile == 3) then ! Our double void profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          kofr = VP%kmax_2 * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - VP%kmax_2 * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
          if(rr.lt.L*VP%L2_1)then
             kofr = kofr + (VP%kmax - VP%kmax_2)* (1.d0 - (rr/(L*VP%L2_1))**a)**b
             dkdr = dkdr - (VP%kmax - VP%kmax_2) * a * b * (1.d0 - (rr/(L*VP%L2_1))**a)**(b-1) * rr**(a-1) / (L*VP%L2_1)**a
          end if
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if (VoidProfile == 4) then ! Our sinvoid profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          kofr = kmax * (1.d0 - (dsin(3.d0*pi*rr/L) + 3.d0*pi*rr/L )/3.d0/pi) + kb
          dkdr = - kmax * ( 1.d0 / L * dcos(3.d0*pi*rr/L) + 1.d0 / L  )
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if (VoidProfile == 5) then ! Our sin4void profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          kofr = kmax * (1.d0 - (dsin(4.d0*pi*(rr/L-0.25d0)) + 4.d0*pi*rr/L )/4.d0/pi) + kb
          dkdr = - kmax * ( 1.d0 / L * dcos(4.d0*pi*(rr/L-0.25d0)) + 1.d0 / L  )
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if (VoidProfile == 6) then ! Our sin4void profile + extension

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          ! 0:(1-a)*L/2  (1-a)*L/2:1-(1-a)*L/2   1-(1-a)*L/2:L    
          if(rr<(1.d0-a)*L*0.5d0)then
             kofr = kmax * (1.d0 - (dsin(4.d0*pi*(rr/((1.d0-a)*L)-0.25d0)) + 4.d0*pi*rr/((1.d0-a)*L) )/4.d0/pi) + kb
             dkdr = - kmax * ( 1.d0 / ((1.d0-a)*L) * dcos(4.d0*pi*(rr/((1.d0-a)*L)-0.25d0)) + 1.d0 / ((1.d0-a)*L)  )
          else if (rr< (1.d0+a)*L*0.5d0) then
             kofr = 0.5d0*kmax + kb
             dkdr = 0.d0
          else
             kofr = kmax * (1.d0 - (dsin(4.d0*pi*((rr-a*L)/((1.d0-a)*L)-0.25d0)) + 4.d0*pi*(rr-a*L)/((1.d0-a)*L) )/4.d0/pi) + kb
             dkdr = - kmax * ( 1.d0 / ((1.d0-a)*L) * dcos(4.d0*pi*((rr-a*L)/((1.d0-a)*L)-0.25d0)) + 1.d0 / ((1.d0-a)*L)  )
          end if
       else
          kofr = kb
          dkdr = 0.d0
       end if
    else if (VoidProfile == 7) then ! Our EXP double void profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VoidInfinity
       if(rr.lt.VP%L)then
          kofr = VP%kmax_2 * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - VP%kmax_2 * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
       else
          kofr = kb
          dkdr = 0.d0
       end if

       kofr = kofr + (VP%kmax - VP%kmax_2)* dexp( - (rr/(L*VP%L2_1))**b )
       dkdr = dkdr + (VP%kmax - VP%kmax_2) * b * ( - (rr**(b-1)/(L*VP%L2_1)**(b)) )* dexp( - (rr/(L*VP%L2_1))**b )

    else if (VoidProfile == 8) then ! Our TANH void profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VoidInfinity
       if(rr.lt.VP%L)then
          kofr = VP%kmax_2 * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - VP%kmax_2 * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
!          if(r.lt.L*VP%L2_1)then
!          end if
       else
          kofr = kb
          dkdr = 0.d0
       end if

       kofr = kofr + (VP%kmax - VP%kmax_2)* (1 - dtanh( (rr/(L*VP%L2_1)  )) )
       dkdr = dkdr + (VP%kmax - VP%kmax_2) * ( - 1.d0 / dcosh(rr/(L*VP%L2_1))**2 )/(L*VP%L2_1)
          
    else if (VoidProfile == 9) then ! Our powerlaw void profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VoidInfinity
       if(rr.lt.VP%L)then
          kofr = VP%kmax_2 * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - VP%kmax_2 * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
!          if(r.lt.L*VP%L2_1)then
!          end if
       else
          kofr = kb
          dkdr = 0.d0
       end if
       
       kofr = kofr + (VP%kmax + VP%kmax_2) /(1+(rr/(L*VP%L2_1))**b) !(1.d0 - (rr/(L*VP%L2_1))**a)**b
       dkdr = dkdr - (VP%kmax + VP%kmax_2) * 2.d0 * rr /(L*VP%L2_1)**b /(1+(rr/(L*VP%L2_1))**b)**2 
          
    else if (VoidProfile == 10) then ! Our powerlaw void profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VoidInfinity
!       if(r.lt.VP%L)then
          kofr = (VP%kmax - VP%kmax_2) / (1.d0 + (rr/L)**2) + kb
          dkdr = - (VP%kmax - VP%kmax_2) * 2 / (1.d0 + (rr/L)**2)**2 * rr / L**2
!          if(r.lt.L*VP%L2_1)then
!          end if
!       else
!          kofr = kb
!          dkdr = 0.d0
!       end if
       
       kofr = kofr + ( VP%kmax_2) /(1+(rr/L))
       dkdr = dkdr - ( VP%kmax_2)  / L / (1+(rr/L))**2 
    else if (VoidProfile == 11) then ! Our overcompensated void profile without central lump

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VoidInfinity

       kofr = (VP%kmax) / (1.d0 + (rr/L)**2) + kb
       dkdr = - (VP%kmax) * 2 / (1.d0 + (rr/L)**2)**2 * rr / L**2
!!$       kofr = (VP%kmax - VP%kmax_2) / (1.d0 + (rr/L)**2) + kb
!!$       dkdr = - (VP%kmax - VP%kmax_2) * 2 / (1.d0 + (rr/L)**2)**2 * rr / L**2
       
       kofr = kofr + ( VP%kmax_2) /(1+(rr/L))**VP%L2_1 * (dtanh((rr - L )/(0.3*L)) - dtanh(-1.d0/3.d-1))
       dkdr = dkdr - ( VP%kmax_2)  / L *VP%L2_1 / (1+(rr/L))**(VP%L2_1+1)  * (dtanh((rr - L )/(0.3*L)) - dtanh(-1.d0/3.d-1))&
            + ( VP%kmax_2) /(1+(rr/L))**VP%L2_1 /  dcosh((rr - L )/(0.3*L))**2 / (0.3*L)
    else if(VoidProfile == 12) then ! Our compensated profile

       ! default in profile 1: 
       ! a=4
       ! b=2
       ! profile 12 = profile 1 with:
       a=2
       b=3
       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          kofr = kmax * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - kmax * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if(VoidProfile == 13) then ! Our compensated profile

       ! default in profile 1: 
       ! a=4
       ! b=2
       ! profile 12 = profile 1 with:
       a=8
       b=2
       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          kofr = kmax * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - kmax * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if(VoidProfile == 14) then ! Our compensated profile

       ! default in profile 1: 
       ! a=4
       ! b=2
       ! profile 12 = profile 1 with:
       a=20
       b=2
       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          kofr = kmax * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - kmax * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if (VoidProfile == 15) then ! Our double void profile

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          kofr = VP%kmax_2 * (1.d0 - (rr/L)**a)**b + kb
          dkdr = - VP%kmax_2 * a * b * (1.d0 - (rr/L)**a)**(b-1) * rr**(a-1) / L**a
          if(rr.lt.L*VP%L2_1)then
             kofr = kofr + (VP%kmax - VP%kmax_2)* (1.d0 - (rr/(L*VP%L2_1))**a)**b
             dkdr = dkdr - (VP%kmax - VP%kmax_2) * a * b * (1.d0 - (rr/(L*VP%L2_1))**a)**(b-1) * rr**(a-1) / (L*VP%L2_1)**a
          end if
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if (VoidProfile == 16) then ! Our cos void: continuous second derivative

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          if(rr.gt.L*VP%L2_1)then
             x=rr/L-VP%L2_1
             y=1._dl-VP%L2_1
             if(VP%d0_2<1._dl)stop 'd0_2 must be > 1 for profile 16.'
             if(VP%L2_1<0. .or. VP%L2_1>1.)stop 'L2_1 must be 0 < L2_1 < 1 for profile 16.'
             mycosterm = 1._dl-x/y + 0.5_dl/pi * Cos( pi * (2._dl*x/y - 0.5_dl)  )
             if(mycosterm<0.)mycosterm=0.d0 ! this only happens for errors O(10^-17).
             kofr = kb + VP%kmax * (mycosterm)**VP%d0_2
             dkdr = - VP%kmax / y / L * (1._dl +  Sin( pi * (2._dl*x/y - 0.5_dl)  ))
             if(VP%d0_2/=1._dl)dkdr = dkdr * VP%d0_2 * (mycosterm)**(VP%d0_2-1._dl)
          else             
             kofr = kb + VP%kmax
             dkdr = 0._dl
          end if
          if(isnan(kofr))then
             write(*,*)'k:',kofr
             write(*,*)'dkdr:',dkdr
             write(*,*)'mycosterm:',mycosterm
             write(*,*)'d0_2',VP%d0_2
             write(*,*)'mycosterm**d0_2',mycosterm**VP%d0_2
             write(*,*)'x,y:',x,y
             write(*,*)'x/y',x/y
             stop 'isnan(kofr)'
          end if
          
          
       else
          kofr = kb
          dkdr = 0.d0
       end if

    else if (VoidProfile == 17) then ! Our cos void: continuous second derivative

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          if(rr/L < 1._dl/3._dl)then
             kofr=kb
             dkdr =0._dl
          else if (rr/L < 2._dl/3._dl - VP%L2_1 / 3._dl)then
             x=-rr/L+2._dl/3._dl- VP%L2_1 / 3._dl
             y=1._dl/3._dl*(1._dl-VP%L2_1)
             if(VP%d0_2<1._dl)stop 'd0_2 must be > 1 for profile 17.'
             if(VP%L2_1<0. .or. VP%L2_1>1.)stop 'L2_1 must be 0 < L2_1 < 1 for profile 17.'
             mycosterm = 1._dl-x/y + 0.5_dl/pi * Cos( pi * (2._dl*x/y - 0.5_dl)  )
             if(mycosterm<0.)mycosterm=0.d0 ! this only happens for errors O(10^-17).
             kofr = kb + VP%kmax_2 * (mycosterm)**VP%d0_2
             dkdr = VP%kmax_2 / y / L * (1._dl +  Sin( pi * (2._dl*x/y - 0.5_dl)  ))
             if(VP%d0_2/=1._dl)dkdr = dkdr * VP%d0_2 * (mycosterm)**(VP%d0_2-1._dl)
             
          else if (rr/L > 2._dl/3._dl + VP%L2_1 / 3._dl) then
             x=rr/L-2._dl/3._dl- VP%L2_1 / 3._dl
             y=1._dl/3._dl*(1._dl-VP%L2_1)
             if(VP%d0_2<1._dl)stop 'd0_2 must be > 1 for profile 17.'
             if(VP%L2_1<0. .or. VP%L2_1>1.)stop 'L2_1 must be 0 < L2_1 < 1 for profile 17.'
             mycosterm = 1._dl-x/y + 0.5_dl/pi * Cos( pi * (2._dl*x/y - 0.5_dl)  )
             if(mycosterm<0.)mycosterm=0.d0 ! this only happens for errors O(10^-17).
             kofr = kb + VP%kmax_2 * (mycosterm)**VP%d0_2
             dkdr = - VP%kmax_2 / y / L * (1._dl +  Sin( pi * (2._dl*x/y - 0.5_dl)  ))
             if(VP%d0_2/=1._dl)dkdr = dkdr * VP%d0_2 * (mycosterm)**(VP%d0_2-1._dl)
             
          else
             kofr=kb + VP%kmax_2
             dkdr=0._dl

          end if
          
          if(isnan(kofr))then
             write(*,*)'k:',kofr
             write(*,*)'dkdr:',dkdr
             write(*,*)'mycosterm:',mycosterm
             write(*,*)'d0_2',VP%d0_2
             write(*,*)'mycosterm**d0_2',mycosterm**VP%d0_2
             write(*,*)'x,y:',x,y
             write(*,*)'x/y',x/y
             stop 'isnan(kofr)'
          end if
          
          
       else
          kofr = kb
          dkdr = 0.d0
       end if
          
    else if (VoidProfile == 18) then ! Our cos void: continuous third derivative

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L
       if(rr.lt.VP%LTB_maxr)then
          if(rr.gt.L*VP%L2_1)then
             x=rr/L-VP%L2_1
             y=1._dl-VP%L2_1
!             if(VP%d0_2<1._dl)stop 'd0_2 must be > 1 for profile 16.'
             if(VP%L2_1<0. .or. VP%L2_1>1.)stop 'L2_1 must be 0 < L2_1 < 1 for profile 18.'
             if(x/y<0.5_dl)then
                kofr = kb + VP%kmax*0.25_dl/pi**2 * ( &
                     1._dl + pi**2*(4._Dl-8._dl*(x/y)**2) - Cos(4._dl*pi*x/y))
                dkdr = VP%kmax / y / L *  ( &
                     -4._dl*x/y + Sin(4._dl*pi*x/y)/pi )
             else
                kofr = kb +VP%kmax*0.25_dl/pi**2 * ( &
                     -1._dl + 8._dl*pi**2*(-1._Dl+(x/y))**2 + Cos(4._dl*pi*x/y))
                dkdr = VP%kmax / y / L *  ( &
                     -4._dl + 4._dl*x/y - Sin(4._dl*pi*x/y)/pi )
             end if
          else             
             kofr = kb +VP%kmax
             dkdr = 0._dl
          end if

          if(isnan(kofr))then
             write(*,*)'k:',kofr
             write(*,*)'dkdr:',dkdr
             write(*,*)'mycosterm:',mycosterm
             write(*,*)'d0_2',VP%d0_2
             write(*,*)'mycosterm**d0_2',mycosterm**VP%d0_2
             write(*,*)'x,y:',x,y
             write(*,*)'x/y',x/y
             stop 'isnan(kofr)'
          end if
          
       else
          kofr = kb
          dkdr = 0.d0
       end if

    else if (VoidProfile == 19) then ! Our exp function, continuous up to n derivatives.

       if(.not.(associated(VP%LTB_maxr )))VP%LTB_maxr => VP%L

       if(rr.lt.VP%LTB_maxr)then
          if(rr.gt.L*VP%L2_1)then
          powern=3
             x=rr/L-VP%L2_1
             y=1._dl-VP%L2_1
             xx=x/y
             if(VP%L2_1<0. .or. VP%L2_1>1.)stop 'L2_1 must be 0 < L2_1 < 1 for profile 19.'
                kofr = kb + VP%kmax* ( 1.d0 - Exp(-(1-xx)**powern/xx)) 
                dkdr = - VP%kmax * (((powern*(1 - xx)**(-1 + powern))/(L*x) + (1 - xx)**powern/(L*x*xx))*Exp(-(1 - xx)**powern/xx))
          else             
             kofr = kb +VP%kmax
             dkdr = 0._dl
          end if

          if(isnan(kofr))then
             write(*,*)'k:',kofr
             write(*,*)'dkdr:',dkdr
             write(*,*)'mycosterm:',mycosterm
             write(*,*)'d0_2',VP%d0_2
             write(*,*)'mycosterm**d0_2',mycosterm**VP%d0_2
             write(*,*)'x,y:',x,y
             write(*,*)'x/y',x/y
             stop 'isnan(kofr)'
          end if
          

       else
          kofr = kb
          dkdr = 0.d0
       end if

    else

       
       stop 'VoidProfile not defined.'
    end if

    if(present(k))k=kofr
    if(present(dk))dk=dkdr
    
    rr = rsign

  end subroutine VoidkofrContainer


  function Voidkofr(r)
    implicit none
    real(dl) :: Voidkofr! , thisk
    real(dl), intent(in) :: r! , thisk

    if(VP%L>0.d0)then

       call VoidkofrContainer(k=Voidkofr, r=r)

    else
       voidkofr=VP%kb_om
    end if

  end function Voidkofr

  function Voiddkdr(r)
    implicit none
    real(dl) :: Voiddkdr! , thisdk
    real(dl), intent(in) :: r! , thisk
    

    if(VP%L>0.d0)then
       
       call VoidkofrContainer(dk=Voiddkdr, r=r)

    else
       voiddkdr=0.d0
    end if

       
    
  end function Voiddkdr
  
