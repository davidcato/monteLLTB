module wlltb_kofr
  use wlltb_constants
  use kofrattributes
  use wlltb_types
  implicit none

  public fk1, fk2
  
contains

  function fk1(r,att)
    implicit none
    real(dl), intent(in) :: r
    type(kofr_attributes), intent(in) :: att
    real(dl) :: fk1
    
    real(dl) :: L, kmax, a, b, kb
    
    L = att%L
    kmax = att%kmax
    a = att%a
    b = att%b
    kb = att%kb
    
    fk1 = kmax * (1 - (r/L)**a )**b + kb
  end function fk1

  function fk1prime(r,att)
    implicit none
    real(dl), intent(in) :: r
    type(kofr_attributes), intent(in) :: att
    real(dl) :: fk1prime
    
    real(dl) :: L, kmax, a, b, kb
    
    L = att%L
    kmax = att%kmax
    a = att%a
    b = att%b
    kb = att%kb
    
    fk1prime = b * a * kmax * (1 - (r/L)**a )**(b-1) * (r/L)**(a-1) + kb
  end function fk1prime

  function fk2(r,att)
    implicit none
    real(dl), intent(in) :: r
    type(kofr_attributes), intent(in) :: att
    real(dl) :: fk2
    
    real(dl) :: L, kmax, a, b, kb
    
    L = att%L
    kmax = att%kmax
    a = att%a
    b = att%b
    kb = att%kb
    
    fk2 = kmax * (1 - (r/L)**a )**b / 1.d3
  end function fk2

  subroutine fk3cont(r,att, kofr, dkdr)
    implicit none
    real(dl), intent(in) :: r
    real(dl), intent(out), optional :: kofr, dkdr
    type(kofr_attributes), intent(in) :: att
    real(dl) :: fk2
    
    real(dl) :: L, kmax, a, b, kb, L2_1, d0_2, mkofr, mdkdr, rr,x , y
    
    L = att%L
    kmax = att%kmax
    a = att%a
    b = att%b
    kb = att%kb
    L2_1 = 0.d0
    d0_2 = 0.d0
    rr = r

    if(rr.lt.L)then
      if(rr.gt.L*L2_1)then
        x=rr/L-L2_1
        y=1._dl-L2_1
!             if(VP%d0_2<1._dl)stop 'd0_2 must be > 1 for profile 16.'
        if(L2_1<0. .or. L2_1>1.)stop 'L2_1 must be 0 < L2_1 < 1 for profile 18.'
        if(x/y<0.5_dl)then
          mkofr = kb + kmax*0.25_dl/pi**2 * ( &
            1._dl + pi**2*(4._Dl-8._dl*(x/y)**2) - Cos(4._dl*pi*x/y))
          mdkdr = kmax / y / L *  ( &
            -4._dl*x/y + Sin(4._dl*pi*x/y)/pi )
        else
          mkofr = kb +kmax*0.25_dl/pi**2 * ( &
            -1._dl + 8._dl*pi**2*(-1._Dl+(x/y))**2 + Cos(4._dl*pi*x/y))
          mdkdr = kmax / y / L *  ( &
            -4._dl + 4._dl*x/y - Sin(4._dl*pi*x/y)/pi )
        end if
      else             
        mkofr = kb +kmax
        mdkdr = 0._dl
      end if

      if(isnan(mkofr))then
        write(*,*)'k:',mkofr
        write(*,*)'dkdr:',mdkdr
!        write(*,*)'mycosterm:',mycosterm
        write(*,*)'d0_2',d0_2
!        write(*,*)'mycosterm**d0_2',mycosterm**d0_2
        write(*,*)'x,y:',x,y
        write(*,*)'x/y',x/y
        stop 'isnan(kofr)'
      end if
          
          
    else
      mkofr = kb
      mdkdr = 0.d0
    end if

    if(present(kofr))kofr = mkofr
    if(present(dkdr))dkdr = mdkdr

  end subroutine fk3cont

  function fk3(r,att)
    implicit none
    real(dl), intent(in) :: r
    type(kofr_attributes), intent(in) :: att
    real(dl) :: fk3
    
    call fk3cont(r,att,kofr=fk3) 
  end function fk3

  function fk3prime(r,att)
    implicit none
    real(dl), intent(in) :: r
    type(kofr_attributes), intent(in) :: att
    real(dl) :: fk3prime
    
    call fk3cont(r,att,dkdr=fk3prime)    
  end function fk3prime
  

end module wlltb_kofr