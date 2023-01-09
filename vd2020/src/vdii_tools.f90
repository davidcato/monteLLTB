!*
! *  vdii_tools.f90
! *  VoidDistancesII
! *
! *  Created by Wessel Valkenburg on 29/01/2012.
! *  Copyright 2012 Wessel Valkenburg. All rights reserved.
! *
! */

! 
module vdii_tools
  use vdii_constants
  implicit none

  private

  public vdii_D1D2

contains

  ! vdii_D1D2 returns the first and second derivative of f(x) in x,
  ! accurate up to order h^4.
  subroutine vdii_D1D2(g,x,d, h)
    implicit none
    real(dl), intent(in) :: x, h
    real(dl), intent(out) :: d(:)
    interface
      function g(xx)
        use vdii_constants
        implicit none
        real(dl) :: g
        real(dl), intent(in) :: xx
      end function g
    end interface
    
    real(dl), target :: fa(5)
    real(dl), pointer :: f, fph, fp2h, fmh, fm2h
    integer, parameter :: nmem=20
    real(dl), parameter :: myfac=2d0
    real(dl) :: hh, phh(nmem), a(nmem,nmem,2), erro(2), nerro(2), thisfac
    real(dl) :: myres(2), pres(nmem,size(myres)), improvemem(nmem,size(myres))
    integer :: count, ii, i, j
    if(size(d)<2)stop 'result array in vdii_D1D2 must have dimension >1'
    
    fm2h => fa(1)
    fmh => fa(2)
    f => fa(3)
    fph => fa(4)
    fp2h => fa(5)

    count = 0
!    if(present(h))then
      hh=h
!    else
!      hh=1.d-16 !* f
!      do while ((x+hh==x) .or. ( x+hh/=x .and. (g(x)==g(x+hh))))
!        hh=hh*1.d1
!        if(hh>1.d200)then
!          write(0,*)'Can find good starting h in vdii_d1d2.'
!            d=0.d0
!            return
!          stop
!        end if
!      end do
!    end if
if(.false.) then
    pres=0.d0
    myres=0.d0
    improvemem=0.d0
    do ii = 1, nmem
      pres(ii,:)=thisres()
      if(ii>1)then
        improvemem(ii,:)=abs(pres(ii,:)-pres(ii-1,:))
        where(abs(2*hh*pres(ii,:))>1.d-15 .and. improvemem(ii,:)/=0.d0)myres=pres(ii,:)
        !write(*,*)improvemem(ii,:),hh,myres,thisres()
      end if
      if(all(abs(2*hh*pres(ii,:))<1.d-15))exit
      hh=hh/myfac
    end do
    if(.not. all(myres==0.d0))then
      hh = (minval(abs(myres),abs(myres)/=0.d0) * 1.d-15)**0.25d0
      if(x+hh/=x)myres = thisres()
    end if

   
       
!    do while (all(x+hh/=x+phh))
!      write(0,*)hh
!      call pushphh(hh,phh)
!      if(.not. all(myres==0.d0))then
!        hh = (minval(abs(myres),abs(myres)/=0.d0) * 1.d-15)**0.25d0
!        if(x+hh/=x)myres = thisres()
!      end if
!    end do
else
    a=0.d0
    a(1,1,:) = thisres()
    erro=1.d30
    do i = 2, nmem
      hh=hh/myfac
      a(1,i,:)=thisres()
      thisfac = myfac**2
      do j = 2, i
        a(j,i,:)=(a(j-1,i,:)*thisfac - a(j-1,i-1,:))/(thisfac-1.d0)
        thisfac=thisfac*myfac**2
        nerro=abs(a(j,i,:)-a(j-1,i-1,:))
        where (abs(a(j,i,:)-a(j-1,i,:))>abs(a(j,i,:)-a(j-1,i-1,:)))
          nerro=abs(a(j,i,:)-a(j-1,i,:))
!        else
        end where
        where (nerro<erro)
          erro=nerro
          myres=a(j,i,:)
        end where
      end do
      if(all(a(i,i,:)-a(i-1,i-1,:)>2*erro))exit
    end do
!    write(*,*)'rel:',hh/myres, fa
end if    
    d(1:2)=myres

    nullify(f, fph, fp2h, fmh, fm2h)
    !write(*,*)'count:',count
  contains
    
    function thisres()
      implicit none
      real(dl) :: thisres(2)
      integer :: mi
  
      do mi = 1, 5
        fa(mi) = g(x+(mi-3)*hh)
        !write(*,*)'fa(I):',fa(i),x+(i-3)*hh
      end do
      
      ! first derivative:
      thisres(1) = (- fp2h + 8*fph - 8*fmh + fm2h) / (12.d0 * hh)
      
      ! second derivative:
      thisres(2) = (- fp2h + 16*fph - 30*f + 16*fmh - fm2h) / (12.d0 * hh**2)
      
      count = count + 1
    end function thisres
  
    subroutine pushphh(hh1,phh1)
      implicit none
      real(dl) :: hh1, phh1(:)
      integer :: i, n
      n=size(phh1)
      do i = n,2,-1
        phh1(i)=phh1(i-1)
      end do
      phh(1)=hh1
      
    end subroutine pushphh
  end subroutine vdii_D1D2


end module vdii_tools
