! WV's owm arbitrary precision multiplication and addition. Double precisions only.
! All internal operations use integer arrays, where each index corresponds to a decimal position.
! So arrays may range from -300 to 300, corresponding to the range of doubles.

module lltb_apprec
use lltb_precision

private

interface APQ
#ifdef NOQUAD
    module procedure APQ_d
#else
    module procedure APQ_d, APQ_q
#endif
end interface

public APQ, TestAP

contains

  subroutine APQ_d(Qout,xyz)
    implicit none
    real(dl), intent(in) :: xyz(3)
    real(ql) :: Qout(2)
    integer(dl), pointer :: xi(:), yi(:), zi(:)
    nullify(xi,yi,zi)

    call APDOubleToInteger(xyz(1),xi)
    call APDOubleToInteger(xyz(2),yi)
    call APDOubleToInteger(xyz(3),zi)

    call APQ_i(Qout,xi,yi,zi)

    deallocate(xi,yi,zi)
    nullify(xi,yi,zi)

  end subroutine APQ_d

  subroutine APQ_q(Qout,xyz)
    implicit none
    real(ql), intent(in) :: xyz(3)
    real(ql) :: Qout(2)
    integer(dl), pointer :: xi(:), yi(:), zi(:)
    nullify(xi,yi,zi)


    call APQuadToInteger(xyz(1),xi)
    call APQuadToInteger(xyz(2),yi)
    call APQuadToInteger(xyz(3),zi)

    call APQ_i(Qout,xi,yi,zi)

    deallocate(xi,yi,zi)
    nullify(xi,yi,zi)

  end subroutine APQ_q

  
  subroutine APQ_i(Qout,xi,yi,zi)
    implicit none
!    real(dl), intent(in) :: xyz(3)
    integer(dl), intent(in), pointer ::  xi(:),yi(:),zi(:)
    real(ql) :: Qout(2)
    integer(dl), pointer :: xi2(:),zi4(:),yi3(:),zi3(:),zi2(:),n81(:),n12(:), n9(:),t1(:),t2(:),t3(:),t1a(:), t2a(:), t3a(:), t12(:)
    character(len=4096) :: numstr

    nullify(xi2,zi4,yi3,zi3,zi2,n81,n12,n9,t1,t2,t3, t1a, t2a, t3a, t12)
    !write(*,'(3ES24.15)')xyz

!    call APDOubleToInteger(81.d0,n81)
!    call APDOubleToInteger(12.d0,n12)
!    call APDOubleToInteger(9.d0,n9)
    allocate(n81(0:1),n12(0:1),n9(0:0))
    n81=(/1,8/)
    n12=(/2,1/)
    n9=(/9/)

    call APPow_i(zi,2_dl,zi2)
    call APMult_i(zi,zi2,zi3)
    call APPow_i(zi2,2_dl,zi4)

    call APPow_i(yi,3_dl,yi3)

!!$    call APwrite(yi3,'y^3:')

    call APPow_i(xi,2_dl,xi2)

    call APMult_i(n81,zi4,t1a)
    call APMult_i(t1a,xi2,t1)
    call APMult_i(n12,zi3,t2a)
    call APMult_i(t2a,yi3,t2)

    call APSum_i(t1,t2,t12)

!    call APWrite(t12,'Sum:')

    call APMult_i(n9,zi2,t3a)
    call APMult_i(t3a,xi,t3)

!    call APWrite(t3,'Last term:')
    
    numstr = APIntegerToString(t12)
    read(numstr,*)Qout(1)
    numstr = APIntegerToString(t3)
    read(numstr,*)Qout(2)
    deallocate(xi2,zi4,yi3,zi3,zi2,n81,n12,n9,t1,t2,t3, t1a, t2a, t3a, t12)

  end subroutine APQ_i


  subroutine APWrite(iarr, str)
    integer(dl), pointer :: iarr(:)
    character(len=*), optional :: str

    if(present(str))then
       write(*,'(A,A)')trim(str),trim(APIntegerToString(iarr))
    else
       write(*,'(A)')trim(APIntegerToString(iarr))
    end if
  end subroutine APWrite

  function APIntegerToString(iarr)
    implicit none
    integer(dl), pointer :: iarr(:)
    character(len=4096) :: APIntegerToString
    integer :: mag, n, blow, sign
    character(len=4) :: ns
    character(len=1) :: myminus
    character(len=32) :: tform
    
    sign=1
    myminus=""
    if(all(iarr<=0))then
       sign=-1
       myminus="-"
    end if

    mag=ubound(iarr,dim=1)
    blow=lbound(iarr,dim=1)
    n=size(iarr)
    write(ns,'(I4)')n
    
    tform='("'//myminus//'",I1,".",'//trim(adjustl(ns))//'I1)'
    write(APIntegertoString,tform)sign*iarr(mag:blow:-1)
    
    write(ns,'(I4)')mag
    APIntegerToString=trim(adjustl(APIntegerToString))//"E"//trim(adjustl(ns))
    
    return

  end function APIntegerToString



  subroutine APDoubleToInteger(x,iarr)
    implicit none
    real(dl) :: x
    integer(dl), pointer :: iarr(:)
    character(len=24) :: xc, xc1, xc2
    integer, pointer :: alld(:,:)
    integer :: offset, i, j, k, blow, bup
    

    if(associated(iarr))then
       stop 'APDoubleToInteger must be called with an empty / unallocated / unassociated pointer.'
    end if

    write(xc,'(E26.17E3)')x
    read(xc(21:24),*)offset
    allocate(alld(17,2))
    do j = 1, 17
       k = j
       read( xc(2+j:2+j) , *) alld(k,1)
       alld(k,2) = offset - j 
       ! since first digit is x in 0.x
    end do

    blow=minval(alld(:,2))
    bup=maxval(alld(:,2))

    allocate(iarr(blow:bup))
    iarr=0
    
    do j = 1, 17
       iarr(alld(j,2))=alld(j,1)
    end do

    if(x<0)then
       iarr=-iarr
    end if

    deallocate(alld)
    
    
  end subroutine APDoubleToInteger



  subroutine APQuadToInteger(x,iarr)
    implicit none
    real(ql) :: x
    integer(dl), pointer :: iarr(:)
    character(len=41) :: xc, xc1, xc2
    integer, pointer :: alld(:,:)
    integer :: offset, i, j, k, blow, bup
    

    if(associated(iarr))then
       stop 'APDoubleToInteger must be called with an empty / unallocated / unassociated pointer.'
    end if
    
    write(xc,'(E43.33E4)')x
    read(xc(37:41),*)offset
    allocate(alld(33,2))
    do j = 1, 33
       k = j
       read( xc(2+j:2+j) , *) alld(k,1)
       alld(k,2) = offset - j 
       ! since first digit is x in 0.x
    end do

    blow=minval(alld(:,2))
    bup=maxval(alld(:,2))

    allocate(iarr(blow:bup))
    iarr=0
    
    do j = 1, 33
       iarr(alld(j,2))=alld(j,1)
    end do

    if(x<0)then
       iarr=-iarr
    end if

    deallocate(alld)
    


  end subroutine APQuadToInteger

  ! When subtracting, 
  ! and result must be positive.
  subroutine APAdd_i(iarr1, iarr2, ires)
    implicit none
    integer(dl), pointer :: iarr1(:), iarr2(:)
    integer(dl), pointer :: ires(:)
    integer(dl), pointer :: alld(:,:), sortd(:,:), mires(:)
    integer :: i, j, n, k, offset, blow,bup, finup, blow1, bup1, blow2, bup2
    character(len=23) :: xc, xc1, xc2

    if(associated(ires))then
       stop 'APSum_i must be called with an empty / unallocated / unassociated pointer for the result.'
    end if

    blow1=lbound(iarr1,dim=1)
    bup1=ubound(iarr1,dim=1)
    blow2=lbound(iarr2,dim=1)
    bup2=ubound(iarr2,dim=1)
    blow=min(blow1,blow2)
    bup=max(bup1,bup2)+2
    n=size(iarr1)+size(iarr2)

    finup=bup

    allocate(sortd(blow:finup,n))
    allocate(mires(blow:finup))
    sortd=0
    mires=0

    do i = blow1,bup1
       do k = 1, n
          if(sortd(i,k)==0)exit
       end do
       sortd(i,k)=iarr1(i)
    end do
    do i = blow2,bup2
       do k = 1, n
          if(sortd(i,k)==0)exit
       end do
       sortd(i,k)=iarr2(i)
    end do
   
    ! Here it goes: sum it all.
    ! Hoping there is still enough space in sortd for new digits
    do i = blow, finup
       do j = 1, n
          if(sortd(i,j)/=0)then
             mires(i)=mires(i)+sortd(i,j)
             ! intermediate overflowing:
             if(mires(i)>89)then
write(*,*)'taking off:',mires(i)
                if(i+1>finup)stop 'estimate of order of magnitude wrong in APSum'
                mires(i) = mires(i)-90
                do k = 1, size(sortd(1,:))!n
                   if(sortd(i+1,k)==0)exit
                end do
                if(sortd(i+1,k)/=0)stop 'sortd too small in APSum'
                sortd(i+1,k)=9
write(*,*)'newval:',mires(i)
             end if
          else
             exit
          end if
       end do
       ! final overflowing:
       if(mires(i)>9)then
          if(i+1>finup)stop 'estimate of order of magnitude wrong in APSum'
          do k = 1,size(sortd(1,:))! n
             if(sortd(i+1,k)==0)exit
          end do
          if(sortd(i+1,k)/=0)stop 'sortd too small in APSum'
          sortd(i+1,k)=int(real(mires(i))/10.)
          mires(i)=mires(i)-10*sortd(i+1,k)
       end if
       
    end do

    ! now mires is a list of digits with there exponents, 
    ! of the final result.

    ! determine size:
    do i = finup, blow, -1
       if(mires(i)/=0)exit
    end do
    n=i-blow+1
    finup=i ! final final upper limit
    ! size = n digits
    ! How many doubles?
    allocate(ires(blow:finup))
    ires=mires(blow:finup)

    deallocate(sortd, mires)

    if(any(ires<0))call APFinalizeSubtraction(ires)

  end subroutine APAdd_i

  subroutine APFinalizeSubtraction(ires)
    implicit none
    integer(dl), intent(inout), pointer :: ires(:)
    integer(dl), pointer :: sortd(:)
    integer :: blow,bup, flow, fup, n, i, j, k

    blow=lbound(ires,dim=1)
    bup=ubound(ires,dim=1)
    flow=blow-2
    fup=bup
    n=size(ires)

    allocate(sortd(flow:fup))
    sortd=0
    sortd(blow:bup)=ires

    do j = 1, n
       do i = flow, fup
          if(sortd(i)>=0)cycle
          sortd(i)=sortd(i)+10
          sortd(i+1)=sortd(i+1)-1
       end do
       if(.not.(any(sortd<0)))exit
    end do
    if(any(sortd<0))stop 'Cannot finish this subtraction. Order right?'

    ! determine size:
    do i = fup, flow, -1
       if(sortd(i)/=0)exit
    end do
    fup=i ! final final upper limit
    do i = flow, fup
       if(sortd(i)/=0)exit
    end do
    flow=i

    deallocate(ires)
    allocate(ires(flow:fup))
    ires=sortd(flow:fup)
    deallocate(sortd)

  end subroutine APFinalizeSubtraction
  

  subroutine APSum_i(xi,yi,ires)
    implicit none
    integer(dl), intent(in), pointer :: xi(:), yi(:)
    integer(dl), intent(inout), pointer :: ires(:)
    integer(dl), pointer :: biggest(:), smallest(:)
    integer :: sign, xsign, ysign, xlow,xup,ylow,yup,i, bsign

    xsign=1
    ysign=1
    if(any(xi<0))xsign=-1
    if(any(yi<0))ysign=-1
    sign = xsign*ysign

    xlow=lbound(xi,dim=1)
    xup=ubound(xi,dim=1)
    ylow=lbound(yi,dim=1)
    yup=ubound(yi,dim=1)

    if(sign==1)then ! means addition
       if((xsign==1).and.(ysign==1))then
          call APAdd_i(xi,yi,ires)
          return
       end if
       allocate(biggest(ylow:yup),smallest(xlow:xup))
       bsign=ysign
       biggest=ysign*yi
       smallest=xsign*xi
    else ! means subtraction
       ! find which is absolute biggest
       if(xup/=yup)then ! simple
          if(xup>yup)then
             allocate(biggest(xlow:xup),smallest(ylow:yup))
             biggest=xsign*xi
             smallest=-ysign*yi
             bsign=xsign
          else
             allocate(biggest(ylow:yup),smallest(xlow:xup))
             biggest=ysign*yi
             smallest=-xsign*xi
             bsign=ysign
          end if
       else
          do i = xup, max(ylow,xlow), -1
             if(xsign*xi(i)/=ysign*yi(i))exit
          end do
          if(xsign*xi(i)==ysign*yi(i))then
             if(size(yi)>size(xi))then
                allocate(ires(ylow:i))
                ires=yi(ylow:i) ! sign is correct already
                return
             else if(size(xi)>size(yi))then
                allocate(ires(xlow:i))
                ires=xi(xlow:i) ! sign is correct already
                return
             else
                allocate(ires(1))
                ires=0
                return
             end if
          end if
          if(xsign*xi(i)>ysign*yi(i))then
             allocate(biggest(xlow:xup),smallest(ylow:yup))
             biggest=xsign*xi
             smallest=-ysign*yi
             bsign=xsign
          else
             allocate(biggest(ylow:yup),smallest(xlow:xup))
             biggest=ysign*yi
             smallest=-xsign*xi
             bsign=ysign
          end if
       end if
    end if

    call APadd_i(biggest,smallest,ires)
    
    ! overall sign? Sign of 'biggest'.
    ires=bsign*ires
    deallocate(biggest, smallest)
    nullify(biggest, smallest)

  end subroutine APSum_i

  subroutine APMult_i(xin,yin,ires)
    implicit none
    integer(dl), pointer :: xin(:),yin(:)
    integer(dl), pointer :: ires(:)
    integer(dl), pointer :: xi(:),yi(:)
    integer(dl), pointer :: alld(:,:,:), sortd(:,:), mires(:)
    integer :: nx, ny, blow, bup, blow1, bup1, blow2, bup2, finup
    integer :: xlist(size(xin)), i, j, k, m, n
    integer :: ysign, xsign
    
    if(associated(ires))then
       stop 'APMult_i must be called with an empty / unallocated / unassociated pointer for the result.'
    end if

    blow1=lbound(xin,dim=1)
    bup1=ubound(xin,dim=1)
    blow2=lbound(yin,dim=1)
    bup2=ubound(yin,dim=1)
    blow=blow1+blow2-2
    bup=max(bup1,bup2)

    allocate(xi(blow1:bup1),yi(blow2:bup2))
    xsign=1
    ysign=1
    if(any(xin<0))xsign=-1
    if(any(yin<0))ysign=-1
    xi(blow1:bup1)=xsign*xin(blow1:bup1)
    yi(blow2:bup2)=ysign*yin(blow2:bup2)

    xlist=(/ (i,i=blow1,bup1) /)

    finup=bup1+bup2+2
    ! write(*,*)blow1,blow2,blow
    ! write(*,*)bup1,bup2,bup, finup
    nx=size(xi)
    ny=size(yi)
!write(*,*)ny,size(yin) ,lbound(xi,dim=1),ubound(xi,dim=1) ,lbound(yi,dim=1),ubound(yi,dim=1)
    allocate(alld(nx,ny,2))
    alld=0

    do i = 1, ny
       alld(1:nx,i,1)=xi(blow1:bup1)*yi(blow2+i-1)
       alld(1:nx,i,2)=xlist+blow2+i-1
       !write(*,'(16I6)')alld(1:nx,i,2)
    end do

    allocate(sortd(blow:finup,2*(nx+ny)))
    sortd=0

    do i = 1, nx
       do j = 1, ny
          k=(i-1)*nx+j
          m=alld(i,j,2)
          do n = 1, nx+ny-1
             if(sortd(m,n)==0)exit
          end do
          if(sortd(m,n)/=0)stop 'sortd too small in APMult'
          sortd(m,n)=alld(i,j,1)
          if(sortd(m,n)>9)then
             do k = 1, size(sortd(1,:))!nx+ny-1
                if(sortd(m+1,k)==0)exit
             end do
             if(sortd(m+1,k)/=0)stop 'sortd too small in APMult'
             sortd(m+1,k)=int(real(sortd(m,n))/10.)
             sortd(m,n)=sortd(m,n)-10*sortd(m+1,k)
          end if
          
       end do
    end do

    allocate(mires(blow:finup))
    mires=0
    ! Here it goes: sum it all.
    ! Hoping there is still enough space in sortd for new digits
    do i = blow, finup
       do j = 1, nx+ny
          if(sortd(i,j)/=0)then
             mires(i)=mires(i)+sortd(i,j)
             ! intermediate overflowing:
             if(mires(i)>89)then
                if(i+1>finup)stop 'estimate of order of magnitude wrong in APMult'
                mires(i) = mires(i)-90
                do k = 1, size(sortd(1,:))!nx+ny-1
                   if(sortd(i+1,k)==0)exit
                end do
                if(sortd(i+1,k)/=0)then
                   write(*,*)i,j,k
                   write(*,'(16I4)')sortd(i+1,:)
                   stop 'sortd too small in APMult'
                end if
                sortd(i+1,k)=9
             end if
          else
             exit
          end if
       end do
       ! final overflowing:
       if(mires(i)>9)then
          if(i+1>finup)stop 'estimate of order of magnitude wrong in APMult'
          do k = 1, size(sortd(1,:))!nx+ny-1
             if(sortd(i+1,k)==0)exit
          end do
          if(sortd(i+1,k)/=0)then
             write(*,*)i,k, sortd(i+1,k)
             write(*,*)blow1,bup1
             write(*,'(64I3)')xi
             write(*,*)blow2,bup2
             write(*,'(64I3)')yi
             do k = blow, finup
                write(*,'(64I3)')sortd(k,:)
             end do
             stop 'sortd too small in APMult'
          end if
          sortd(i+1,k)=int(real(mires(i))/10.)
          mires(i)=mires(i)-10*sortd(i+1,k)
       end if
       
    end do
  
    ! determine size:
    do i = finup, blow, -1
       if(mires(i)/=0)exit
    end do
    finup=i ! final final upper limit
    do i = blow, finup
       if(mires(i)/=0)exit
    end do
    blow=i

    n=finup-blow+1

    ! size = n digits
    ! How many doubles?
    allocate(ires(blow:finup))
    ires=xsign*ysign*mires(blow:finup)

    deallocate(sortd, mires, alld,xi,yi)
    nullify(sortd, mires, alld,xi,yi)
  end subroutine APMult_i

  subroutine APCp_i(a1,a2)
    implicit none
    integer(dl), pointer :: a1(:), a2(:)
    integer :: low, up

    low=lbound(a1,dim=1)
    up=ubound(a1,dim=1)
    
    allocate(a2(low:up))
    a2=a1

  end subroutine APCp_i

  subroutine APPow_i(xi,n,ires)
    implicit none
    integer(dl), pointer :: xi(:)
    integer(dl), intent(in) :: n
    integer(dl), pointer :: ires(:)
    integer(dl), pointer :: xi2(:), xin(:), xim(:), xi3(:)
    integer :: i
    
    nullify(xi2, xin, xim, xi3)

    if(associated(ires))then
       stop 'APPow_i must be called with an empty / unallocated / unassociated pointer for the result.'
    end if


    ! What is the most clever way? I don't know...

    call APMult_i(xi,xi,xi2)

    if(n==2)then
!       ires => xi2
       call APcp_i(xi2,ires)
       deallocate(xi2)
       return
    end if

    
    if(n/=4)then
       call APMult_i(xi,xi2,xi3)       
    end if
    if(n==3)then
!       ires=>xi3
       call APcp_i(xi3,ires)
       deallocate(xi2, xi3)
       return
    end if

    if(n==4)then
       call APMult_i(xi2,xi2,ires)
       deallocate(xi2)
       return
    end if
    
    if(n==5)then
       call APMult_i(xi2,xi3,ires)
       deallocate(xi2,xi3)
       return
    end if
    
    if(n==6)then
       call APMult_i(xi3,xi3,ires)
       deallocate(xi2,xi3)
       return
    end if
   
    call APMult_i(xi2,xi2,xim)
    call APMult_i(xi3,xim,xin)
    
    if(n==7)then
       ires=>xin
       deallocate(xi2,xi3,xim)
       return
    end if

    deallocate(xim,xi2,xi3)

    do i = 8,n
       call APMult_i(xin,xi,xim)
       deallocate(xin)
       allocate(xin(lbound(xim,dim=1):ubound(xim,dim=1)))
       xin=xim
       deallocate(xim)
    end do

!    ires=>xin
    call APcp_i(xin,ires)
    deallocate(xin)

  end subroutine APPow_i

  subroutine TestAP()
    implicit none
    real(dl), pointer ::  x(:), y(:)
    integer(dl), pointer :: xi(:), yi(:), mi(:), ni(:), yi4(:),yi5(:),res(:),res2(:)
    integer(dl) :: i
    nullify(xi,yi,mi,ni,res,res2,yi4,yi5)

    allocate(x(2),y(2))
    x(:)=0.1234567890123456_dl
    y(:)=0.1234567890123456_dl

    call APDoubleToInteger(x(1),xi)
    call APAdd_i(xi,xi,yi)
    write(*,*)trim(APIntegerToString(yi))
    write(*,*)2.d0*x(1)
    
    yi=-yi
    do i = 2, 12
       call APPow_i(yi,i,mi)
       write(*,'(A)')trim(APIntegerToString(mi))
       deallocate(mi)
       nullify(mi)
    end do

    call APPow_i(yi,4_dl,yi4)
    call APPow_i(yi,5_dl,yi5)

    call APSum_i(yi4,yi5,res)
    write(*,'(A)')trim(APIntegerToString(res))
    yi5=-yi5
    yi4=-yi4
    call APSum_i(yi4,yi5,res2)
    write(*,'(A)')trim(APIntegerToString(res2))
    deallocate(res2)
    nullify(res2)
    yi5=yi5
    yi4=-yi4
    call APSum_i(yi4,yi5,res2)
    write(*,'(A)')trim(APIntegerToString(res2))
    deallocate(res2)
    nullify(res2)
    yi5=-yi5
    yi4=-yi4
    call APSum_i(yi4,yi5,res2)
    write(*,'(A)')trim(APIntegerToString(res2))
    deallocate(res2)
    nullify(res2)
stop

  end subroutine TestAP


end module lltb_apprec
