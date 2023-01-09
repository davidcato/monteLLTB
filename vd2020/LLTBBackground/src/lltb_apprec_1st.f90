! WV's owm arbitrary precision multiplication and addition. Double precisions only.
module lltb_apprec
use lltb_precision

interface APIntDigits
   module procedure APIntDigits_a, APIntDigits_i
end interface

contains

  subroutine APCheckPSize1D(pp,n)
    implicit none
    real(dl), pointer, intent(inout) :: pp(:)
    integer, intent(in) :: n

    if(associated(pp))then
       if(size(pp)/= n)then
          deallocate(pp)
          allocate(pp(n))
       end if
    else
       allocate(pp(n))
    end if
  end subroutine APCheckPSize1D

  subroutine APIntegerToDouble(iarr,res)
    implicit none
    integer(dl), intent(in) :: iarr(2)
    real(dl), intent(out) :: res
    character(len=24) :: xc, xc2, xc3

    write(xc,*)iarr(1)
    xc=trim(adjustl(xc))
    xc=xc(1:1)//'.'//xc(2:)
    
    write(xc2,*)iarr(2)

    xc3=trim(xc)//'E'//trim(adjustl(xc2))
    read(xc3,*)res
write(*,*)iarr
write(*,*)res
  end subroutine APIntegerToDouble

  subroutine APSplit_i(x,x_parts_i)
    implicit none
    real(dl) :: x(:)
    integer(dl) :: x_parts_i(:,:,:)
    character(len=23) :: xc, xc1, xc2
    integer :: i, n, xc1offset, xc2offset, j

    n=size(x)

    if(size(x)/=size(x_parts_i(:,1,1)))stop 'arrays of different size in APSplit'
    if(size(x_parts_i(1,:,1))/=2)stop 'array x_parts of wrong size in APSplit'
    if(size(x_parts_i(1,1,:))/=2)stop 'array x_parts of wrong size in APSplit'

    ! each entry in x has 16 digits
    do i = 1, n
       write(xc,'(E23.16E3)')x(i)
       ! first two characters: "0."
       ! last five: "E+000" or any number..
!       xc1='0.'//xc(3:10)//xc(19:23)
!       xc2='0.00000000'//xc(11:18)//xc(19:23)
       xc1=xc(3:10)
       xc2=xc(11:18)
       read(xc1,*)x_parts_i(i,1,1)
       xc1offset=0
       j=1
       do while (xc1(j:j)=="0")
          xc1offset=xc1offset+1
          j=j+1
          if(j==9)exit
          if(xc1(j:j)/="0")exit         
       end do
       read(xc2,*)x_parts_i(i,2,1)
       xc2offset=0
       j=1
       do while (xc2(j:j)=="0")
          xc2offset=xc2offset+1
          j=j+1
          if(j==9)exit
          if(xc1(j:j)/="0")exit
       end do
       read(xc(20:23),*)x_parts_i(i,1,2)
       x_parts_i(i,1,2)=x_parts_i(i,1,2)-1-xc1offset
       x_parts_i(i,2,2)=x_parts_i(i,1,2)-8-xc2offset
    end do

  end subroutine APSplit_i

  subroutine APSplit_d(x,x_parts)
    implicit none
    real(dl) :: x(:)
    real(dl) :: x_parts(:,:)
    character(len=23) :: xc, xc1, xc2
    integer :: i, n

    n=size(x)

    if(size(x)/=size(x_parts(:,1)))stop 'arrays of different size in APSplit'
    if(size(x_parts(1,:))/=2)stop 'array x_parts of wrong size in APSplit'

    ! each entry in x has 16 digits
    do i = 1, n
       write(xc,'(E23.16E3)')x(i)
       ! first two characters: "0."
       ! last five: "E+000" or any number..
       xc1='0.'//xc(3:10)//xc(19:23)
       xc2='0.00000000'//xc(11:18)//xc(19:23)
       read(xc1,*)x_parts(i,1)
       read(xc2,*)x_parts(i,2)
    end do

  end subroutine APSplit_d
  
  function APIntDigits_a(iarr)
    implicit none
    integer(dl) :: iarr(:)
    integer :: APIntDigits_a(size(iarr))
    integer :: i

    do i = 1, size(iarr)
       APIntDigits_a(i) = APIntDigits_i(iarr(i))
    end do

  end function APIntDigits_a

  function APIntDigits_i(iin)
    implicit none
    integer(dl) :: iin
    integer :: APIntDigits_i
    character(len=64) :: xc
    
    write(xc,*) iin
    
    APIntDigits_i = len(trim(adjustl(xc)))

  end function APIntDigits_i

  ! multiply two numbers, x and y, each can be arbitrary array of doubles.
  ! result is new array res, arbitrary size.
  subroutine APMult(x,y,res)
    implicit none
    real(dl), intent(in) :: x(:),y(:)
    real(dl), pointer, intent(inout) :: res(:)
    integer :: n, nx, ny, i, j, k, m, p, q
    integer(dl), pointer :: x_parts_i(:,:,:), y_parts_i(:,:,:),res_parts_i(:,:,:)
    real(dl), pointer :: x_parts(:,:), y_parts(:,:), res_parts(:,:)
    real(dl), pointer :: res_parts2(:,:), resp2(:)
    
    nx=size(x)
    ny=size(y)
    n=nx+ny
    
!    call CheckPSize1D(res,n)
    allocate(x_parts(nx,2),y_parts(ny,2),res_parts(nx*ny,4))
    allocate(res_parts2(nx*ny,4),resp2(nx*ny*4))
    allocate(x_parts_i(nx,2,2),y_parts_i(ny,2,2),res_parts_i(nx*ny,4,2))
    
    call APSplit_I(x,x_parts_i)
    call APSplit_i(y,y_parts_i)
    call APSplit_d(x,x_parts)
    call APSplit_d(y,y_parts)
    
    ! i: x
    ! j: y
    do i = 1, nx
       do j = 1, ny 
          
          m=(i-1)*ny+j

          res_parts(m,1:2)=x_parts(i,:)*y_parts(j,:)
          res_parts(m,3:4)=x_parts(i,1:2)*y_parts(j,2:1:-1)
          !write(*,*)x_parts(i,1),y_parts(j,1),x_parts(i,1)*y_parts(j,1)
          ! multiply integers:
          res_parts_i(m,1:2,1)=x_parts_i(i,:,1)*y_parts_i(j,:,1)
          res_parts_i(m,3:4,1)=x_parts_i(i,1:2,1)*y_parts_i(j,2:1:-1,1)
          ! add exponents:
          res_parts_i(m,1:2,2)=x_parts_i(i,:,2)+y_parts_i(j,:,2)
          res_parts_i(m,3:4,2)=x_parts_i(i,1:2,2)+y_parts_i(j,2:1:-1,2)

          !write(*,*)x_parts_i(i,1,1),y_parts_i(j,1,1)
          !write(*,*)res_parts_i((i-1)*ny+j,1,:)
          do k = 1, 4
             if(k<=2)then
                p=k
                q=k
             else
                p=-k+5
                q=k-2
             end if

             if(APIntDigits(res_parts_i(m,k,1))>APIntDigits(x_parts_i(i,q,1))+APIntDigits(y_parts_i(j,p,1))-1)then
                !write(*,*)APIntDigits(res_parts_i(m,k,1)),APIntDigits(x_parts_i(i,q,1))+APIntDigits(y_parts_i(j,p,1))-1
                !write(*,*)res_parts_i(m,k,1),x_parts_i(i,q,1),y_parts_i(j,p,1)
                res_parts_i(m,k,2)=res_parts_i(m,k,2)+1
             end if
             call APIntegerToDouble(res_parts_i(m,k,:),res_parts2(m,k))
          end do
       end do
    end do

    ! now total result = sum(res_parts).
    ! each entry in res_parts has 16 significant digits.
    do i = 1, 4
       resp2((i-1)*nx*ny+1:nx*ny*i) = res_parts2(:,i)
    end do

write(*,*)res_parts
write(*,*)res_parts2
write(*,*)resp2

    call APSum(resp2,res)


    deallocate(x_parts,y_parts,res_parts)
    deallocate(x_parts_i,y_parts_i,res_parts_i, res_parts2,resp2)
  end subroutine APMult


  ! sum over arbitrary array of doubles, result is arbitrary array res.
  subroutine APSum(array, res)
    implicit none
    real(dl), intent(in) :: array(:)
    real(dl), pointer, intent(inout) :: res(:)
    integer, pointer :: alld(:,:), sortd(:,:), ires(:)
    integer :: i, j, n, k, offset, blow,bup, finup
    character(len=23) :: xc, xc1, xc2

    n=size(array)*16
    allocate(alld(n,2))
    ! Make list of single digits with their exponent:
    ! alld: list of (digit, decimal position)
    do i = 1, size(array)
       write(xc,'(E23.16E3)')array(i)
       ! first two characters: "0."
       ! last five: "E+000" or any number..

       read(xc(20:23),*)offset

       do j = 1, 16
          k = (i-1)*16+j
          read( xc(2+j:2+j) , *) alld(k,1)
          alld(k,2) = offset - j 
          ! since first digit is x in 0.x
       end do
    end do
    
    blow=minval(alld(:,2))
    bup=maxval(alld(:,2))

    ! estimate rough upper limit of result:
    finup = nint(log10(sum(array)))+1

    allocate(ires(blow:finup))
    ires=0

    allocate(sortd(blow:finup,n))

    sortd=0
    
    ! Make huge array with exponents as indices, collecting digits there.
    do i = 1, n
       j = alld(i,2)
       do k = 1, n
          if(sortd(j,k)==0)exit
       end do
       sortd(j,k)=alld(i,1)
    end do

    
    ! Here it goes: sum it all.
    ! Hoping there is still enough space in sortd for new digits
    do i = blow, finup
       do j = 1, n
          if(sortd(i,j)/=0)then
             ires(i)=ires(i)+sortd(i,j)
             ! intermediate overflowing:
             if(ires(i)>89)then
write(*,*)'taking off:',ires(i)
                if(i+1>finup)stop 'estimate of order of magnitude wrong in APSum'
                ires(i) = ires(i)-90
                do k = 1, n
                   if(sortd(i+1,k)==0)exit
                end do
                if(sortd(i+1,k)/=0)stop 'sortd too small in APSum'
                sortd(i+1,k)=9
write(*,*)'newval:',ires(i)
             end if
          else
             exit
          end if
       end do
       ! final overflowing:
       if(ires(i)>9)then
          if(i+1>finup)stop 'estimate of order of magnitude wrong in APSum'
          do k = 1, n
             if(sortd(i+1,k)==0)exit
          end do
          if(sortd(i+1,k)/=0)stop 'sortd too small in APSum'
          sortd(i+1,k)=int(real(ires(i))/10.)
          ires(i)=ires(i)-10*sortd(i+1,k)
       end if
       
    end do

    ! now ires is a list of digits with there exponents, 
    ! of the final result.

    ! determine size:
    do i = finup, blow, -1
       if(ires(i)/=0)exit
    end do
    n=i-blow+1
    finup=i ! final final upper limit
    ! size = n digits
    ! How many doubles?

    k=int(real(n)/16.)+1
    if(n==(k-1)*16)k=k-1 ! in case n=int*16

    call APCheckPSize1D(res,k)
    res=0.d0

    ! And now loop over ires...
    do i = 1, k
       offset=(i-1)*16
       write(xc,'(I1,".",15I1)')ires(finup-offset:max(finup-offset-15,blow):-1)
       write(xc1,*)finup-offset
       xc=trim(adjustl(xc))//'E'//trim(adjustl(xc1))
       read(xc,*)res(i)
write(*,*)'xc:',xc
    end do

    deallocate(alld)

  end subroutine APSum


  subroutine TestAP()
    implicit none
    real(dl), pointer :: res(:), x(:), y(:), res2(:)


    allocate(x(2),y(2))
    x(:)=0.1234567890123456_dl
    y(:)=0.1234567890123456_dl

    call APmult(x,y,res)
    call APmult(res,res,res2)
    
    write(*,*)res
    write(*,*)sum(x)*sum(y)

    write(*,*)res2
    write(*,*)(sum(x)*sum(y))**2

stop

  end subroutine TestAP

end module lltb_apprec
