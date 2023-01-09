  ! WV
  ! asplint (array splint): splint as in numerical recipes
  ! over time-parameter x, but array yy 
  ! (e.g. several functions that depend on same parameter x.)
  subroutine asplint(x,y,xx,yy,ddy)
    implicit none
    integer :: klo,khi,k,n
    real(dl) :: x(:),y(:,:),ddy(:,:)
    real(dl) :: xx,yy(:)
    real(dl) :: a,b,h
    integer :: counter
    logical :: nosuccess

    n=size(x)

    counter=0

    klo=1
    khi=n


    ! Check if this loop gets stalled:
    !    if (khi==klo) stop 'STOP - gT_splint called with empty array?'
    ! nosucces: not found lowest possible answer
    nosuccess=.true.
    do while (nosuccess)
       do while (khi-klo.gt.1)
          k=int((khi+klo)/2)

          if(x(k).gt.xx) then
             khi=k
          else
             klo=k
          endif

          counter=counter+1
          if(counter>2*n .or. counter.gt.MaxSteps) then ! stop 'STOP - gT_splint cannot find right interval.'
             ! This may perhaps happen if x(k) is not real.
             ! In principal it is impossible that this happens, since 'else'
             ! always applies when x is not real.
             write(*,*) 'gT_splint cannot find right interval.'
             !error=.true.
             !return
             !          stop
             VoidError=.true.
             return
          end if
       end do
       ! Next lines: check is the found k corresponds to the highest option.
       nosuccess=.false.
       if(any(x(khi:size(x)).le.xx))then
          if (size(x)-khi .gt. 0) then
             klo=smallest(x(khi:size(x)))+klo
             khi=size(x)
             nosuccess=.true.
             counter =0
          else
             khi=size(x)
             klo=khi-1
             nosuccess=.false.
          end if
       else
          nosuccess=.false.
       end if

    end do


    h=x(khi)-x(klo)
    if (h.eq.0.) then
       write(*,*)'Problem: h=0 in splint of rr(z)'
!       stop
       VoidError = .true.
       return
    end if
    a=(x(khi)-xx)/h
    b=(xx-x(klo))/h
    yy=a*y(:,klo)+b*y(:,khi)+((a**3-a)*ddy(:,klo)+(b**3-b)*ddy(:,khi))*(h**2)/6.d0
    
    if(any(abs(y(:,khi)-y(:,klo))>1.d-30))then
       do k = 1, size(yy)
!   write(*,'(4ES25.16,I4)')yy(k),y(k,khi),y(k,klo),(yy(k)-y(k,klo))-(y(k,khi)-y(k,klo)),k
          if((yy(k)>y(k,khi)) .or. (yy(k)<y(k,klo)))then
!!$             write(*,*)'ddy broken in asplint: result out of bounds.'
!!$             !          write(*,*)'ddy1:',ddy(1,:)
!!$             !          write(*,*)'ddy2:',ddy(2,:)
!!$             write(*,'(4ES25.16,I4)')yy(k),y(k,khi),y(k,klo),(yy(k)-y(k,klo))-(y(k,khi)-y(k,klo)),k
!!$             write(*,'(20ES25.16,I4)')a,b,x(khi),xx,x(klo)
             yy(k)=a*y(k,klo)+b*y(k,khi)
!!$             write(*,'(4ES25.16,I4)')yy(k),y(k,khi),y(k,klo),(yy(k)-y(k,klo))-(y(k,khi)-y(k,klo)),k
!!$             do n = 2, size(x)
!!$                write(*,'(20ES25.16)')x(n), x(n)-x(n-1)
!!$             end do
             !          stop
             !write(*,*)'returning very rough estimate. This should not happen in your good fit models.'
          end if
       end do
    end if
    if(VoidFEedback>0)then
       if(any(isnan(yy)))then
          call VoidErrorReport()
          write(*,*)'isnan(yy) in asplint. Not always a problem.'
       end if
    end if
  end subroutine asplint
!!$




  ! WV
  ! We only have to do this once for each model.
  ! Reports all indices k for which x(k) is a local
  ! minimum (kminima(nminima))
  ! or  maximum (kmaxima(nmaxima))
  subroutine GetLocalMinMax(x,kmaxima,kminima,nmaxima,nminima)
    implicit none
    real(dl) :: x(:)
    integer :: kmaxima(:),kminima(:),nmaxima,nminima,k, n, i
    logical :: plateau

    n=size(x)
    if(size(kmaxima).ne.n   .or.   size(kminima).ne. n) stop 'Wrong sizes of kmaxima / kminima in GetLocalMinMax'

    if(n==1)then
       kmaxima = 1
       kminima = 1
       nmaxima = 1
       nminima = 0
       write(*,*)'GetLocalMinMax called with array of size 1.'
       return
    end if

    ! First find all local maxima and minima
    nmaxima=0
    nminima=0
    
    if( x(1).lt.x(2) )then
       nminima = nminima+1
       kminima(nminima) = 1
    else if( x(1).gt.x(2) )then
       nmaxima = nmaxima+1
       kmaxima(nmaxima) = 1
    end if

    do k=2,n-1
       if( (x(k-1).lt.x(k)) .and. (x(k+1).lt.x(k)) )then ! local maximum
          nmaxima = nmaxima+1
          kmaxima(nmaxima) = k
       else if( (x(k-1).gt.x(k)) .and. (x(k+1).gt.x(k)) )then ! local minimum
          nminima = nminima+1
          kminima(nminima) = k
       end if
    end do

    plateau = .true.
    i = n
    do while (plateau)
       
       if(i.le.0)then
          VoidError = .true.
          return
       end if

       if(x(i) .eq. x(i-1))then
          plateau = .true.
          i = i - 1
          cycle
       else if( x(i).lt.x(i-1) )then
          nminima = nminima+1
          kminima(nminima) = n
       else if( x(i).gt.x(i-1) )then
          nmaxima = nmaxima+1
          kmaxima(nmaxima) = n
       end if
       plateau = .false.
       exit
    end do
    
  end subroutine GetLocalMinMax
  

  ! WV
  ! masplint: multivalued asplint. 
  ! returns all possible outcomes if there are more than one.
  subroutine masplint(x,y,xx,yy,ddy,inkmaxima,inkminima)
    implicit none
    !    integer, save :: klom,khim
    real(dl) :: x(:),y(:,:),ddy(:,:)
    integer :: klo,khi,k,n, nmaxima, nminima, ndom
    integer, allocatable :: kmaxima(:), kminima(:), kmaxima2(:), kminima2(:), domains(:,:)
    integer, pointer, optional :: inkmaxima(:), inkminima(:)
    real(dl) :: xx
    real(dl), pointer :: yy(:,:)
    real(dl), allocatable :: tempyy(:,:)
    real(dl) :: a,b,h
    integer :: counter, thisk, domcount, resultcount, nvar
    logical ::  up, redomaxmin, dothisloop

    stop 'do not use masplint'

    n=size(x)
    nvar = size(y(:,1))

    redomaxmin = .false.
    dothisloop = .true.


    counter=0

       
    if(present(inkmaxima).and.present(inkminima) .and. .not. (redomaxmin))then
       allocate(kmaxima(size(inkmaxima)))
       allocate(kminima(size(inkminima)))
       kmaxima=inkmaxima
       kminima=inkminima
       nmaxima = size(inkmaxima)
       nminima = size(inkminima)
    else
       allocate(kmaxima2(n),kminima2(n))
       call GetLocalMinMax(x,kmaxima2,kminima2,nmaxima,nminima)
       allocate(kmaxima(nmaxima),kminima(nminima))
       kmaxima=kmaxima2
       kminima=kminima2
       deallocate(kmaxima2,kminima2)
    end if
    
       ! Now we have kmin1 < kmax1 < kmin2 < kmax2 etc. etc.
       ! call that     -- dom1 -- dom2 -- dom3   etc. etc.
       ! number of domains = 2 nmaxima - 1 if highest kmax == n  && lowest kmin = 1
       ! number of domains = 2 nmaxima  if highest kmin == n  && lowest kmin = 1
       ! number of domains = 2 nmaxima - 2 if highest kmax == n  && lowest kmax = 1
       ! number of domains = 2 nmaxima - 1 if highest kmin == n  && lowest kmax = 1
       
!!$       if(kmaxima(nmaxima)==n .and. kminima(1)==1)then
!!$          ndom = 2 * nmaxima - 1
!!$          dothisloop = .false.
!!$          exit
!!$       else if(kminima(nminima)==n .and. kminima(1)==1)then
!!$          ndom = 2 * nmaxima
!!$          dothisloop = .false.
!!$          exit
!!$       else if(kmaxima(nmaxima)==n .and. kmaxima(1)==1)then
!!$          ndom = 2 * nmaxima - 2
!!$          dothisloop = .false.
!!$          exit
!!$       else if(kminima(nminima)==n .and. kmaxima(1)==1)then
!!$          ndom = 2 * nmaxima - 1
!!$          dothisloop = .false.
!!$          exit
!!$       else
!!$          write(*,*)'Bug in masplint occurred:'
!!$          write(*,*)kmaxima,'.'
!!$          write(*,*)kminima,'.'
!!$          write(*,*)n,'.'
!!$          if(redomaxmin)stop 'Bug in masplint.'
!!$          redomaxmin = .true.
!!$          deallocate(kmaxima,kminima,stat=mystat)
!!$       end if
!!$    end do

    if(max(maxval(kminima),maxval(kmaxima)).ne.n)then
       write(*,*)'Bug in masplint.'
       write(*,*)kmaxima,'.'
       write(*,*)kminima,'.'
       write(*,*)n,'.'
       write(*,*)x,'.'
    end if

    ndom = nminima + nmaxima - 1

    ! Collect all domain border in domains(:,2)
    allocate(domains(ndom,2))
    if(kminima(1)==1)then
       do k=1,ndom
          if(modulo(k,2).ne.0)then ! odd k: ascending
             thisk= k/2 +1
             domains(k,1)=kminima(thisk)
             domains(k,2)=kmaxima(thisk)
          else ! even k: descending
             thisk= k/2 
             domains(k,1)=kmaxima(thisk)
             thisk= k/2 +1
             domains(k,2)=kminima(thisk)
          end if
       end do
    else
       do k=1,ndom
          if(modulo(k,2).ne.0)then ! odd k: descending
             thisk= k/2 +1
             domains(k,1)=kmaxima(thisk)
             domains(k,2)=kminima(thisk)
          else ! even k: ascending
             thisk= k/2 
             domains(k,1)=kminima(thisk)
             thisk= k/2 +1
             domains(k,2)=kmaxima(thisk)
          end if
       end do
    end if


    ! Loop over domains
    allocate(tempyy(1:nvar,1:size(x)))
    counter = 0
    resultcount = 0
    do domcount = 1, ndom
       if( (xx.lt.min(x(domains(domcount,1)),x(domains(domcount,2)))) .or. (xx.gt.max(x(domains(domcount,1)),x(domains(domcount,2)))) ) cycle
       ! Continue only if x is actually in domain
       klo=domains(domcount,1)
       khi=domains(domcount,2)
       
       if(x(klo).lt.x(khi))then
          up = .true.
       else
          up = .false.
       end if
       
       ! Check if this loop gets stalled:
       !    if (khi==klo) stop 'STOP - gT_splint called with empty array?'
       ! nosucces: not found lowest possible answer
       do while (abs(khi-klo).gt.1)
          k=int((khi+klo)/2)
          
          if(x(k).gt.xx) then
             if(up)then
                khi=k
             else
                klo=k
             end if
          else
             if(up)then
                klo=k
             else
                khi=k
             end if
          endif
          
          counter=counter+1
          if(counter>2*n .or. counter.gt.MaxSteps) then ! stop 'STOP - gT_splint cannot find right interval.'
             ! This may perhaps happen if x(k) is not real.
             ! In principal it is impossible that this happens, since 'else'
             ! always applies when x is not real.
             write(*,*) 'gT_splint cannot find right interval.'
             !error=.true.
             !return
             !          stop
             VoidError=.true.
             !return
             exit
          end if
       end do



       h=x(khi)-x(klo)
       if (h.eq.0.) then
          write(*,*)'Problem: h=0 in splint of rr(z)'
          VoidError = .true.
          !return
          exit
       end if
       a=(x(khi)-xx)/h
       b=(xx-x(klo))/h
       
       resultcount = resultcount + 1
       tempyy(:,resultcount)=a*y(:,klo)+b*y(:,khi)+((a**3-a)*ddy(:,klo)+(b**3-b)*ddy(:,khi))*(h**2)/6.d0

       if(any((tempyy(:,resultcount)-y(:,klo))/(y(:,khi)-y(:,klo))>1._dl))then
          write(*,*)'ddy broken in masplint: result out of bounds.'
!          write(*,*)'ddy1:',ddy(1,:)
!          write(*,*)'ddy2:',ddy(2,:)
          tempyy(:,resultcount)=a*y(:,klo)+b*y(:,khi)
!write(*,*)tempyy(:,resultcount),y(:,khi),y(:,klo)
!          stop
          write(*,*)'returning very rough estimate. This should not happen in your good fit models.'
       end if
    end do

!    deallocate(yy,stat=mystat)
    if(.not.voiderror .and. resultcount>0)then
       deallocate(yy)
       allocate(yy(1:nvar,1:resultcount))
       yy(1:nvar,1:resultcount)=tempyy(1:nvar,1:resultcount)
    end if

    deallocate(kminima,kmaxima,domains,tempyy)

!!$    ! testing it:
!!$    if(resultcount.gt.1)then
!!$       write(*,*)'Found more than one solution:'
!!$       do k = 1,resultcount
!!$          write(*,*)'--------------------------------------------------------------'
!!$          write(*,'(A,2I4)')'klo,khi:',alllos(k),allhis(k)
!!$          write(*,'(A,3E14.5)')' z:',x(alllos(k)),xx,x(allhis(k))
!!$          write(*,'(A,3E14.5)')' y(1):',y(1,alllos(k)),yy(1,k),y(1,allhis(k))
!!$          write(*,'(A,3E14.5)')' y(2):',y(2,alllos(k)),yy(2,k),y(2,allhis(k))
!!$          write(*,*)'--------------------------------------------------------------'
!!$       end do
!!$       stop
!!$    end if
    
    if(resultcount.eq.0)then !stop 'masplint cannot extrapolate.'
!       deallocate(yy, stat=mystat)
       deallocate(yy)
       allocate(yy(1:nvar,1))
       
       call asplint(x,y,xx,yy(1:nvar,1),ddy)
    end if
    
  end subroutine masplint
!!$

  function largest(ar)
    implicit none
    integer :: largest
    integer :: i
    real(dl) :: ar(:), this
    
    this=ar(1)
    largest = size(ar)
    if(size(ar).eq.1)return

    do i=1,size(ar)
       if(ar(i).gt.this)then
          largest=i
          this=ar(i)
       end if
    end do
    
  end function largest

  function smallest(ar)
    implicit none
    integer :: smallest
    integer :: i
    real(dl) :: ar(:), this
    
    this=ar(1)
    smallest = 1
    
    if(size(ar).eq.1)return

    do i=2,size(ar)
       if(ar(i).lt.this)then
          smallest=i
          this=ar(i)
       end if
    end do
    
  end function smallest

  ! Spline from Numerical Recipes
  ! input: x(n), y(n), a
  ! output:: ddy(n)
  Subroutine aspline(xin,yin,ddyout)
    implicit none
    integer :: n,m,i!,j
    real(dl), intent(in), target :: xin(:),yin(:,:)
    real(dl), intent(out), target :: ddyout(:,:)
    real(dl), allocatable :: newx(:), newy(:,:), newddy(:,:)
    integer, allocatable :: copymem(:)
    real(dl), pointer :: x(:),y(:,:)
    real(dl), pointer :: ddy(:,:)
    real(dl) :: sig_spline
    real(dl), allocatable :: u_spline(:,:),p_spline(:)
    logical :: didallocation
    integer :: ntot, ccounter

    didallocation = .false.

    n = size(xin)
    m=size(yin(:,1))

    x => xin
    y => yin
    ddy => ddyout

    ntot = n

    
    if(any(x(2:n)-x(1:n-1)==0.d0) .or. any(x(3:n)-x(1:n-2)==0.d0) )then
       didallocation = .true.
       allocate(newx(n),newy(m,n),newddy(m,n),copymem(n))
       ntot = 1
       newx(ntot) = x(ntot)
       newy(:,ntot) = y(:,ntot)
       newddy(:,ntot) = ddy(:,ntot)
       ccounter = 0
       copymem = 0
       do i = 2, n
          if((x(i)-x(i-1)==0.d0).or.(x(i)-x(max(i-2,1))==0.d0))then
             ccounter = ccounter + 1
             copymem(ccounter) = i
             cycle
          else
             ntot = ntot + 1
             newx(ntot) = x(i)
             newy(:,ntot) = y(:,i)
             newddy(:,ntot) = ddy(:,i)
          end if
       end do
       nullify(x,y,ddy)
       allocate(x(ntot), y(m,ntot), ddy(m,ntot))

       x=newx(1:ntot)
       y=newy(:,1:ntot)
       ddy=newddy(:,1:ntot)

       deallocate(newx,newy,newddy)
    end if

    n=ntot


    allocate(u_spline(m,n))
    !    allocate(sig_spline(m))
    allocate(p_spline(m))



    ddy(:,1)=0.d0
    u_spline(:,1)=0.d0
    do i=2,n-1
       sig_spline=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p_spline=sig_spline*ddy(:,i-1)+2.d0
       ddy(:,i)=(sig_spline-1.)/p_spline
       u_spline(:,i)=(6.d0*((y(:,i+1)-y(:,i))&
            /(x(i+1)-x(i))-(y(:,i)-y(:,i-1))&
            /(x(i)-x(i-1)))/(x(i+1)-x(i-1))&
            -sig_spline*u_spline(:,i-1))/p_spline
!       write(*,'(20ES25.16)')u_spline(:,i),y(:,i+1)-y(:,i),x(i+1)-x(i),y(:,i)-y(:,i-1),x(i)-x(i-1),x(i+1)-x(i-1),p_spline,sig_spline
    end do
    ddy(:,n)=0.d0
    do i=n-1,1,-1
       ddy(:,i)=ddy(:,i)*ddy(:,i+1)+u_spline(:,i)
    end do


    deallocate(u_spline)
    deallocate(p_spline)

    if(any(isnan(ddy)))then
       call VoidErrorReport()
       write(*,*)x
       do i = 1, size(y(:,1))
          write(*,*)y(i,:)
       end do
       do i = 2, n
          write(*,*)x(i)-x(i-1), x(min(n,i+1))-x(i-1)
       end do
       do i = 1, size(ddy(:,1))
          write(*,*)ddy(i,:)
       end do
       if(didallocation)write(*,*)copymem(1:ntot)
       call wlltb_error('isnan in aspline')
       stop 'isnan in aspline'
    end if

    
    if(didallocation)then
       n=size(xin)
       ! First copy one to one.
       ddyout(:,1:ntot)=ddy(:,1:ntot)
       ! No more use: deallocate straight away.
       deallocate(x,y,ddy)
       ! Then shift for each duplicate entry
       do i = 1, n
          if(copymem(i)==0)exit
          ddyout(:,copymem(i):ntot+i)=ddyout(:,copymem(i)-1:ntot+i-1)
       end do

    else
       nullify(x,y,ddy)
    end if


  end Subroutine aspline


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



  FUNCTION voidsplint(xa,ya,y2a,x)
    IMPLICIT NONE	
    REAL(dl), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
    REAL(dl), INTENT(IN) :: x
    REAL(dl) :: voidsplint
    INTEGER :: khi,klo,n, i
    REAL(dl) :: a,b,h
!    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    n=size(xa)
!    klo=max(min(locate(xa,x),n-1),1)
    ! WV:
    do i = 1,n
       if(xa(1).lt.xa(2))then
          if(xa(i).ge.x)exit
       else
          if(xa(i).lt.x)exit
       end if
    end do
!    khi=klo+1
    if(i>n)i=n
    khi = i
    klo = i-1
    if(klo.lt.1)then
       klo = 1
       khi = 2
    end if
    ! End WV
    h=xa(khi)-xa(klo)
!    if (h == 0.0)
!    call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    voidsplint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp

  END FUNCTION voidsplint

!----------------------------------------------------!
! Non mathematical utils:
!----------------------------------------------------!

  subroutine Voidsort(x1,x2)
    implicit none
    real(dl) :: x1(:),x2(:)
    integer :: i, n, nswaps, counter
    n = size(x1)
    counter = 0
    nswaps = 1
    do while(nswaps.ne.0)
       nswaps = 0
       do i = 1, n-1
          if(x1(i).gt.x1(i+1))then
             call voidswap(x1(i),x1(i+1))
             call voidswap(x2(i),x2(i+1))
             nswaps = nswaps+1
          end if
       end do
       counter = counter + 1
       if(counter.gt.n**2)then
          voiderror = .true.
          if(voidfeedback>1)write(*,*)'Voidsort was hanging.'
          return
       end if
    end do
  end subroutine Voidsort

  subroutine voidswap(x1,x2)
    implicit none
    real(dl) :: x1,x2,xmem
    xmem = x1
    x1=x2
    x2=xmem
  end subroutine voidswap


   function GetNearestIndex(array2,val2)
     implicit none
     integer :: GetNearestIndex, kmin, kmax, k
     real(dl) :: array2(:),array(size(array2)), val, val2
     kmin = 1
     kmax = size(array2)
     
     array = dabs(array2)
     val = dabs(val2)
     
     do while (kmax - kmin .gt. 1)
        k = int((kmax+kmin)/2)
        if(array(k).gt.val)then
           kmax=k
        else
           kmin=k
        end if
     end do
     
     GetNearestIndex=kmax
     
   end function GetNearestIndex
   
   
   ! Count the number of space-separated items in mystr
   function CountItemsIn(mystr2)
     implicit none
     integer :: CountItemsIn, i, j
     character(len=*), intent(in) :: mystr2
     character(len=len(mystr2)+2) :: mystr
     
     mystr=trim(adjustl(mystr2))//'  '
     
     j=0
     do i=1,len(trim(mystr))+2
        if(mystr(i:i).eq.' '.or.mystr(i:i).eq.CHAR(9))then
           if(i.eq.1)then ! This should never happen
              j=j+1
           else if((i.gt.1).and.(mystr(i-1:i-1).ne.' '))then
              j=j+1
           end if
        end if
     end do
     CountItemsIn = j 
   end function CountItemsIn
   
   subroutine VoidRedMessage(msg)
    implicit none
    character(len=*) :: msg
    
!    Black      0;30       Dark Gray    1;30 
!    Red        0;31       Bold Red     1;31 
!    Green      0;32       Bold Green   1;32 
!    Yellow     0;33       Bold Yellow  1;33 
!    Blue       0;34       Bold Blue    1;34  
!    Purple     0;35       Bold Purple  1;35 
!    Cyan       0;36       Bold Cyan    1;36 
!    Light Gray 0;37       White        1;37
    
    write(0,'(3A)')char(27)//'[01;31m',trim(msg),char(27)//'[01;0m'

    
  end subroutine VoidRedMessage
   
   
   subroutine VoidBigAlert(charstring)
     character(len=*) :: charstring
     character(len=512), allocatable :: copyallalerts(:)
     if(.not.allocated(allalerts))allocate(allalerts(1))
     if(.not.any(allalerts==charstring))then
        allocate(copyallalerts(size(allalerts)))
        copyallalerts = allalerts
        deallocate(allalerts)
        allocate(allalerts(size(copyallalerts)+1))
        allalerts(1:size(allalerts)-1) = copyallalerts
        allalerts(size(allalerts))=charstring
        
        write(*,*)''
        write(*,'(A)')'***************** ---------------------------------------------------'
        write(*,'(A)')'***************** BIG VOID ALERT: '//charstring
        write(*,'(A)')'***************** ---------------------------------------------------'
        write(*,*)''
        
        deallocate(copyallalerts)
     end if
     
     return
     
   end subroutine VoidBigAlert

   
   subroutine VoidAlertAtMostOnceEveryModel(charstring)
     character(len=*), optional :: charstring
     character(len=512), allocatable :: copyallalerts(:)

  if(present(charstring))then
     if(.not.allocated(oncealerts))allocate(oncealerts(1))
     if(.not.any(oncealerts==charstring))then
        allocate(copyallalerts(size(oncealerts)))
        copyallalerts = oncealerts
        deallocate(oncealerts)
        allocate(oncealerts(size(copyallalerts)+1))
        oncealerts(1:size(oncealerts)-1) = copyallalerts
        oncealerts(size(oncealerts))=charstring
        
!        write(*,*)''
        write(*,'(A)')'* ---------------------------------------------------'
        write(*,'(A)')'* '//charstring
        write(*,'(A)')'* ---------------------------------------------------'
!        write(*,*)''
        
        deallocate(copyallalerts)
     end if
  else
    if(allocated(oncealerts))deallocate(oncealerts)
  end if
     return
     
   end subroutine VoidAlertAtMostOnceEveryModel

   subroutine VoidErrorReport(CMB,GetMsg)
     type(FakeCMBParams), optional :: CMB
     type(FakeCMBParams), save :: MyCMB
     character(len=1) :: mynl
     character(len=*), optional :: GetMsg
     character(len=10000) :: msg
     character(len=35) :: thisnum
     mynl = NEW_LINE('A')

     if(present(CMB))then
        MyCMB = CMB
        return
     end if

     msg = '***********************************************'//mynl
     msg = trim(msg)//'------------- VOID ERROR ----------------------'//mynl


     write(thisnum,'(ES31.21)')myCMB%VCP%delta0
     msg = trim(msg)//'CMB%VCP%delta0: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%H0
     msg = trim(msg)//'CMB%H0: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%VCP%zB
     msg = trim(msg)//'CMB%VCP%zB: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%omk
     msg = trim(msg)//'CMB%omk: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%omdmh2
     msg = trim(msg)//'CMB%omdmh2: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%ombh2
     msg = trim(msg)//'CMB%ombh2: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%VCP%alpha
     msg = trim(msg)//'CMB%VCP%alpha: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%VCP%beta
     msg = trim(msg)//'CMB%VCP%beta: '//trim(thisnum)//mynl

     msg = trim(msg)//mynl

     write(thisnum,'(ES31.21)')myCMB%omch2
     msg = trim(msg)//'CMB%omch2: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%omnuh2
     msg = trim(msg)//'CMB%omnuh2: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%zre
     msg = trim(msg)//'CMB%zre: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%nufrac
     msg = trim(msg)//'CMB%nufrac: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%reserved(1)
     msg = trim(msg)//'CMB%reserved(1) (tau): '//trim(thisnum)//mynl

     write(thisnum,'(ES31.21)')myCMB%omv
     msg = trim(msg)//'CMB%omv: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')myCMB%w
     msg = trim(msg)//'CMB%w: '//trim(thisnum)//mynl



     msg = trim(msg)//mynl

     write(thisnum,'(ES31.21)')myCMB%ombh2
     msg = trim(msg)//'param[omegabh2] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl
     write(thisnum,'(ES31.21)')myCMB%omdmh2
     msg = trim(msg)//'param[omegadmh2] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl
     write(thisnum,'(ES31.21)')myCMB%H0
     msg = trim(msg)//'param[H0] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl
     write(thisnum,'(ES31.21)')myCMB%w
     msg = trim(msg)//'param[w] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl
     if(myCMB%omv/=0.e0 .and. myCMB%omk==0.d0)then
        write(thisnum,'(ES31.21)')1.e-16
     else 
        write(thisnum,'(ES31.21)')myCMB%omk
     end if
     msg = trim(msg)//'param[omegak] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl
     write(thisnum,'(ES31.21)')myCMB%VCP%zB
     msg = trim(msg)//'param[zb] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl
     write(thisnum,'(ES31.21)')myCMB%VCP%delta0
     msg = trim(msg)//'param[delta0] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl
     write(thisnum,'(ES31.21)')myCMB%VCP%L2_1
     msg = trim(msg)//'param[L21] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl
     write(thisnum,'(ES31.21)')myCMB%VCP%d0_2
     msg = trim(msg)//'param[spare] = '//trim(thisnum)//' '//trim(thisnum)//' '//trim(thisnum)//' 0 0'//mynl

     msg = trim(msg)//mynl


     write(thisnum,'(ES31.21)')VP%zB
     msg = trim(msg)//'VP%zB: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')VP%L
     msg = trim(msg)//'VP%L: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')VP%delta0
     msg = trim(msg)//'VP%delta0: '//trim(thisnum)//mynl
     write(thisnum,'(ES31.21)')VP%kmax
     msg = trim(msg)//'VP%kmax: '//trim(thisnum)//mynl
     

     msg = TRIM(MSG)//'-----------------------------------------------'//mynl
     msg = TRIM(MSG)//'***********************************************'//mynl

     if(present(GetMsg))then
        GetMsg=''
        GetMsg=msg
     else
!        write(*,'(A)')trim(msg)
        write(0,'(A)')trim(msg)
     end if

   end subroutine VoidErrorReport


   subroutine maspline(x,y,ddy,inkmaxima,inkminima)
     real(dl) :: x(:),y(:,:),ddy(:,:)
     integer :: klo,khi,k,n, nmaxima, nminima, ndom
     integer, allocatable :: kmaxima(:), kminima(:), kmaxima2(:), kminima2(:), domains(:,:)
     integer, pointer, optional :: inkmaxima(:), inkminima(:)
     real(dl) :: xx
!     real(dl), pointer :: yy(:,:)
!     real(dl), allocatable :: tempyy(:,:)
     real(dl) :: a,b,h
     integer :: counter, thisk, domcount, resultcount, nvar, direction
     logical ::  up, redomaxmin, dothisloop
     
     stop 'Do not use maspline.'

     n=size(x)
     nvar = size(y(:,1))
     
     counter=0
     
     
     if(present(inkmaxima).and.present(inkminima))then
        allocate(kmaxima(size(inkmaxima)))
        allocate(kminima(size(inkminima)))
        kmaxima=inkmaxima
        kminima=inkminima
        nmaxima = size(inkmaxima)
        nminima = size(inkminima)
     else
        allocate(kmaxima2(n),kminima2(n))
        call GetLocalMinMax(x,kmaxima2,kminima2,nmaxima,nminima)
        allocate(kmaxima(nmaxima),kminima(nminima))
        kmaxima=kmaxima2
        kminima=kminima2
        deallocate(kmaxima2,kminima2)
     end if
     
     ndom = nminima + nmaxima - 1

     ! Collect all domain border in domains(:,2)
     allocate(domains(ndom,2))
     if(kminima(1)==1)then
        do k=1,ndom
           if(modulo(k,2).ne.0)then ! odd k: ascending
              thisk= k/2 +1
              domains(k,1)=kminima(thisk)
              domains(k,2)=kmaxima(thisk)
           else ! even k: descending
              thisk= k/2 
              domains(k,1)=kmaxima(thisk)
              thisk= k/2 +1
              domains(k,2)=kminima(thisk)
           end if
        end do
     else
        do k=1,ndom
           if(modulo(k,2).ne.0)then ! odd k: descending
              thisk= k/2 +1
              domains(k,1)=kmaxima(thisk)
              domains(k,2)=kminima(thisk)
           else ! even k: ascending
              thisk= k/2 
              domains(k,1)=kminima(thisk)
              thisk= k/2 +1
              domains(k,2)=kmaxima(thisk)
           end if
        end do
     end if
     

     do k = 1, ndom
        if(x(domains(k,1))<x(domains(k,2)))then
           klo=domains(k,1)
           khi=domains(k,2)
           direction=1
        else
           klo=domains(k,2)
           khi=domains(k,1)
           direction=-1
        end if

        call aspline(x(klo:khi:direction),y(:,klo:khi:direction),ddy(:,klo:khi:direction))

     end do
    
!     deallocate(domains)
    deallocate(kminima,kmaxima,domains)

   end subroutine maspline
   

   subroutine MultiValuedInvert(thisvpnr,z,tr)
     type(vpnr_type) :: thisvpnr
     real(dl) :: z
     integer :: mystat
     real(dl), pointer :: tr(:,:)
     real(dl), pointer :: mytr(:,:)
     integer :: k, ndom, thisk, rescount, thisdomain(2)
     integer, allocatable :: domains(:,:)
     real(dl) :: thesebounds(2), thistr(2)

     if(thisvpnr%monotonic_in_z)then
        if(size(tr(1,:))/=1)then
           deallocate(tr)
           allocate(tr(2,1),stat=mystat)
        end if
        call asplint(thisvpnr%ztr(1,:),thisvpnr%ztr(2:3,:),z,tr(:,1),thisvpnr%ddztr(2:3,:))
        return
     end if


     ! else:
     ! when redshift is not monotonic, the second derivatives of t and r
     ! with respect to redshift z are not correct (they can blow up).
     ! So we have to to an inversion...
     
     ndom = thisvpnr%nmaxz + thisvpnr%nminz - 1

     allocate(domains(ndom,2))
     if(thisvpnr%kminz(1)==1)then
        do k=1,ndom
           if(modulo(k,2).ne.0)then ! odd k: ascending
              thisk= k/2 +1
              domains(k,1)=thisvpnr%kminz(thisk)
              domains(k,2)=thisvpnr%kmaxz(thisk)
           else ! even k: descending
              thisk= k/2 
              domains(k,1)=thisvpnr%kmaxz(thisk)
              thisk= k/2 +1
              domains(k,2)=thisvpnr%kminz(thisk)
           end if
        end do
     else
        do k=1,ndom
           if(modulo(k,2).ne.0)then ! odd k: descending
              thisk= k/2 +1
              domains(k,1)=thisvpnr%kmaxz(thisk)
              domains(k,2)=thisvpnr%kminz(thisk)
           else ! even k: ascending
              thisk= k/2 
              domains(k,1)=thisvpnr%kminz(thisk)
              thisk= k/2 +1
              domains(k,2)=thisvpnr%kmaxz(thisk)
           end if
        end do
     end if

     rescount = 0
     do k = 1, ndom
        thesebounds(:) = thisvpnr%rtz(3,domains(k,:))
        thisdomain = (/minval(domains(k,:)),maxval(domains(k,:))/)
        if((z<=maxval(thesebounds)).and.(z>=minval(thesebounds)))then
           rescount = rescount + 1
           if(rescount>1)then
              allocate(mytr(2,size(tr(1,:))))
              mytr=tr
              deallocate(tr)
              allocate(tr(2,1+size(mytr(1,:))))
              tr(:,1:size(mytr(1,:)))=mytr
              deallocate(mytr)
           end if
           call ztrsplinvert(thisvpnr%rtz(1,thisdomain(1):thisdomain(2)),thisvpnr%rtz(2:3,thisdomain(1):thisdomain(2)),z,tr(:,rescount),thisvpnr%ddrtz(2:3,thisdomain(1):thisdomain(2)))
        end if
     end do


     deallocate(domains)


   end subroutine MultiValuedInvert


   subroutine ztrsplinvert(x,y,xx,yy,ddy)
     real(dl), intent(in) :: x(:),y(:,:),xx, ddy(:,:)
     real(dl), intent(out) :: yy(:)
     !- -
     real(dl) :: thisr, rtop, rbot
     real(dl) :: thisz
     real(dl) :: thist
     real(dl) :: thistz(2), thistr(2)
     integer :: direction, counter, maxit
     ! goal:
     thisz = xx
     ! this domain is monotonic:
     rtop=maxval(x)
     rbot=minval(x)
     direction=int(sign(1.d0,y(2,size(y(2,:)))-y(2,1)))

     thisr=(rtop+rbot)*0.5_dl
     maxit=60 ! 2^60 = 1e18
     counter = 0
     do while(dabs(rtop-rbot)>1.e-15_dl)
        call asplint(x,y,thisr,thistz,ddy)
        if(dble(direction)*(thistz(2)-xx)>0._dl)then
           rtop=thisr
        else
           rbot=thisr
        end if
        thisr=0.5_dl*(rtop+rbot)
        counter = counter + 1
        if(counter > maxit)exit ! because then the domain already shrank by 1e-18, so
                                ! the answer won't get any better. Good enough.
     end do

     thistr(1)=thistz(1)
     thistr(2)=thisr
     yy=thistr

   end subroutine ztrsplinvert

   subroutine deallocateThisVPNR(tVPNR)
    implicit none
    type(vpnr_type) :: tVPNR
    integer :: mystat
   
    if(associated(tVPNR%rtz))deallocate(tVPNR%rtz,stat=mystat)
    if(associated(tVPNR%ddrtz))deallocate(tVPNR%ddrtz,stat=mystat)

    if(associated(tVPNR%ztr))deallocate(tVPNR%ztr,stat=mystat)
    if(associated(tVPNR%ddztr))deallocate(tVPNR%ddztr,stat=mystat)
    if(associated(tVPNR%ztrerr))deallocate(tVPNR%ztrerr,stat=mystat)
    if(associated(tVPNR%kmaxz))deallocate(tVPNR%kmaxz,stat=mystat)
    if(associated(tVPNR%kminz))deallocate(tVPNR%kminz,stat=mystat)
  
  end subroutine deallocateThisVPNR


  subroutine copyVPNRToThisVPNR(sVPNR,tVPNR)
    implicit none
    type(vpnr_type) :: sVPNR,tVPNR ! source, target
    
    allocate(tVPNR%rtz(size(sVPNR%ztr(:,1)),size(sVPNR%ztr(1,:))))
    allocate(tVPNR%ddrtz(size(sVPNR%ddztr(:,1)),size(sVPNR%ddztr(1,:))))

    allocate(tVPNR%ztr(size(sVPNR%ztr(:,1)),size(sVPNR%ztr(1,:))))
    allocate(tVPNR%ddztr(size(sVPNR%ddztr(:,1)),size(sVPNR%ddztr(1,:))))
    allocate(tVPNR%ztrerr(size(sVPNR%ztr(:,1)),size(sVPNR%ztr(1,:))))
    allocate(tVPNR%kmaxz(sVPNR%nmaxz))
    allocate(tVPNR%kminz(sVPNR%nminz))
    
    tVPNR%rtz = sVPNR%rtz
    tVPNR%ddrtz = sVPNR%ddrtz

    tVPNR%ztr = sVPNR%ztr
    tVPNR%ddztr = sVPNR%ddztr
    tVPNR%ztrerr = sVPNR%ztrerr
    tVPNR%kmaxz = sVPNR%kmaxz
    tVPNR%kminz = sVPNR%kminz
    tVPNR%nmaxz = sVPNR%nmaxz
    tVPNR%nminz = sVPNR%nminz

    tVPNR%monotonic_in_z = sVPNR%monotonic_in_z

  end subroutine copyVPNRToThisVPNR

  subroutine associateVPNRToThisVPNR(sVPNR,tVPNR)
    implicit none
    type(vpnr_type) :: sVPNR,tVPNR ! source, target
    
    nullify(tVPNR%rtz)
    nullify(tVPNR%ddrtz)

    nullify(tVPNR%ztr)
    nullify(tVPNR%ddztr)
    nullify(tVPNR%ztrerr)
    nullify(tVPNR%kmaxz)
    nullify(tVPNR%kminz)
    
    tVPNR%rtz => sVPNR%rtz
    tVPNR%ddrtz => sVPNR%ddrtz

    tVPNR%ztr => sVPNR%ztr
    tVPNR%ddztr => sVPNR%ddztr
    tVPNR%ztrerr => sVPNR%ztrerr
    tVPNR%kmaxz => sVPNR%kmaxz
    tVPNR%kminz => sVPNR%kminz
    tVPNR%nmaxz = sVPNR%nmaxz
    tVPNR%nminz = sVPNR%nminz

    tVPNR%monotonic_in_z = sVPNR%monotonic_in_z

  end subroutine associateVPNRToThisVPNR

  subroutine nullifyThisVPNR(tVPNR)
    implicit none
    type(vpnr_type) :: tVPNR ! source, target
    
    nullify(tVPNR%rtz)
    nullify(tVPNR%ddrtz)

    nullify(tVPNR%ztr)
    nullify(tVPNR%ddztr)
    nullify(tVPNR%ztrerr)
    nullify(tVPNR%kmaxz)
    nullify(tVPNR%kminz)
    
  end subroutine nullifyThisVPNR

