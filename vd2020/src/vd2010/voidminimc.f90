subroutine vdi_minimc(p,ranges,proposals,temperature,chi2)
  implicit none
  real(dl), intent(inout) :: p(:)
  real(dl), intent(in) :: ranges(:),proposals(:)
  real(dl), intent(in) :: temperature
	INTERFACE
		FUNCTION chi2(pin)
			use precision
		IMPLICIT NONE
		REAL(dl), DIMENSION(:), INTENT(IN) :: pin
		REAL(dl) :: chi2
		END FUNCTION chi2
	END INTERFACE
  integer, parameter :: maxsamples=1000, maxtries = 100, updateprop=25!maxsamples/10
  real(dL) :: chain(size(p)+2,maxsamples), currsamp(size(p)+2), nextsamp(size(p)+2)
  real(dl) :: propvecs(size(p),size(p)), evals(size(p)), proposesigmas(size(p)), startp(size(p))
  integer :: i, n, j, usesamples, k, v
  real(dl) :: nextgauss(size(p)), sddevs(size(p)), totsamps, xi(size(p),size(p)), meanchi, chivar

  n=size(p)
  startp = p
  propvecs=0.d0
  do i = 1, n
    propvecs(i,i) = 1.d0!proposals(i)!0.01*(ranges(2:2*n:2)-ranges(1:2*n-1:2)) ! fifth of the prior width
    evals(i) = proposals(i)**2
    proposesigmas(i) = proposals(i)
  end do

  ! start with rough guess, just off-center:
 ! p =ranges(1:2*n-1:2) + 0.4 * (ranges(2:2*n:2)-ranges(1:2*n-1:2))
  ! let input give best guess.
   xi=0.d0
  do i = 1, n
    xi(i,i) = 1.d0!p(i)*0.01
  end do
  call void_nr_powell(p,xi,1.d-4,i,currsamp(2),chi2)

  currsamp(3:) = p
  currsamp(2) = chi2(p) / temperature
  currsamp(1) = 1.d0
  do i = 1, maxsamples
   do k = 1, maxtries
    do j = 1, maxtries
      call vdi_NextGauss(nextgauss)
      nextgauss = nextgauss * proposesigmas
      nextsamp(3:) = currsamp(3:) +  matmul(propvecs,nextgauss)
      nextsamp(3:) = min(nextsamp(3:),ranges(2:2*n:2))
      nextsamp(3:) = max(nextsamp(3:),ranges(1:2*n-1:2))
      nextsamp(2) = chi2(nextsamp(3:)) / temperature
      nextsamp(1) = 1.d0
!write(*,*)'delta:',matmul(propvecs,nextgauss)
!write(*,*)'next:',nextsamp
!write(*,*)'curr:',currsamp
!write(*,*)'prop:',proposesigmas,i
      if(vdi_accept(currsamp(2),nextsamp(2)))then
        chain(:,i) = currsamp
        exit
      else
        currsamp(1) = currsamp(1) + 1
      end if
    end do

    if(j>=maxtries .and. .not. all(chain(:,i) == currsamp))then ! not the case if we accepted
      ! smaller steps
      proposesigmas = proposesigmas * 0.5
      proposesigmas = max(proposesigmas,1.d-7*proposals)
!if(i>1)write(*,*)'prop:',proposesigmas,i,nint(currsamp(1))
    else
      chain(:,i) = currsamp
      currsamp = nextsamp
      exit
    end if
   end do

!    if(modulo(i,updateprop)==0)then
    if(i==updateprop .or. (modulo(i,updateprop)==0 .and. any(abs(proposesigmas**2/abs(evals)-1.d0)>1.d-14)))then
      call vdi_UpdateProposal(chain(:,i-updateprop+1:i),propvecs,evals,size(p))
      proposesigmas = max(1.d-7*proposals, sqrt(abs(evals))) ! abs because -0.d0 may occur.
!      write(*,*)proposesigmas
    end if

!
!      totsamps = sum(chain(1,:i))
!      meanchi = sum(chain(1,:i)*chain(2,:i))/totsamps
!      chivar = sqrt(sum(chain(1,:i)*(meanchi-chain(2,:i))**2)/totsamps)
!      do j = 1, n
!        p(j) = sum(chain(1,:i)*chain(2+j,:i))/totsamps
!        sddevs(j) = sqrt(sum(chain(1,:i)*(p(j)-chain(2+j,:i))**2)/totsamps)
!      end do
!      write(*,'(I6,20ES14.5)')i,chivar,p,sddevs,currsamp(2+2)
!write(*,'(I6,20ES14.5)')i,chain(:,i)
  end do

  usesamples = maxsamples/2
  totsamps = sum(chain(1,usesamples:))
  do i = 1, n
    p(i) = sum(chain(1,usesamples:)*chain(2+i,usesamples:))/totsamps
!write(*,*)i,sum(chain(1,usesamples:)*chain(2+i,usesamples:))
    sddevs(i) = sqrt(sum(chain(1,usesamples:)*(p(i)-chain(2+i,usesamples:))**2)/totsamps)
  end do

!  write(*,*)'mcmc'
!  write(*,*)'p:',p
!  write(*,*)'sig:',sddevs
!  write(*,*)'ndev of start:',(startp-p)/sddevs
!  write(*,*)'totsamps:',totsamps
!stop

contains

  function vdi_accept(curr,next)
    IMPLICIT NONE
    real(dl), intent(in) :: curr, next
    logical :: vdi_accept
    real(dl) :: thisrnd, p
    if(next<curr)then
      vdi_accept = .true.
    else
      call random_number(thisrnd)
      p = exp(-next*0.5+curr*0.5)
      if(p>thisrnd)then
        vdi_accept = .true.
      else
        vdi_accept = .false.
      end if
    end if

  end function vdi_accept

  subroutine vdi_NextGauss(ng)
    implicit NONE
    real(dl), intent(out) :: ng(:)
    integer :: n, halfn,np, ii
    real(dL) :: u,v,g1,g2

    n = size(ng)
    halfn = n/2
    np=n
    if(n/=2*halfn)np=n+1
    halfn = np/2

    ! http://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
    do ii = 1, halfn
     call random_number(u)
     call random_number(v)
      g1 = sqrt(-2 * log(u)) * cos(2*pi*v)
      g2 = sqrt(-2 * log(u)) * sin(2*pi*v)
      ng(2*ii-1) = g1
      if(2*ii<n) ng(2*ii) = g2
    end do

  end subroutine vdi_NextGauss

  subroutine vdi_UpdateProposal(chainin,outmat,mevals,nn)
    implicit NONE
    real(dl), intent(in) :: chainin(:,:)
    integer, intent(in) :: nn
    real(dl), intent(out) :: outmat(nn,nn),mevals(nn)
    integer :: ii, jj
    real(dl) :: nsamp, inmat(nn,nn)

    nsamp = sum(chainin(1,:))
    do ii = 1, nn
      do jj = ii, nn
        outmat(ii,jj) = (sum(chainin(2+ii,:)*chainin(2+jj,:))/nsamp)
        outmat(jj,ii) = outmat(ii,jj)
      end do
    end do

!    write(*,'(3ES14.5)')outmat
    inmat = outmat
    call Matrix_Diagonalize(outmat,mevals,nn)
    return

write(*,*)
    write(*,'(3ES14.5)')outmat
write(*,*)
    write(*,'(3ES14.5)')mevals
do ii = 1, nn
write(*,*)
    write(*,'(3ES14.5)')matmul(inmat,outmat(:,ii))/outmat(:,ii),sum(outmat(:,ii)**2)
    write(*,'(3ES14.5)')matmul(inmat,outmat(ii,:))/outmat(ii,:),sum(outmat(ii,:)**2)
end do
stop

  end subroutine vdi_UpdateProposal

  ! literal copy form cosmomc
  subroutine Matrix_Diagonalize(M, diag, n)
  !Does m = U diag U^T, returning U in M
      integer, intent(in) :: n
      real(dl), intent(inout):: m(n,n)
      real(dl), intent(out) :: diag(n)
      integer ierr, tmpsize
      real(dl), allocatable, dimension(:) :: tmp
      integer, external :: ILAENV
 
      tmpsize =  max( (ILAENV(1,'DSYTRD','U',n,n,n,n)+2)*N,max(1,3*n-1))  !3*n**2
      allocate(tmp(tmpsize));
      call DSYEV('V','U',n,m,n,diag,tmp,tmpsize,ierr) !evalues and vectors of symmetric matrix

      if (ierr /= 0) call MpiStop('Error in Matrix_Diagonalize')
      deallocate(tmp)

  end subroutine Matrix_Diagonalize

end subroutine vdi_minimc


