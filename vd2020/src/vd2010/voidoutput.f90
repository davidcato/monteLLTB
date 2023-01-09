
subroutine VoidIniBestFitOutput(rootname)
  implicit none
  character(len=*), intent(in) :: rootname

  open(222,file=trim(rootname)//'best_kofr.txt', status='replace')
!  open(123,file=trim(rootname)//'best_mpk.txt', status='replace')
  open(124,file=trim(rootname)//'best_cls.txt', status='replace')
  open(125,file=trim(rootname)//'best_sn.txt', status='replace')
  open(126,file=trim(rootname)//'best_bao.txt', status='replace')
  


end subroutine VoidIniBestFitOutput

subroutine VoidFinishBestFitOutput()
  implicit none
  type(VoidLikesMem), pointer :: myVLM(:)
  integer :: i, n

  close(222)
  close(124)
  close(125)
  close(126)

  call VoidSaveLikes(VLMpoint=myVLM)
  n=size(myVLM)
  do i = 1, n
     write(*,*)trim(myVLM(i)%name),myVLM(i)%like
  end do
  stop 'Bestfit output done.'

end subroutine VoidFinishBestFitOutput

subroutine VoidBestFitOutputCMB(acl,cmblike)
  implicit none
  real, intent(in) :: acl(:,:), cmblike
  integer :: cols
  integer :: ls
  character(len=64) :: myformat, intchar
  integer :: i
 
  cols = size(acl(1,:))
  ls = size(acl(:,1))

  write(intchar,*)cols
  write(myformat,*)'(I5,'//trim(adjustl(intchar))//'E14.5)'
  
  do i = 2, size(acl(:,1))
     write(124,myformat)i,acl(i,:)*i*(i+1)*0.5d0/pi
  end do

  call VoidSaveLikes("CMB",cmblike)
  
end subroutine VoidBestFitOutputCMB

subroutine VoidBestFitOutputSN() !vdai,snlike)
    write(*, *) "VoidBestFitOutputSN does nothing."
    return
!  use voiddistsnchi2
!  implicit none
!  real, intent(in) :: snlike
!!  real(dl), intent(in) :: logh
!  type(voiddistinfo), intent(in) :: vdai
!
!  call voidsnoutput(125,vdai,-real(vdai%avvdiff))
!
!  call VoidSaveLikes("SN",snlike)
!
end subroutine VoidBestFitOutputSN


subroutine VoidBestFitOutputHST(theory,obs,sigma,hstlike)
  implicit none
  real(dl), intent(in) :: theory,obs,sigma,hstlike
  
  call VoidSaveLikes("HST",real(hstlike))

  write(*,*)'HST:'
  write(*,*)'Observed distance:',1.d0/obs
  write(*,*)'Theoretical distance:',1.d0/theory
  write(*,*)'sigma:',1.d0/obs**2 * sigma

end subroutine VoidBestFitOutputHST

subroutine VoidSaveLikes(name,like,VLMpoint)
  implicit none
  character(len=*), intent(in), optional :: name
  real, intent(in), optional :: like
  type(VoidLikesMem), pointer, optional :: VLMpoint(:)
  type(VoidLikesMem), pointer, save :: VLM(:)
  type(VoidLikesMem), pointer :: VLMtmp(:)
  logical, save :: firstcall = .true.
  integer :: inisize, endsize

  if(present(VLMpoint))then
     if(firstcall)stop 'Attempting to retrieve saved lnlikes, but no likes were saved.'
     VLMpoint => VLM
     return
  end if

  if(.not.(present(name).and.present(like)))stop 'Calling VoidSaveLikes with either name or like missing.'

  if(.not.firstcall)then
     inisize = size(VLM)
     endsize = inisize + 1
     allocate(VLMtmp(inisize))
     VLMtmp = VLM
     deallocate(VLM)
     allocate(VLM(endsize))
     VLM(1:inisize) = VLMtmp
     deallocate(VLMtmp)
  else
     firstcall = .false.
     inisize = 0
     endsize = 1
     allocate(VLM(endsize))
  end if

  VLM(endsize)%name = name
  VLM(endsize)%like = like

end subroutine VoidSaveLikes
