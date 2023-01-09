   subroutine voidsnoutput(fnum, vdai, mylogh)
    use voiddistSNchi2
     use ode_path

!     use camb
     integer :: fnum
     real :: mylogh
     type(voiddistinfo) :: vdai
     real(dl) :: minz, maxz, testr, dtm_tr(2), thisz, thismu, thisda, rstart, maxr, minr
     integer :: i, testimax, j
     integer, parameter :: nres = 5
     real(dl), allocatable :: ddyy(:,:), testy(:)
     character(len=32) :: snfile(17)
     character(len=256) :: myline
     integer :: mystat, nmax

!     minz = minval(vdai%snz)
     minz=0._dl
     maxz = maxval(vdai%snz)
     !testimax = nres*kount
     testimax = 5000
!     if(dabs(VP%zB).gt.1.d-20)then
     if(.not.VP%is_pure_flrw)then
!        allocate(testy(size(yp(:,1))), ddyy(size(yp(:,1)),kount))
!        call aspline(xp(1:kount), yp(:,1:kount), ddyy)
        nmax=size(VPNR%rtz(:,1))
        allocate(testy(nmax-1))
        
        ! get rstart:
        rstart = NR(minz)
        minr = rstart!NR(minz)
        maxr = NR(maxz)

        do i = 1, testimax
!           testr = rstart + (xp(kount) - rstart - xp(1))/(dble(testimax-1))*dble(i-1)
           testr = (maxr-minr)*dble(i-1)/dble(testimax-1) + minr
!           call asplint(xp(1:kount),yp(:,1:kount),testr,testy,ddyy)
           call asplint(VPNR%rtz(1,:),VPNR%rtz(2:nmax,:),testr,testy,VPNR%ddrtz(2:nmax,:))
           
           dtm_tr(1) = testy(1)
           dtm_tr(2) = testr
           
           
           thisz = testy(2) !-testy(2)-1.d0
!           if(thisz.lt.minz)cycle
!           if(thisz.gt.maxz)exit
           
           thisda =  DA(0.d0, dtm_tr)
           !          thismu=5.d0*log10((1.d0+thisz)**2*thisda)+25.d0 
           
           !          write(fnum,'(7E14.5)')-testy(2)-1.d0,thisda,thismu,thismu+mylogh
           call writeline()
        end do

!!$        ! Add last point:
!!$        testr = NR(maxz)
!!$        call asplint(xp(1:kount),yp(:,1:kount),testr,testy,ddyy)
!!$        dtm_tr(1) = testy(1)
!!$        dtm_tr(2) = testr
!!$        thisda =  DA(0.d0, dtm_tr)
!!$        call writeline()

        
        deallocate(testy)!,ddyy)
        
        
     else
        allocate(testy(2))
        do i = 1, testimax
           thisz = minz + (maxz-minz)/(dble(testimax-1))*dble(i-1)
           testy(2) = -(1.d0+thisz)
           call writeline()
        end do
        deallocate(testy)

        mystat = 0
        open(133,file='../data_share/MLCS.FITRES',status='OLD')
        open(132,file='renormmoduli.txt', status='REPLACE')
        do i = 1, size(vdai%moduli)+20
           read(133,*,iostat = mystat)snfile
           if(mystat/=0)exit
           if(trim(snfile(1))=='')cycle
!           if(i<3)cycle
           snfile(3)=adjustl(snfile(3))
           if(ichar(snfile(3)(1:1))<ichar('0') .or. ichar(snfile(3)(1:1))>ichar('9'))cycle 
           read(snfile(3),*)thisz
!           thisz = vdai%snz(i)
!           thisda = DA(thisz)
           thismu=5.d0*log10((1.d0+thisz)**2*thisda)+25.d0 +mylogh
           read(snfile(5),*)thisz
           thisz = thisz - thismu
           myline = ''
           do j = 1, size(snfile)
              myline = trim(myline)//' '//trim(snfile(j))
           end do
           write(132,'(A,2E15.4)')trim(myline),thisz,thismu
 !          write(132,'(7E14.5)')thisz,thismu,vdai%moduli(i),vdai%moduli(i)-thismu
           
        end do
        close(132)
        close(133)
     end if

     return
     
     contains

       subroutine writeline()

          thismu=5.d0*log10((1.d0+thisz)**2*thisda)+25.d0 

!          write(fnum,'(7E14.5)')-testy(2)-1.d0,thisda,thismu,thismu+mylogh
          write(fnum,'(7E14.5)')thisz,thismu,thismu+mylogh

       end subroutine writeline

   end subroutine voidsnoutput


   subroutine VoidSNMock()
     implicit none
!     character(len=10), parameter :: SNSET="SDSSMLCS"
!     character(len=10), parameter :: SNSET="SDSSSALT"
!     character(len=10), parameter :: SNSET="Union"
     character(len=10), parameter :: SNSET(3)=(/"Union     ","SDSSSALT  ","SDSSMLCS  "/)
     character(len=1), parameter :: whitespaces(2)=(/' ',achar(9)/)
     character(len=256) :: inputfile
     character(len=256) :: outputfile, outputfilesubscript
     character(len=256) :: inline, thiszstr, thismustr
 !    integer, parameter :: mucoll(6) = (/14, 16, 18, 20, 22, 24/)
     integer, allocatable :: mucoll(:)
     integer :: zcoll
     integer :: i, n, fnum, fnum2, mystat, ii
     real(dl) :: thisz, thismu, thisda
     
     do ii = 1, size(SNSET)
        voiderror=.false.
        
        call dothemagic
        
        ! old:     write(outputfilesubscript,'(F8.6)')abs(1._dl-VP%coldspot)
        write(outputfilesubscript,'(F8.3)')VP%coldspotmk
        outputfilesubscript=trim(adjustl(outputfilesubscript))
        
        if(trim(SNSET(ii))=="SDSSMLCS")then
           inputfile="/home/wv347191/cosmomc/data_share/SDSS_SN_2009/MLCS.FITRES"
           outputfile="/home/wv347191/cosmomc/data_share/SDSS_SN_2009/MYMOCKMLCS_"//trim(outputfilesubscript) // ".FITRES"
           allocate(mucoll(1))
           mucoll = (/5/)
           zcoll=3
        else if(trim(SNSET(ii))=="SDSSSALT")then
           inputfile="/home/wv347191/cosmomc/data_share/SDSS_SN_2009/SALT2.FITRES"
           outputfile="/home/wv347191/cosmomc/data_share/SDSS_SN_2009/MYMOCKSALT_"//trim(outputfilesubscript) // ".FITRES"
           allocate(mucoll(6))
           mucoll = (/14, 16, 18, 20, 22, 24/)
           zcoll=3
        else if (trim(SNSET(ii))=="Union")then
           inputfile="/home/wv347191/cosmomc/data_share/Union2/sn_z_mu_dmu_union2.txt"
           outputfile="/home/wv347191/cosmomc/data_share/Union2/MYMOCKUnion2_"//trim(outputfilesubscript) // ".txt"
           allocate(mucoll(1))
           mucoll = (/3/)
           zcoll=2
        else
           stop 'SN Undefined in mock sn.'
        end if
        write(*,'(A)')' Making mock SN based on '//trim(SNSET(ii))
        write(*,'(A)')' Reading SN from '//trim(inputfile)
        writE(*,'(A)')' Writing mock SN to '//trim(outputfile)
        
        
        fnum=MyNewIONum()
        open(fnum,file=inputfile,status='old')!,readonly)
        fnum2=MyNewIONum()
        open(fnum2,file=outputfile,status='replace')
        
        mystat=0
        do while (mystat==0)
           read(fnum,'(A)',iostat=mystat)inline
           if(mystat/=0)exit
           inline=adjustl(inline)
           if((inline(1:3)=='SN:').or.(trim(SNSET(ii))=='Union'))then ! generate our mock data
              thiszstr = MyGetColl(inline,zcoll)
              read(thiszstr,*)thisz
              thisda = DA(thisz)
              thismu=5.d0*log10((1.d0+thisz)**2*thisda)+25.d0 
              write(thismustr,'(F18.9)')thismu
              !           thismustr=adjustl(thismustr)
              call MyReplaceColls(inline,thismustr,mucoll)
           end if
           write(fnum2,'(A)')trim(adjustl(inline))
        end do
        
        deallocate(mucoll)
        
        close(fnum)
        close(fnum2)
        
     end do
     
   contains
     
     subroutine MyReplaceColls(string,newstr,colns)
       implicit none
       character(len=*), intent(inout) :: string
       character(len=*), intent(in) :: newstr
       integer, intent(in) :: colns(:)
       integer :: i, n
       
       n=size(colns)
       do i = 1, n
          call MyReplaceColl(string,newstr,colns(i))
       end do
       
     end subroutine MyReplaceColls
     
     subroutine MyReplaceColl(string,newstr,coln)
       implicit none
       character(len=*), intent(inout) :: string
       character(len=*), intent(in) :: newstr
       character(len=256) :: part1, part2
       integer, intent(in) :: coln
       integer :: i, n, poss(2)
       character(len=256) :: dummy
       
       dummy = MyGetColl(string,coln,poss)
       part1 = string(1:poss(1)-1)
       part2 = string(poss(2)+1:)
       string = trim(part1)//trim(newstr)//trim(part2)
       
     end subroutine MyReplaceColl
     
     function MyGetColl(string,coln, poss)
       implicit none
       character(len=*), intent(in) :: string
       integer, intent(in) :: coln
       integer, intent(out), optional :: poss(2)
       character(len=256) :: MyGetColl, mystring
       integer :: mi, mn, imem1, imem2, lastcharpos, itemcount
       
       mystring = adjustl(string)
       mn = len(trim(mystring))
       
       imem1=1
       imem2=0
       lastcharpos=1
       itemcount = 0
       do i = 1, mn
          if(all(mystring(i:i)/=whitespaces))then
             if(lastcharpos < i-1)then ! previous was space
                itemcount=itemcount + 1
                if(itemcount==coln)then
                   imem2=i-1 ! = ' '
                   exit
                else
                   imem1=i
                end if
             end if
             lastcharpos=i
          end if
       end do
       
       MyGetColl = trim(mystring(imem1:imem2))
       !imem2 = imem1 + len(trim(MyGetColl))-1
       imem2 = imem1 + mytrimlen(MyGetColl)-1
       MyGetColl = mystring(imem1:imem2)
       if(present(poss))poss=(/imem1,imem2/)
       
     end function MyGetColl
     
     function mytrimlen(strin)
       implicit none
       integer :: mytrimlen
       character(len=*) :: strin
       integer :: ii
       
       mytrimlen=len(trim(strin))
       do ii = len(trim(strin)), 1, -1
          mytrimlen=ii
          if(all(strin(ii:ii)/=whitespaces))exit
       end do
       
     end function mytrimlen
     
     function MyNewIONum()
       integer :: MyNewIONum, i
       logical :: notfree
       integer, parameter :: maxit = 10000
       
       notfree = .true.
       i=20
       
       do while (notfree)
          i=i+1
          inquire(i,opened=notfree)
          if(i.gt.maxit)stop 'MyNewIONum could not find free filenumber.'
       end do
       MyNewIONum = i
       
     endfunction MYNewIONum
   end subroutine VoidSNMock
   
   
   Subroutine MyReplaceChars(string,in,out)
     implicit none
     character(len=*) :: string, in, out
     integer :: i, n, score, m
     integer, allocatable :: poslist(:), poslist2(:)
     
     n = len(trim(string))
     m = len(in)
     allocate(poslist(1))
     score = 0
     
     ! Count items:
     do i = 1, n
        if(string(i:i+m-1)==in)then
           score = score+1
           if(score>size(poslist))then
              allocate(poslist2(score))
              poslist2(1:score-1)=poslist
              deallocate(poslist)
              allocate(poslist(score))
              poslist = poslist2
              deallocate(poslist2)
           end if
           poslist(score)=i       
        end if
     end do
     
     if(score>0)then
        do i = 1, score
           !         write(*,*)trim(string),poslist
           call MyReplaceString(string,in,out,poslist(i))
           poslist = poslist - m + len(out)
        end do
        
     end if
     
     
   end Subroutine MyReplaceChars
   
   subroutine MyReplaceString(string,in,out,pos)
     implicit none
     character(len=*) :: string,in,out
     integer :: pos
     character(len=len(string)) :: part1, part2
     
     part1 = string(1:pos-1)
     part2 = string(pos+len(in):len(string))
     
     string = trim(part1)//out//part2
     
     
   end subroutine MyReplaceString
   
