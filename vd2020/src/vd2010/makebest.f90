program MakeBestfitIni
  implicit none 
  integer, parameter :: as = 50
  integer, parameter :: mylen = 256
  character(len=mylen) :: inifile, bfinifile, fileroot, fullfileroot
!  character(len=7) :: param(as)
  character(len=mylen), allocatable :: likestats(:)
  character(len=mylen), allocatable :: otherpars(:)
  character(len=mylen), pointer :: namesarray(:)

!  call prepparams(param)

  call GetFlag(inifile)

  call GetFileRoot(inifile, fileroot, fullfileroot)
  
  call MakeBFname(inifile,bfinifile)

  call InquireParamNames(inifile, namesarray)

  call ReadLikestats(fileroot, likestats, namesarray)

  call Setotherpars(otherpars,fullfileroot)

  call MakeBF(inifile,bfinifile,likestats,otherpars)

  call DoRun(bfinifile)

  deallocate(likestats,otherpars)

contains
  function NewIONum()
    integer :: NewIONum, i
    logical :: notfree
    integer, parameter :: maxit = 10000

    notfree = .true.
    i=20
    
    do while (notfree)
       i=i+1
       inquire(i,opened=notfree)
       if(i.gt.maxit)stop 'NewIONum could not find free filenumber.'
    end do
    NewIONum = i
    
  endfunction NewIONum

  subroutine gettempfile(thisname)
    character(len=*) :: thisname
    logical :: exists
    integer :: i
    integer, parameter :: maxit = 10000
    character(len=12) :: thisnum
    character(len=12), parameter :: theroot="getflgtmp"

    exists=.true.
    i=0
    do while (exists)
       write(thisnum,*)i
       thisname = trim(theroot)//trim(adjustl(thisnum))
       inquire(file=thisname, exist=exists)
       if(.not.exists)exit
       i=i+1
    end do
    
    call system('touch '//trim(thisname))

  end subroutine gettempfile

  subroutine InquireParamNames(ifile,narray)
    character(len=*) :: ifile
    character(len=*), pointer :: narray(:)
    character(len=mylen) :: comm, fname, thisline
    integer :: fnum, mystat, i, flen, j, k, m

    call GetTempFile(fname)

    comm = 'grep "param\[" '//trim(ifile)//' > '//trim(fname)
    call system(comm)

    fnum = NewIOnum()
    
    open(fnum,file=fname)
    mystat = 0
    i = 0
    do while (mystat==0)
       read(fnum,'(A)',iostat=mystat)thisline
       if(mystat/=0)exit
       i=i+1
    end do
   
    close(fnum)
    comm = 'rm -f '//trim(fname)
    call system(comm)


    if (i>0)then ! names
       open(fnum,file="params_CMB.paramnames")
       mystat=0
       i=0
       do while (mystat==0)
          read(fnum,'(A)',iostat=mystat)thisline
          if(mystat/=0)exit
          i=i+1
       end do
       close(fnum)
       flen=i

       allocate(narray(flen))
       narray = 'XXXX'
       
       open(fnum,file="params_CMB.paramnames")
       mystat=0
       i=0
       m=0
       do k = 1, flen
          i=k-m
!write(*,*)flen,i,k,m
          read(fnum,'(A)')narray(i)
          narray(i) = adjustl(narray(i))
          if(narray(i)(1:1)=='#')then
             m=m+1
             cycle
          end if
          do j = 1, len(trim(narray(i)))
             if(narray(i)(j:j)==' ' .or. narray(i)(j:j)==char(9))exit             
          end do
          if(j/=len(trim(narray(i))))j=j-1
          narray(i)=narray(i)(1:j)//'                                 '
          narray(i) = 'param['//trim(adjustl(narray(i)))//']'
       end do

       close(fnum)
    else ! numbers
       allocate(narray(100))
       do i = 1, 100
          write(narray(i),*)i
          narray(i) = 'param'//trim(adjustl(narray(i)))
       end do
    end if
       

  end subroutine InquireParamNames

 
  subroutine DoRun(ifile)
    implicit none
    character(len=*)ifile
    character(len=mylen) :: comm

    comm = './cosmomc '//trim(adjustl(ifile))
    
    call system(comm)
    
  end subroutine DoRun

  subroutine MakeBF(ifile,bfifile,stats,opars)
    implicit none
    character(len=*) :: ifile, bfifile, stats(:), opars(:)
    character(len=mylen) :: thisline
    character(len=mylen), allocatable :: alllines(:)    
    integer :: fnum(2) = (/125,126/), mystat, i
    logical :: donepar(size(opars)), isthere
    inquire(file=ifile,exist=isthere)
    if(.not.isthere)then
       write(*,*)'File '//trim(ifile)//' not found.'
       stop
    end if

    open(fnum(1),file=ifile,status='OLD')
    open(fnum(2),file=bfifile,status='REPLACE')

    allocate(alllines(size(opars)+size(stats)))
    alllines(1:size(opars))=opars
    alllines(size(opars)+1:size(opars)+size(stats))=stats
    alllines = adjustl(alllines)

    mystat = 0
    donepar = .false.
    do while (mystat==0)
       read(fnum(1),'(A)',iostat = mystat)thisline
       if(mystat/=0)exit
       
       do i = 1, size(alllines)
          ! check first three characters for simplicity        
          if(thisline(1:3)==alllines(i)(1:3)) call ReplaceIfSame(thisline,alllines(i))
          if(thisline==alllines(i) .and. i.le.size(opars))donepar(i) = .true.
       end do
       write(fnum(2),'(A)')trim(thisline)
    end do
    do i = 1, size(opars)
       if(.not.donepar(i))write(fnum(2),'(A)')trim(opars(i))
    end do

    close(fnum(1))
    close(fnum(2))

    deallocate(alllines)
    
  end subroutine MakeBF
  
  subroutine ReplaceIfSame(line,newline)
    implicit none
    character(len=*) :: line,newline
    character(len=mylen) :: mystr(4)

    mystr(1) = line
    call split(mystr(1),mystr(2),' ')
    mystr(3) = newline
    call split(mystr(3),mystr(4),' ')
    
    ! Now mystr(2) and mystr(4) contain first word of line
    if(trim(adjustl(mystr(2)))==trim(adjustl(mystr(4))))line = newline


  end subroutine ReplaceIfSame

  subroutine Setotherpars(pars,fullroot)
    implicit none
    character(len=*) :: fullroot
    character(len=*), allocatable :: pars(:)

    if(allocated(pars))deallocate(pars)

    allocate(pars(7))

    pars(1) = 'Do_bestfit_output = T'
    pars(2) = 'VoidTesting = F'
    pars(3) = 'VoidTestingIntegrationOutput = F'
    pars(4) = trim(adjustl(fullroot))//'_bf'
    pars(5) = 'checkpoint = F'
    pars(6) = 'propose_matrix = '
    pars(7) = 'temperature = 1'


  end subroutine Setotherpars

  subroutine ReadLikestats(root,stats, narray)
    implicit none
    character(len=*), allocatable :: stats(:)
    character(len=mylen), allocatable :: mystats(:)
    character(len=*) :: root
    character(len = len(root)+9) :: likefile
    character(len=mylen) :: myline(7), outline(5)
    character(len=1) :: thischar
    real(8) :: thisval(5)
    real(8), parameter :: myspace = 1.d-6
    integer :: fnum, mystat, thispar, npar, i
    logical :: isthere
    character(len=mylen), pointer :: narray(:)
    if(allocated(stats))deallocate(stats)

    likefile = trim(root)//'.likestats'
    fnum = 124
    inquire(file=likefile,exist=isthere)
    if(.not.isthere)then
       likefile = 'distout/'//trim(root)//'.likestats'
       inquire(file=likefile,exist=isthere)
       if(.not.isthere)then
          write(*,*)'File '//trim(likefile)//' not found.'
          stop
       end if
    end if

    open(fnum,file=likefile)
    mystat = 0

    allocate(mystats(as))
    mystats=''
    
    npar = 0

    do while (mystat==0)
       read(fnum,*,iostat = mystat)myline
       if(mystat/=0)exit
       thischar = trim(adjustl(myline(1)))
       if(ichar(thischar)>=ichar('0') .and. ichar(thischar)<=ichar('9'))then  
          ! We have a parameter line
          npar = npar + 1
          read(myline(1),*)thispar
          read(myline(2),*)thisval(1)
          if(thisval(1)>0)then
             thisval(2) = thisval(1)*(1.d0 - myspace)
             thisval(3) = thisval(1)*(1.d0 + myspace)
          else
             thisval(3) = thisval(1)*(1.d0 - myspace)
             thisval(2) = thisval(1)*(1.d0 + myspace)
          end if
          thisval(4) = 0.d0
          thisval(5) = dabs(thisval(1))*myspace * 1.d-3
          do i = 1, 5
             write(outline(i),'(ES19.10)')thisval(i)
          end do

!          mystats(npar)='param'//trim(adjustl(myline(1)))//' = '!//trim(adjustl(outline(2)))//' '//trim(adjustl(myline(2)))//' '//trim(adjustl(myline(2)))
          mystats(npar)=trim(adjustl(narray(thispar)))//' ='!//trim(adjustl(outline(2)))//' '//trim(adjustl(myline(2)))//' '//trim(adjustl(myline(2)))
          do i = 1, 5
             mystats(npar) = trim(adjustl(mystats(npar)))//' '//trim(adjustl(outline(i)))
          end do

       end if
    end do
    close(fnum)

    allocate(stats(npar))
    stats = mystats(1:npar)

  end subroutine ReadLikestats

  subroutine MakeBFname(str1,str2)
    implicit none
    character(len=*) :: str1, str2
    character(len=len(str1)) :: mystr, mystr2
    
    mystr = str1

    call split(mystr, mystr2,'.ini')
   
    str2 = trim(adjustl(mystr2))//'_bf.ini'
    
  end subroutine MakeBFname


  subroutine GetFileRoot(inifile,fileroot, fullroot)
    implicit none
    character(len=*) :: inifile, fileroot, fullroot
    character(len=mylen) :: thisline
    integer :: filenum, i, mystat, nslash
    logical :: isthere

    filenum = 123
    inquire(file=inifile,exist=isthere)
    if(.not.isthere)then
       write(*,*)'File '//trim(inifile)//' not found.'
       stop
    end if

    open(filenum, file=inifile)

    thisline = ''
    fileroot = ''
    mystat = 0
    do while (mystat == 0)
       read(filenum,'(A)',iostat=mystat)thisline
       if(mystat /=0 ) exit
       thisline = trim(adjustl(thisline))
       if(thisline(1:9)=='file_root')then
          fullroot = thisline
          call split(thisline,fileroot,'=')
          fileroot = trim(adjustl(thisline))
          exit
       end if
    end do

    close(filenum)
    
    if(trim(fileroot)=='')then
       stop 'file_root not found'
    end if

    nslash = 0
    do i = 1, len(trim(fileroot))
       if(fileroot(i:i)=='/')nslash = nslash+1
    end do
    
    do i = 1, nslash
       call split(fileroot, thisline, '/')
    end do
       
    write(*,*)'Fileroot:',trim(fileroot)
    
    

  end subroutine GetFileRoot

  subroutine GetFlag(mystring)
    implicit none
    character(len=*) :: mystring
    integer :: nargs, gi
    nargs = iargc()
    
    if(nargs==0)call DisplayHelp()
!    if(allocated(allargs))deallocate(allargs)
!    allocate(allargs(nargs))
!    do gi = 1, nargs
    call GetArg(1,mystring)
!       if(trim(allargs(gi)) == '-h')call DisplayHelp
!    end do

    if(mystring(1:7)/='params_')mystring='params_'//trim(adjustl(mystring))
    gi=len(trim(mystring))
    if(mystring(gi-3:gi)/='.ini')mystring=trim(adjustl(mystring))//'.ini'

  end subroutine GetFlag

  subroutine DisplayHelp()
    write(*,*)'What parameter file?'
    stop
  end subroutine DisplayHelp

  subroutine prepparams(arr)
    implicit none
    character(len=*) :: arr(:)
    character(len=3) :: nstring
    integer :: i

    do i = 1, size(arr)
       write(nstring,'(I3)')i
       nstring = adjustl(nstring)
       arr(i) = 'param'//trim(nstring)
    end do

  end subroutine prepparams

  
  subroutine split(string,word,seperator)
    character(len=*)::string,word,seperator  
    integer :: i,sep
    string = adjustl(string)
    sep = len(seperator)-1
    do i = 1,len(trim(string))!+2 
       if(string(i:i+sep).eq.seperator)then
          word=trim(string(1:i-1))
          string=trim(string(i+sep+1:))
          exit
        end if
     end do
   end subroutine split

end program MakeBestfitIni
