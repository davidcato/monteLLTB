program TestEllipticLTB
use elliptics
use lltb_precision
use lltb_params
use lltb_prefuncs
!use lltb
use lltb_apprec
use lltb_background_ht
use lltb_background_aofht
use lltb_background_errors

implicit none

	integer :: i, imax, j, k
	real(dl) :: myt, mya, inline(6), averr, thiserr,averrlog
	character(len=128), parameter :: filename="../test.txt"
	integer, parameter :: fnum=10
	integer :: mystat, discarded, errorcounter(2*dl),thetime(3)
character(len=32) :: dummmmmy
	type(voidparams) :: vp
	logical :: prev_bad
	character(len=256) :: outputstring
	character(len=1) :: testkind
	real(dl) :: trueval, ourval, sigma
        logical :: testonemodel

!        call testap()

	testkind="a"
        testonemodel=.true.
!        testonemodel=.false.

	if(.not.testonemodel)open(fnum,file=filename,status="old")

	vp%xyz(1) = 8._dl / 3._dl * pi * vp%mtilde**2
	vp%xyz(2) = 2._dl * vp%mtilde**2 * vp%k
	vp%xyz(3) = vp%lambda / 3._dl


	lfeedback=3
!	lfeedback=1
	if(lfeedback>0)call system_clock(count=thetime(1),count_rate=thetime(3))

	imax = 100
	mystat=0
	i=0
	discarded=0
	averr = 0._dl
	averrlog = 0._dl
	prev_bad=.false.
	errorcounter = 0
	sigma=0._dl
	do while (mystat==0)
           if(.not. testonemodel)then
              read(fnum,*,iostat=mystat)inline
              if(mystat/=0)exit
           else
!x
inline(1) = 0.142857142857142876968E+00_dl
 !y
inline(2) = -0.511012821070062073225E+00_dl
 !z
inline(3) = 0.333333333333333314830E+00_dl
 !a
inline(4) = 0.000000000000000000000E+00_dl
 !t
inline(5) = 0.148201896080525003695E-04_dl
           end if
                vp%xyz(1)=inline(1)
		vp%xyz(2)=inline(2)
		vp%xyz(3)=inline(3)
		mya = inline(4)
		myt = inline(5)

                if(.not.testonemodel)then
                   if(inline(3)<0.)cycle ! Don't care about negativa lamba
                   if(inline(1)<1.d-4)cycle ! Don't care about negativa omega_m
                   !		if(  (inline(3)>1.5_dl)  .and.   ( abs(inline(2))>1.5_dl ) ) cycle ! too extreme
                   if(  (abs(inline(3))+abs(inline(2))>3._dl) ) cycle ! too extreme
                end if
		if(testkind=="t")then
			trueval = myt
			call lltbbackgroundHtofa(vp,a=mya,t=myt)
			ourval = myt
		else if (testkind == "a") then
			trueval = mya
			call lltbbackgroundaofHt(vp,a=mya,t=myt, sigma=sigma)
			ourval = mya
		else
			stop 'Specify testkind: a or t.'
		end if


                if(.not.testonemodel)then
                   if(myt<-.99_dl)then
                      if(.not.(myt<-1.99_dl))discarded = discarded + 1
                      cycle ! Because unperturbative case.
                   end if
                end if
!		if(testkind=="t")then
!			thiserr=abs((myt-inline(5))/(0.5_dl*(abs(myt)+abs(inline(5)))))
!		else if (testkind == "a") then
!			thiserr=abs((mya-inline(4))/(0.5_dl*(abs(mya)+abs(inline(4)))))
!		else
!			stop 'Specify testkind: a or t.'
!		end if
		thiserr=abs((ourval-trueval)/(0.5_dl*(abs(ourval)+abs(trueval))))
		i=i+1
! writE(*,'(I6,10ES19.10,I6)')i,inline(1:5),trueval,ourval,thiserr,averr, averrlog,discarded
!		if(i>88407)stop 'check it'
		if((thiserr>1.d-3) .or. (prev_bad .and. (thiserr>1.d-4)))then
			if(abs(inline(1))<1.d-5 .and. abs(inline(2))<1.d-3) then
				if(.not.testonemodel)cycle
			else
  !			call MyLittleDebugger(write=.true.)
                           if(.not.testonemodel)then
                              if(inline(4)>vp%smallestrealroot)cycle ! turning point!
                              if(inline(4)>vp%smallestlinroot)cycle ! turning point!
                           end if
				writE(*,'(I6,10ES19.10,I6)')i,inline(1:5),trueval,ourval,thiserr,averr, averrlog,discarded
				write(*,*)'type something to continue'
				call lltb_error('Error returned in test driver too big.')
				read(*,*)dummmmmy
			end if
!			stop 'check it'
		else if (thiserr>1.d-4) then
		 	prev_bad=.true.
		 else
		 	prev_bad=.false.
		end if
		if(.not.isnan(thiserr))then
			if(abs(thiserr/2._dl-1._dl)>1.e-30_dl)averr=averr*dble(i-1)/dble(i)+ thiserr / dble(i)
			if(thiserr/=0._dl .and. abs(thiserr/2._dl-1._dl)>1.e-30_dl)	averrlog=averrlog*dble(i-1)/dble(i)+ log10(thiserr) / dble(i)
			j = - int(log10(thiserr) - 0.5)
			if(j>2*dl)j=2*dl
			if(j<1) then
				if(abs(thiserr)<1.d-16)then
					j=2*dl
				else
					write(*,*)thiserr
					write(*,*)log10(thiserr)
					write(*,*)j
					j=1
					stop
				end if
			end if
			errorcounter( j ) = errorcounter( j ) + 1
		end if

		if(lfeedback>0)then
			if(modulo(i,10000)==0)then
				call system_clock(count=thetime(2))
				write(*,'(I9,A,1F9.5,A)')i,' points in ',dble(thetime(2)-thetime(1))/dble(thetime(3)),' seconds.'
				write(*,'(16I9)')errorcounter
				call system_clock(count=thetime(1),count_rate=thetime(3))
			end if
		end if
                if(testonemodel)exit
	end do

	if(.not.testonemodel)close(fnum)

! Report results:

	write(*,*)'Finished all.'
	write(*,*)'Last record:'
	writE(*,'(I6,9ES19.10,I6)')i,inline(1:4),inline(5),myt,thiserr,averr, averrlog,discarded
	write(*,*)
	write(*,*)'Error counts (relative error is 10^{xx}):'
	write(*,'(A)')'       -1       -2       -3       -4       -5       -6       -7       -8       -9      -10      -11      -12      -13      -14      -15      -16'
	outputstring = '     10'
	do j = 1, 2*dl
 		outputstring = trim(adjustl(outputstring))//'       10'
	end do
	write(*,'(A)')trim('      '//outputstring)
	write(*,'(16I9,A,I9)')errorcounter,' absolute counts out of total:',i
	write(*,'(16F9.1,A)')dble(errorcounter)/dble(i)*100.,'  %'

end	program TestEllipticLTB
