program TestEllipticLTB
use elliptics
implicit none
integer, parameter :: dl = kind(1.d0)
real(dl) :: res1, inline(6), reldiff(2),avvdiff(2)
complex(dl) :: res2
logical :: noteof
integer :: myiostat, counter
real(dl) :: x
	real(dl), parameter :: pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068_dl


!
!res1=ellf(.9d0,1.d0)
!write(*,'(1ES50.42)') res1
!res2=ellf(dcmplx(-1.d0,1.0d0),dcmplx(0.6d0,0.8d0))
!write(*,"('(',ES40.32,',',ES40.32,')')")res2
!res2=ellf(dcmplx(-1.d0,1.0d0 + 4.d-3),dcmplx(0.6d0,0.8d0))
!write(*,"('(',ES40.32,',',ES40.32,')')")res2
!res2=ellf(dcmplx(-1.d0,1.0d0),dcmplx(0.6d0,0.8d0))
!write(*,"('(',ES40.32,',',ES40.32,')')")res2

!write(*,*)

write(*,*)1,2,0,rf(1._4,2._4,0._4)
write(*,*)1,2,0,rf(1._8,2._8,0._8)
write(*,*)-1,2,0,rf(dcmplx(-1.,0.),dcmplx(2,0),dcmplx(0,0))
!write(*,*)1,-2,0,rf(dcmplx(1.,0.),dcmplx(-2,0),dcmplx(0,0))
!write(*,*)1,-2,0,'norm',abs(rf(dcmplx(1.,0.),dcmplx(-2,0),dcmplx(0,0)))
!write(*,*)-1,-2,0,rf(dcmplx(-1.,0.),dcmplx(-2,0),dcmplx(0,0))
res2=rf(dcmplx(0,1),dcmplx(0,-1),dcmplx(0,0))
write(*,*)'i,-i,0',rf(cmplx(0,1),cmplx(0,-1),cmplx(0,0))
write(*,*)'i,-i,0',res2
!res2=rf_c2(dcmplx(0,1),dcmplx(0,-1),dcmplx(0,0))
!write(*,*)'i,-i,0',res2

write(*,*)'i-1,i,0',rf(dcmplx(-1,1),dcmplx(0,1),dcmplx(0,0))
write(*,*)'i-1,i,0',rf(cmplx(-1,1),cmplx(0,1),cmplx(0,0))

write(*,*)2,3,4,rf(2._4,3._4,4._4)
write(*,*)2,3,4,rf(2._8,3._8,4._8)

write(*,*)'i,-i,2',rf(cmplx(0,1),cmplx(0,-1),cmplx(2,0))
write(*,*)'i,-i,2',rf(dcmplx(0,1),dcmplx(0,-1),dcmplx(2,0))
!write(*,*)'i,-i,2',rf_c2(dcmplx(0,1),dcmplx(0,-1),dcmplx(2,0))

write(*,*)'i-1,i,1-i',rf(cmplx(-1,1),cmplx(0,1),cmplx(1,-1))
write(*,*)'i-1,i,1-i',rf(dcmplx(-1,1),dcmplx(0,1),dcmplx(1,-1))
 write(*,*)'----------------------'

write(*,*)0,0.25,rc(0._4,0.25_4),rc(0._4,0.25_4)-pi
write(*,*)0,0.25,rc(0._8,0.25_8),rc(0._8,0.25_8)-pi
write(*,*)0,0.25,rc(cmplx(0.,0.),cmplx(0.25,0.)),rc(cmplx(0.,0.),cmplx(0.25,0.))-pi
write(*,*)0,0.25,rc(dcmplx(0.,0.),dcmplx(0.25,0.)),rc(dcmplx(0.,0.),dcmplx(0.25,0.))-pi

write(*,*)9./4.,2.,rc(2.25d0,2._8),rc(2.25d0,2._8)-log(2.)

write(*,*)'0, i',rc(cmplx(0.,0.),cmplx(0.,1.))
write(*,*)'0, i',rc(Dcmplx(0.,0.),Dcmplx(0.,1.))

write(*,*)'-i, i',rc(cmplx(0.,-1.),cmplx(0.,1.))
write(*,*)'-i, i',rc(Dcmplx(0.,-1.),Dcmplx(0.,1.))

write(*,*)'0.25, -2',rc(0.25d0,-2.d0)

write(*,*)'i, -1',rc(Dcmplx(0.,1.),Dcmplx(-1.,0.))
 write(*,*)'----------------------'

write(*,*)'0,1,2,3:',rj(0.d0,1.d0,2.d0,3.d0)

write(*,*)'2,3,4,5:',rj(2.e0,3.e0,4.e0,5.e0)
write(*,*)'2,3,4,5:',rj(2.d0,3.d0,4.d0,5.d0)
write(*,*)'2,3,4,5:',rj(cmplx(2.d0,0),cmplx(3.d0,0),cmplx(4.d0,0),cmplx(5.d0,0))
write(*,*)'2,3,4,5:',rj(dcmplx(2.d0,0),dcmplx(3.d0,0),dcmplx(4.d0,0),dcmplx(5.d0,0))


write(*,*)'2,3,4,-1+i:',rj(cmplx(2.,0),cmplx(3,0),cmplx(4,0),cmplx(-1,1))

res2=rj(dcmplx(2.,0),dcmplx(3,0),dcmplx(4,0),dcmplx(-1,1))
write(*,*)'2,3,4,-1+i:',res2
!res2=rj_cnr(dcmplx(2.,0),dcmplx(3,0),dcmplx(4,0),dcmplx(-1,1))
!write(*,*)'2,3,4,-1+i:',res2
res2=rj(dcmplx(2.,0),dcmplx(3,0),dcmplx(4,0),dcmplx(-0.5,0))
write(*,*)'2,3,4,-0.5:',res2
res2=rj(dcmplx(2.,0),dcmplx(3,0),dcmplx(4,0),dcmplx(-5,0))
write(*,*)'2,3,4,-5:',res2

write(*,*)'i,-i,0,2:',rj(cmplx(0.,1),cmplx(0,-1),cmplx(0,0),cmplx(2,0))
write(*,*)'i,-i,0,2:',rj(dcmplx(0.,1),dcmplx(0,-1),dcmplx(0,0),dcmplx(2,0))

write(*,*)'-1+i,-1-i,1,2:',rj(cmplx(-1.,1),cmplx(-1,-1),cmplx(1,0),cmplx(2,0))
write(*,*)'-1+i,-1-i,1,2:',rj(dcmplx(-1.,1),dcmplx(-1,-1),dcmplx(1,0),dcmplx(2,0))

write(*,*)'i,-i,0,1-i:',rj(cmplx(0,1),cmplx(0,-1),cmplx(0,0),cmplx(1,-1))
write(*,*)'i,-i,0,1-i:',rj(dcmplx(0,1),dcmplx(0,-1),dcmplx(0,0),dcmplx(1,-1))

res2=rj(cmplx(-1.,1),cmplx(-1,-1),cmplx(1,0),cmplx(-3,1))
write(*,*)'-1+i,-1-i,1,-3+i:',res2
res2=rj(dcmplx(-1.,1),dcmplx(-1,-1),dcmplx(1,0),dcmplx(-3,1))
write(*,*)'-1+i,-1-i,1,-3+i:',res2

write(*,*)'-1+i,-2-i,-i,-1+i:',rj(cmplx(-1,1),cmplx(-2,-1),cmplx(0,-1),cmplx(-1,+1))
write(*,*)'-1+i,-2-i,-i,-1+i:',rj(dcmplx(-1,1),dcmplx(-2,-1),dcmplx(0,-1),dcmplx(-1,+1))


write(*,*)'----------------------'
write(*,*)'0,2,1',rd(cmplx(0,0),cmplx(2,0),cmplx(1,0))
write(*,*)'0,2,1',rd(dcmplx(0,0),dcmplx(2,0),dcmplx(1,0))
write(*,*)'2,3,4',rd(cmplx(2,0),cmplx(3,0),cmplx(4,0))
write(*,*)'2,3,4',rd(dcmplx(2,0),dcmplx(3,0),dcmplx(4,0))
write(*,*)'i,-i,2',rd(dcmplx(0,1),dcmplx(0,-1),dcmplx(2,0))
write(*,*)'0,i,-i',rd(dcmplx(0,0),dcmplx(0,1),dcmplx(0,-1))
write(*,*)'0,i-1,i',rd(dcmplx(0,0),dcmplx(-1,1),dcmplx(0,1))
write(*,*)'-2-i,-i,-1+i',rd(dcmplx(-2,-1),dcmplx(0,-1),dcmplx(-1,1))


write(*,*)'done..'

stop 'hier'

write(*,*)'Starting test routine'
myiostat=0

open(21,file="/Users/Valkenburg/Physics/Coding/LuminosityDistance/testRf.tsv")
!write(*,*)isnan(res1),res1
counter = 0
avvdiff=0.d0
do while (myiostat==0)
     read(21,*,iostat=myiostat)inline
	 if(myiostat/=0)exit
	counter = counter + 1
	 res2=rf(dcmplx(inline(1),inline(2)),dcmplx(inline(3),inline(4)),dcmplx(inline(1),inline(2)))
	 reldiff=(/ inline(5)- Real(res2) ,inline(6) - aImag(res2)  /)
	 write(*,'(4ES40.32)')reldiff
	if(any(isnan(reldiff)))then
		if(real(res2)/=inline(5) .and. aimag(res2)/=inline(6))then
			write(*,'(6ES40.32,1L5)')inline,any(isnan(reldiff)),res2
			write(*,*)real(res2)==inline(5),real(res2),inline(5)
			stop
		else
			reldiff=0.d0
		end if
	end if
	 avvdiff = avvdiff + abs(reldiff)
	write(*,*)avvdiff
	if(any(avvdiff>1.))then
		write(*,*)inline
		write(*,*)res2
		stop

	end if
end	do
write(*,*)'Average error',avvdiff/dble(counter)

closE(21)


end	program TestEllipticLTB
