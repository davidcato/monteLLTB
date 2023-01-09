!----

	real(pr) :: q, m4

	real(pr), parameter :: f1=0.3_pr, f2=1._pr/7._pr, &
		 f3=3._pr/8._pr, f4=9._pr/22._pr, f5=159._pr/208._pr, f6=9._pr/8._pr
	logical :: specialkees

	specialkees = .false.

	specialkees =  ((mreal(y) < 0.0) .and. (maimag(y)==0.))

	if(.not.specialkees) then
		xm=x
		ym=y
	else

		xm = x-y
		ym = -y
	end if

	a0=(xm+2._pr*ym)/3._pr
	q=(precision(pr)*3._pr)**(-1._pr/8._pr) * abs(a0-xm)

	am=a0

	m4=1._pr

	pm = ym

	do while (m4 * q >= abs(am))

		! m=0
!		lm = sqrt(xm)*sqrt(ym)+sqrt(xm)*sqrt(zm)+sqrt(ym)*sqrt(zm)
!		dm = (sqrt(pm)+sqrt(xm))*(sqrt(pm)+sqrt(ym))*(sqrt(pm)+sqrt(zm))
		sqx = sqrt(xm)
		sqy = sqrt(ym)

		lm = 2._pr*sqx*sqy + ym

		! m=m+1
		m4 = m4 * 0.25_pr

		xm = 0.25_pr * ( xm + lm )
		ym = 0.25_pr * ( ym + lm )

		am = 0.25_pr * ( am + lm )

!		if(m4 * q < abs(am)) exit
	end do

	ZZ = (pm - a0) * m4 / am

	rrcc = (1._pr + f1 * ZZ**2 + f2 * ZZ**3 + f3 * ZZ**4 + f4 * ZZ**5 + f5 * ZZ**6 + f6 * ZZ**7)/sqrt(am)
	if(specialkees)then
		rrcc = rrcc * sqrt(x / (x-y))
	end if

	if (abs(maimag(rrcc)/abs(rrcc)) .lt. precision(pr)) rrcc = killimag(rrcc)
	if (abs(mreal(rrcc)/abs(rrcc)) .lt. precision(pr)) rrcc = killreal(rrcc)

!----
