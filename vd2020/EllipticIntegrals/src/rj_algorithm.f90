!----

	real(pr) :: q, m4

	real(pr), parameter :: f1=3._pr/14._pr, f2=1._pr/6._pr, &
		 f3=9._pr/88._pr, f4=3._pr/22._pr, f5=9._pr/52._pr, f6=3._pr/26._pr
	logical :: specialkees

	specialkees = .false.

	specialkees =  ((real(p) <= 0.0) .and. all((/ real(x)>0.,real(y)>0.,real(z)>0.  /)).and. all((/ maimag(x)==0.,maimag(y)==0.,maimag(z)==0.  /)))

	if(.not.specialkees) then
		xm=x
		ym=y
		zm=z
		pm=p
	else

		! sort such that (z-y)(y-x)>0
		! real also takes non-complex arguments, aimag doesn't, maimag does.
		xm=x
		if(real(y)<real(xm))xm=y
		if(real(z)<real(xm))xm=z
		zm=x
		if(real(y)>real(zm))zm=y
		if(real(z)>real(zm))zm=z
		ym=x+y+z-xm-zm

		spec_q = - p
		spec_z = zm
		spec_y = ym
		spec_p = (zm-ym)*(ym-xm)/(ym+spec_q) + ym

		spec_pq = -spec_q * spec_p / ym
		spec_xypq = xm*zm/ym
		spec_fac = 3._pr

		pm = spec_p

	end if

	a0=(xm+ym+zm+pm+pm)*0.2_pr
	delta=(pm-xm)*(pm-ym)*(pm-zm)
	q=(precision(pr)*0.25_pr)**(-1._pr/6._pr) * max(abs(a0-xm),abs(a0-ym),abs(a0-zm),abs(a0-pm))

	am=a0

	m4=1._pr
	sum1 = 0._pr

	rc_arg1 = 0._pr
	rc_arg1 = rc_arg1 + 1._pr

	do while (m4 * q > abs(am))

		! m=0
!		lm = sqrt(xm)*sqrt(ym)+sqrt(xm)*sqrt(zm)+sqrt(ym)*sqrt(zm)
!		dm = (sqrt(pm)+sqrt(xm))*(sqrt(pm)+sqrt(ym))*(sqrt(pm)+sqrt(zm))
		sqx = sqrt(xm)
		sqy = sqrt(ym)
		sqz = sqrt(zm)
		sqp = sqrt(pm)
!		lm = sqx*sqy+sqx*sqz+sqy*sqz
		lm = sqx*(sqy+sqz)+sqy*sqz
		dm = (sqp+sqx)*(sqp+sqy)*(sqp+sqz)
		em = m4**3 * delta / dm**2

		rc_arg2 = 1._pr + em

		sum1 = sum1 + m4 / dm * rc(rc_arg1, rc_arg2)

		! m=m+1
		m4 = m4 * 0.25_pr

		xm = 0.25_pr * ( xm + lm )
		ym = 0.25_pr * ( ym + lm )
		zm = 0.25_pr * ( zm + lm )
		pm = 0.25_pr * ( pm + lm )

		am = 0.25_pr * ( am + lm )

!		if(m4 * q < abs(am)) exit
	end do

	XX = (a0 - x) * m4 / am
	YY = (a0 - y) * m4 / am
	ZZ = (a0 - z) * m4 / am
	PP = -0.5_pr * (XX + YY + ZZ)

	E2 = XX*YY + XX*ZZ + YY*ZZ - 3._pr * PP**2
	E3 = XX*YY*ZZ + 2._pr * E2 * PP + 4._pr * PP**3
	E4 = (2._pr * XX*YY*ZZ + E2*PP + 3._pr*PP**3) * PP
	E5 = XX*YY*ZZ*PP**2

	rrjj = 6._pr * sum1 + m4 * am**(-1.5_pr) * (1._pr - f1 * E2 + f2 * E3 + f3 * E2**2 &
			- f4 * E4 - f5 * E2 * E3 + f6 * E5 )

	if(specialkees)then
		rrjj = (( (spec_p - spec_y) * rrjj - 3._pr * rf(xm,ym,zm) + spec_fac * rc(spec_xypq,spec_pq) )) / (spec_y + spec_q)
	end if

	if (abs(maimag(rrjj)/abs(rrjj)) .lt. precision(pr)) rrjj = killimag(rrjj)
	if (abs(mreal(rrjj)/abs(rrjj)) .lt. precision(pr)) rrjj = killreal(rrjj)

!----
