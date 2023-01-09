
		real(pr), parameter :: f1=0.1_pr, f2=1._pr/14._pr, f3=1._pr/24._pr, f4=3._pr/44._pr, f00=1._pr/3._pr, f01=1._pr/6._pr
		real(pr), parameter :: f5=5._pr/208._pr, f6=3._pr/104._pr, f7=1._pr/16._pr
		real(pr) :: m4, q

		a0 = (x+y+z)*f00
		am=a0
		xm=x
		ym=y
		zm=z
!		q=(1._pr/precision(pr)/3._pr)**(1._pr/2._pr) * max(abs(a0-xm),abs(a0-ym),abs(a0-zm)) ! don't understand why, but something's wrong: set precision very high now.
		q=(1._pr/precision(pr)/3._pr)**(1._pr/8._pr) * max(abs(a0-xm),abs(a0-ym),abs(a0-zm))
		m4=1._pr

		do while (m4*q>=abs(am))
!write(*,*)q,am,m4*q,abs(am)
			sqx=sqrt(xm)
			sqy=sqrt(ym)
			sqz=sqrt(zm)
!write(*,*)sqx,sqy,sqz,'sqx,sq,sqz'

! this step may have rounding errors: ...?
!			lm = sqx*sqy + sqx*sqz + sqy*sqz
			lm = sqx*(sqy + sqz) + sqy*sqz


			am = 0.25_pr * (am + lm)
			xm = 0.25_pr * (xm + lm)
			ym = 0.25_pr * (ym + lm)
			zm = 0.25_pr * (zm + lm)
!			am = f00 * (xm+ym+zm)
			m4 = m4*0.25_pr
		end do


		XX = (a0 - x) * m4 / am
		YY = (a0 - y) * m4 / am
		XX=(am-xm)/am
		YY=(am-ym)/am
		ZZ = -XX-YY

		E2 = XX*YY - ZZ**2
		E3 = XX*YY*ZZ

		rrff = 1._pr/sqrt(am) * (1._pr - f1*E2 + f2*E3 + f3*E2**2 - f4*E2*E3 - f5 * E2**3 + f6*E3**2 + f7*E2**2 * E3)

		if (abs(maimag(rrff)/abs(rrff)) .lt. precision(pr)) rrff = killimag(rrff)
		if (abs(mreal(rrff)/abs(rrff)) .lt. precision(pr)) rrff = killreal(rrff)

