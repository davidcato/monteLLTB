
		real(pr), parameter :: f1=3._pr/14._pr, f2=1._pr/6._pr, f3=9._pr/88._pr
		real(pr), parameter :: f4=3._pr/22._pr, f5=9._pr/52._pr, f6=3._pr/26._pr
		real(pr), parameter :: f7=1._pr/16._pr, f8=3._pr/40._pr, f9=3._pr/20._pr
		real(pr), parameter :: f10=45._pr/272._pr, f11=9._pr/68._pr
		real(pr), parameter :: f01=1._pr/3._pr
		real(pr) :: m4, q

		a0 = (x+y+3._pr*z)*0.2_pr
		am=a0
		xm=x
		ym=y
		zm=z
		q=(4._pr/precision(pr))**(1._pr/6._pr) * max(abs(a0-xm),abs(a0-ym),abs(a0-zm)) ! don't understand why, but something's wrong: set precision very high now.
		m4=1._pr
		sum1 = cmplx(0._pr,0._pr)

		do while (m4*q>=abs(am))
!write(*,*)q,am,m4*q,abs(am)
			sqx=sqrt(xm)
			sqy=sqrt(ym)
			sqz=sqrt(zm)

!			lm = sqx*sqy + sqx*sqz + sqy*sqz
			lm = sqx*(sqy + sqz) + sqy*sqz

			!m=0
			sum1 = sum1 + 3._pr * m4 / (sqz * (zm + lm) )


			!m=m+1
			am = 0.25_pr * (am + lm)
			xm = 0.25_pr * (xm + lm)
			ym = 0.25_pr * (ym + lm)
			zm = 0.25_pr * (zm + lm)
!			am = f00 * (xm+ym+zm)
			m4 = m4*0.25_pr

			! exit when condition is met:
			if(m4*q<abs(am))exit
		end do
		XX = (a0 - x) * m4 / am
		YY = (a0 - y) * m4 / am
!		XX=(am-xm)/am
!		YY=(am-ym)/am
		ZZ = (-XX-YY)*f01

		E2 = XX*YY - 6._pr*ZZ**2
		E3 = (3._pr*XX*YY - 8._pr*ZZ**2)*ZZ
		E4 = 3._pr*(XX*YY - ZZ**2)*ZZ**2
		E5 = XX*YY*ZZ**3
		!rrdd = m4/(am**(1.5_pr)) * (1._pr - f1*E2 + f2*E3 + f3*E2**2 - f4**E4 - f5*E2*E3 + f6 * E5  &
		!       - f7 * E2**3 + f8 * E3**2 + f9*E2*E4 + f10 * E2**2 * E3 - f11 * (E3*E4 + E2*E5)) + sum1


		rrdd = m4 * am**(-1.5_pr) * (1._pr - f1 * E2 + f2 * E3 + f3 * E2**2 &
				- f4 * E4 - f5 * E2 * E3 + f6 * E5 )+ sum1

		if (abs(maimag(rrdd)/abs(rrdd)) .lt. precision(pr)) rrdd = killimag(rrdd)
		if (abs(mreal(rrdd)/abs(rrdd)) .lt. precision(pr)) rrdd = killreal(rrdd)

