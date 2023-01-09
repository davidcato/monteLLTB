! This file is included in lltb_background.f90


	subroutine lltb_ricci(H0_inf, Lambda, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, t, theta, &
						R4, R3, R4dd, R3dd,Metricdd,RikRkj,RklRkl,R_Rij)
		! ---------------------------------------------------------- !
		! INPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(in) :: H0_inf, Lambda
		interface
			function kofr(rr)
				use lltb_precision
				real(dl) :: kofr
				real(dl), intent(in) :: rr
			end function kofr
		end interface
		interface
			function dkdr(rr)
				use lltb_precision
				real(dl) :: dkdr
				real(dl), intent(in) :: rr
			end function dkdr
		end interface
		interface
			function tbbofr(rr)
				use lltb_precision
				real(dl) :: tbbofr
				real(dl), intent(in) :: rr
			end function tbbofr
		end interface
		interface
			function dtbbdr(rr)
				use lltb_precision
				real(dl) :: dtbbdr
				real(dl), intent(in) :: rr
			end function dtbbdr
		end interface
		real(dl), intent(in) :: r, t, Mtilde
    real(dl), intent(in), optional :: theta
		! ---------------------------------------------------------- !
		! OUTPUT:
		! ---------------------------------------------------------- !
		real(dl), intent(out), optional :: R4, R3, R4dd(4,4), R3dd(3,3), Metricdd(4,4)
    real(dl), intent(out), optional :: RikRkj(3,3), RklRkl, R_Rij(3,3)
		! ---------------------------------------------------------- !
		! LOCAL:
		! ---------------------------------------------------------- !
    real(dl) :: M2, thisk, thisdkdr, Rltb, Rp, mtheta, a, H, add, ad, Rpd, Rd, Rdotdot, S, ap
    real(dl) :: mRikRkj(3), mRklRkl, mR_Rij(3)
    logical, save :: firstcall = .true.
    integer :: i
    
!    if(firstcall)then
!      write(0,*)"Ricci subroutine has not been tested. But should be fine, as it uses Mathematica's FortranForm."
!      firstcall = .false.
!    end if
    
    if(present(theta))then
      mtheta = theta
    else
      mtheta = pi * 0.5d0
    end if
    
    call lltb_functions(H0_inf, Lambda, kofr, dkdr, tbbofr, dtbbdr, Mtilde, &
						r, t, &
            Rltb=Rltb,Rp=Rp,a=a,H=H,add=add,Rpd=Rpd, S=S)
    
    thisk = kofr(r)
    thisdkdr = dkdr(r)
    M2 = Mtilde**2
    ad = a * H
    Rd = r * ad
    Rdotdot = r * add
    
    
    ! Long live Mathematica OutputForm[FortranForm[%]] !
!    if(present(R3))R3 = (-4*M2*r*(r*Rltb*thisdkdr + 2*Rltb*thisk + r*Rp*thisk))/(Rltb**2*Rp)
    if(present(R3))R3 = (-4*M2*(Rltb*thisdkdr + 2*a*thisk + Rp*thisk))/(a**2*Rp)

    if(present(R3dd))then
      R3dd = 0.d0
    
!      R3dd(1,1) = (-2*M2*r*Rp*(r*thisdkdr + 2*thisk))/(Rltb*(1 + 2*M2*r**2*thisk))
      R3dd(1,1) = (-2*M2*Rp*(r*thisdkdr + 2*thisk))/(a*(1 + 2*M2*r**2*thisk))
      R3dd(2,2) = -((M2*r*(r*Rltb*thisdkdr + 2*Rltb*thisk + 2*r*Rp*thisk))/Rp)
      R3dd(3,3) = -((M2*r*(r*Rltb*thisdkdr + 2*Rltb*thisk + 2*r*Rp*thisk)*dSin(mtheta)**2)/Rp)
    end if
    
! Need to make for R4(dd): Rdotdot, Rdot    
    if(present(R4))R4 = (2*(Rd**2*Rp + 2*Rdotdot*Rltb*Rp + 2*Rd*Rltb*Rpd + Rltb**2*Rpdd - 2*M2*r**2*Rltb*thisdkdr - 4*M2*r*Rltb*thisk - 2*M2*r**2*Rp*thisk))/(Rltb**2*Rp)
    if(present(R4dd))then
      R4dd = 0.d0
    
      R4dd(1,1) = -((2*Rdotdot*Rp + Rltb*Rpdd)/(Rltb*Rp))
      R4dd(2,2) = -((Rp*(-2*Rd*Rpd - Rltb*Rpdd + 2*M2*r**2*thisdkdr + 4*M2*r*thisk))/(Rltb*(1 + 2*M2*r**2*thisk)))
      R4dd(3,3) = -((-(Rd**2*Rp) - Rdotdot*Rltb*Rp - Rd*Rltb*Rpd + M2*r**2*Rltb*thisdkdr + 2*M2*r*Rltb*thisk + 2*M2*r**2*Rp*thisk)/Rp)
      R4dd(4,4) = -(((-(Rd**2*Rp) - Rdotdot*Rltb*Rp - Rd*Rltb*Rpd + M2*r**2*Rltb*thisdkdr + 2*M2*r*Rltb*thisk + 2*M2*r**2*Rp*thisk)*dsin(mtheta)**2)/Rp)
    end if
    
    if(present(Metricdd))then
      Metricdd = 0.d0
      Metricdd(1,1) = -1.d0
      Metricdd(2,2) = S**2
      Metricdd(3,3) = R**2
      Metricdd(4,4) = R**2 * dsin(mtheta)**2
      
    end if
    
    
    if(present(RikRkj))then
      mRikRkj =(/(4*M2**2*(r*thisdkdr + 2*thisk)**2)/(a**2*(1 + &
          2*M2*r**2*thisk)),(M2**2*r**2*(2*ap*r*thisk + a*(r*thisdkdr + &
          4*thisk))**2)/(a**2*(a + ap*r)**2),(M2**2*r**2*(2*ap*r*thisk + &
          a*(r*thisdkdr + 4*thisk))**2*Sin(mtheta)**2)/(a**2*(a + ap*r)**2)/)
      RikRkj=0.d0
      do i = 1, 3
        RikRkj(i,i)=mRikRkj(i)
      end do
    end if
    if(present(RklRkl))then
      mRklRkl = (2*M2**2*(4*ap**2*r**2*thisk**2 + 4*a*ap*r*thisk*(r*thisdkdr + &
          4*thisk) + a**2*(3*r**2*thisdkdr**2 + 16*r*thisdkdr*thisk + &
          24*thisk**2)))/(a**4*(a + ap*r)**2)
      RklRkl = mRklRkl
    end if
    if(present(R_Rij))then
      mR_Rij = (/(8*M2**2*r**2*(r*thisdkdr + 2*thisk)*(r*Rp*thisk + &
          Rltb*(r*thisdkdr + 2*thisk)))/(Rltb**3*(1 + &
          2*M2*r**2*thisk)),(4*M2**2*r**2*(r*Rltb*thisdkdr + 2*(Rltb + &
          r*Rp)*thisk)*(r*Rp*thisk + Rltb*(r*thisdkdr + &
          2*thisk)))/(Rltb**2*Rp**2),(4*M2**2*r**2*(r*Rltb*thisdkdr + 2*(Rltb + &
          r*Rp)*thisk)*(r*Rp*thisk + Rltb*(r*thisdkdr + &
          2*thisk))*Sin(mtheta)**2)/(Rltb**2*Rp**2)/)
      R_Rij=0.d0
      do i = 1, 3
        R_Rij(i,i) = mR_Rij(i)
      end do
    end if
  end subroutine lltb_ricci
		