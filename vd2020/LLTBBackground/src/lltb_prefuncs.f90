module lltb_prefuncs
use lltb_precision
use lltb_complex_tools
use lltb_params
use lltb_background_errors

contains


	subroutine void_uvw(xyz,uvw)
		implicit none
		real(dl), intent(in) :: xyz(3)
		complex(dl), intent(out) :: uvw(3)
		real(dl):: smallestposrealroot
		complex(dl) :: Q, Q13, u, v, w, Q23, Qpart1, Qpart2, Qavv
		integer :: posrealroots
		integer :: i, realroots
		real(dl) :: x, y, z
		real(dl), parameter :: f2p23=2._dl**(2._dl/3._dl), f2p13=2._dl**(1._dl/3._dl), f2p43=2._dl**(4._dl/3._dl)
		real(dl), parameter :: f3p23=3._dl**(2._dl/3._dl), f3p13=3._dl**(1._dl/3._dl), f3p32=3._dl**(1.5_dl)
		real(dl), parameter :: f23p13 = (2._dl/3._dl)**(1._dl/3._dl)
		real(dl), parameter :: sq3=sqrt(3._dl), third=1._dl/3._dl
        real(dl) :: myscale
		integer :: uvwcase

              
                x = xyz(1)
                y = xyz(2)
                z = xyz(3)
                ! This seems to improve accuracy in Q, preventing overflow:
                myscale = x+y+z! = H0^2
                x = xyz(1)/myscale
                y = xyz(2)/myscale
                z = xyz(3)/myscale

		if(y==0.d0)then	
			uvwcase = 2
		else
			Q = (dcmplx(1._dl,0._dl) * dcmplx(81._dl * z**4 + 12._dl * y**3 * z**3 / x**2 ,0._dl)**0.5_dl -dcmplx(9._dl*z**2 ,0._dl) )*x ! x = always real and positive
            if(abs(Q)/abs(9._dl*x*z**2)<sqrt(accuracy(dl))*1.d2) then
				uvwcase = 2
			else
				uvwcase = 1
			end if
		end if



		if(uvwcase==1)then
			if(abs(Q)/abs(9._dl*x*z**2)<1.d-6) then
				write(*,*)(third / dabs(z)**(2._dl*third) / dabs(x)**(third) * y)**2, sqrt(accuracy(dl)*acceptance)
				write(*,*)'xyz:',xyz
				call lltb_error('Wanting to call MyQuadUVW, but that routine should be rudimentary: does the expansion of uvw in y not suffice?')
				stop 'change expansion conditions in uvw.'
				call MyQuadUVW(xyz,u,v,w) ! because than Q is accurate at less than 3 digits.
			else
				Q13 = Q**third
				Q23 = Q**(2._dl/3._dl)

				u = Q13/z/f2p13/f3p23 - f23p13 * y / Q13
				v = dcmplx(1._dl,sq3) * y / Q13 / f2p23 / f3p13 - dcmplx(1._dl,-sq3) * Q13 / z / f2p43 / f3p23
				if(all((/mreal(v),maimag(v)/)<accuracy(dl)*1e5_dl))	v = dcmplx(1._dl,sq3) * y / (Q *12._dl)**(third) - dcmplx(1._dl,-sq3) / z * (Q /144._dl)**third
				w = dcmplx(1._dl,-sq3) * y / Q13 / f2p23 / f3p13 - dcmplx(1._dl,sq3) * Q13 / z / f2p43 / f3p23
				if(all((/mreal(w),maimag(w)/)<accuracy(dl)*1e5_dl))  w = dcmplx(1._dl,-sq3) * y / (Q *12._dl)**(third) - dcmplx(1._dl,sq3) / z * (Q /144._dl)**third
			end if
		else if (uvwcase == 2) then


			u = - dcmplx(x,0)**third / dcmplx(z,0)**third
			v = dcmplx(-1,0)**third * dcmplx(x,0)**third / dcmplx(z,0)**third
			w = - dcmplx(-1,0)**(2._dl*third) * dcmplx(x,0)**third / dcmplx(z,0)**third


			! linear expansion in small y:
			if ( y /=0.d0)then
				u = u + third / dcmplx(z,0)**(2._dl*third) / dcmplx(x,0)**(third) * y
				v = v + dcmplx(-1,0)**(2._dl*third) *third / dcmplx(z,0)**(2._dl*third) / dcmplx(x,0)**(third) * y
				w = w - dcmplx(-1,0)**(third) *third / dcmplx(z,0)**(2._dl*third) / dcmplx(x,0)**(third) * y
			end if

		else 

			stop 'Wrong uvwcase.'

		end if


		uvw(1) = u
		uvw(2) = v
		uvw(3) = w

		! put cc's in end, real soln in begin:
                ! IMPORTANT: 
                ! WV new 21/02/2011: added norm 'sum(abs(uvw))'
                ! which should take care of some residual numerical noise in
                ! complex numbers, but has not been properly tested.
		call SortComplexes(uvw,sum(abs(uvw)))


		realroots=0
		do i = 1, 3
			if (abs( abs(aimag(uvw(i)))/abs((uvw(i)))) < accuracy(dl) ) then ! is real, and imaginary part should be killed in sortcomplexes already.
				realroots = realroots+1
			end if
		end do
		if(realroots==2)then ! numerical error: the third must be real too then!
			do i = 1, 3
				uvw(i)=killimag(uvw(i))
			end do
			! redo the sorting
			call SortComplexes(uvw)
		end if




	end subroutine void_uvw



	subroutine lltb_getallroots(vpi, ThisSubCase)
		implicit none
		type(voidparams), intent(inout) :: vpi
		character(len=1), intent(in) :: ThisSubCase(3)
		real(dl) :: thisxyz(3), dum_r, this_limit
		integer :: dum_i, i
                real(dl) :: checknansfromreals(3)

		type getroots_cache
			real(dl) :: xyz(3), smallestrealroot, smallestlinroot
			integer :: posrealroots, smallestrootindex, smallestlinrootindex
			complex(dl) :: uvw(3), uvw_linear(3)
			character(len=1) :: TSC(3)
		end type getroots_cache
		type(getroots_cache), save :: grc
		integer, parameter :: localfeedback=0
		

                real(dl) :: tot_mag
			call void_makeworkxyz(vpi,1._dl)
			vpi%crossedbranchcut = .false.

			! Caching uvw:
			if(.not. all(vpi%xyz==GRC%xyz) )then
				if(localfeedback>0)write(*,*)'NOT using cached roots.'
				call void_uvw(vpi%workxyz,vpi%uvw)!, vpi%posrealroots,vpi%smallestrealroot)
				grc%xyz=vpi%xyz
				grc%uvw=vpi%uvw
			else
				if(localfeedback>0)write(*,*)'Using cached roots.'

				vpi%uvw=grc%uvw

				vpi%smallestrealroot = grc%smallestrealroot
				vpi%posrealroots=grc%posrealroots
				vpi%smallestrootindex = grc%smallestrootindex

				if(all(ThisSubCase==grc%TSC) )then
					vpi%smallestlinroot = grc%smallestlinroot
                                        vpi%smallestlinrootindex = grc%smallestlinrootindex
					vpi%uvw_linear=grc%uvw_linear
                                        if(localfeedback>0)write(*,*)'Using cached linear roots. Returning.'
! RETURN IF CACHE IS USABLE:
					return
				end if

			end if

			vpi%uvw_linear = vpi%uvw


			if(localfeedback>0)write(*,*)'ThisSubCase=',ThisSubCase

			if(all(ThisSubCase=="L"))then

				continue

			!if(TSCstr=="LSL" .or. TSCstr=="LZL")then ! B
			else if( (ThisSubCase(2)/="L") .and. (ThisSubCase(1)=="L") .and. (ThisSubCase(3)=="L") ) then

				thisxyz=vpi%xyz
				thisxyz(2)=0._dl
				call void_uvw(thisxyz,vpi%uvw_linear)
				
#ifdef FMLIB
#define DOHIGHPRECISION
#endif
#ifndef DOHIGHPRECISION
#ifdef NOQUAD
#undef DOHIGHPRECISION		
#endif
#endif

#ifndef DOHIGHPRECISION
                this_limit = 1.d-4
#else
                this_limit = 1.d-14
#endif

            if(localfeedback>0)write(*,*)'condition to copy:',dabs(vpi%xyz(2))/(dabs(vpi%xyz(1))+dabs(vpi%xyz(3)))<this_limit,dabs(vpi%xyz(2))/(dabs(vpi%xyz(1))+dabs(vpi%xyz(3))),'<',this_limit
                                if(dabs(vpi%xyz(2))/(dabs(vpi%xyz(1))+dabs(vpi%xyz(3)))<this_limit)then
                                   vpi%uvw = vpi%uvw_linear
                                   if(localfeedback>0)write(*,*)'copied uvw_linear to uvw for LSL.'
                                end if

!			elseif(TSCstr=="LLS" .or. TSCstr=="LLZ")then ! C
			elseif( (ThisSubCase(3)/="L") .and. all(ThisSubCase(1:2)=="L") ) then

				vpi%uvw_linear(1)=dcmplx(-vpi%workxyz(1)/vpi%workxyz(2),0.)
				vpi%uvw_linear(2:3)=dcmplx(0.,0.)

!			elseif(TSCstr=="SLL" .or. TSCstr=="ZLL")then ! D
			elseif( (ThisSubCase(1)/="L") .and. all(ThisSubCase(2:3)=="L") ) then

				vpi%uvw_linear(1)=sqrt(dcmplx(-vpi%workxyz(2)/vpi%workxyz(3),0.))
				vpi%uvw_linear(2)=-vpi%uvw(1)
				vpi%uvw_linear(3)=dcmplx(0.,0.)

!			elseif(TSCstr=="SLS")then ! E
			elseif( (ThisSubCase(2)=="L") .and. (ThisSubCase(1)/="L") .and. (ThisSubCase(3)/="L")  ) then

				vpi%uvw_linear = dcmplx(-1._dl,0._dl)
				continue ! non-perturbative case

			!elseif(TSCstr=="SSL")then ! F
			elseif( (ThisSubCase(3)=="L") .and. all(ThisSubCase(1:2)/="L") ) then

				vpi%uvw_linear = dcmplx(-1._dl,0._dl)
				continue ! non-perturbative case

!			elseif(TSCstr=="LSS")then ! G
			elseif( (ThisSubCase(1)=="L") .and. all(ThisSubCase(2:3)/="L") ) then

				vpi%uvw_linear = dcmplx(-1._dl,0._dl)

                                if(dabs(vpi%xyz(2))/(dabs(vpi%xyz(1))+dabs(vpi%xyz(3)))<accuracy(dl)*acceptance_pure)then
                                   vpi%uvw = vpi%uvw_linear
                                   if(localfeedback>0)write(*,*)'copied uvw_linear to uvw for LSS.'
                                end if
			else

				vpi%uvw_linear = dcmplx(-1._dl,0._dl)

			end if

                        ! Now the linear roots are set. This is the moment to check for NaNs in 
                        ! full uvw (possible for LSL and others):
                        call lltb_set_error(vpi)
                        checknansfromreals=real(vpi%uvw)
                        if(any(isnan(checknansfromreals)))then
                           checknansfromreals=real(vpi%uvw_linear)
                           if(.not.any(isnan(checknansfromreals)))then
                              vpi%uvw = vpi%uvw_linear
                              if(localfeedback>0)write(*,*)'copied uvw_linear to uvw for NaNs.'
                           else
                              call lltb_error('All roots are NaN.')
                           end if
                        end if

			if(any(ThisSubCase=="Z"))vpi%uvw=vpi%uvw_linear ! This works, because non of the LSZ cases has a turning point in the near future..

			vpi%smallestrealroot = Get_Smallest_posrealval(vpi%uvw,vpi%smallestrootindex)
			vpi%smallestlinroot = Get_Smallest_posrealval(vpi%uvw_linear,vpi%smallestlinrootindex)

                        if(vpi%smallestlinroot > vpi%smallestrealroot)then
                           vpi%smallestlinroot=vpi%smallestrealroot
                           vpi%smallestlinrootindex=vpi%smallestrootindex
                        end if

                        tot_mag = sum(abs(vpi%uvw))

			vpi%posrealroots=0
			do i = 1, 3
        if(abs((vpi%uvw(i)))/=0.d0)then
          if (abs( abs(aimag(vpi%uvw(i)))/abs((vpi%uvw(i)))) < accuracy(dl) ) then ! is real, and imaginary part should be killed in sortcomplexes already.
  !				if (abs( abs(aimag(vpi%uvw(i)))/tot_mag) < accuracy(dl) ) then ! is real, and imaginary part should be killed in sortcomplexes already.
            if(mreal(vpi%uvw(i))>0._dl)then
              vpi%posrealroots = vpi%posrealroots+1
            end if
          end if
        end if
			end do

			! caching:
			grc%smallestrealroot = vpi%smallestrealroot
			grc%posrealroots=vpi%posrealroots
                        grc%smallestrootindex = vpi%smallestrootindex

			grc%smallestlinroot = vpi%smallestlinroot
                        grc%smallestlinrootindex = vpi%smallestlinrootindex
			grc%uvw_linear=vpi%uvw_linear

                        if(any(grc%uvw/=vpi%uvw))then ! Then we overwrote uvw with linear ones, cannot cache:
                           grc%xyz=-1.d30
                        end if
!			grc%uvw=vpi%uvw ! this line is not redundant: sometimes uvw gets set to uvw_linear, see above.

      if(grc%posrealroots>0)call lltb_doublecheckposrealroot(vpi%xyz,vpi%smallestrealroot)

	end subroutine lltb_getallroots

        subroutine CopyUVWtoLinear(vpi)
		implicit none
		type(voidparams), intent(inout) :: vpi

                vpi%uvw_linear = vpi%uvw
                vpi%smallestlinroot=vpi%smallestrealroot
                vpi%smallestlinrootindex=vpi%smallestrootindex

        end subroutine CopyUVWtoLinear

	subroutine void_makeworkxyz(vpi, ain) ! rescale omega's to most accurate ratios.
		implicit none
		type(voidparams), intent(inout) :: vpi
		real(dl) ::  minpar, maxpar, absxyz(3)
		real(dl), intent(in) :: ain
		real(dl), parameter :: third = 1._dl/3._dl


		vpi%workxyz = vpi%xyz

		! disable this routine:
		vpi%workxyznonzero=vpi%workxyz
		vpi%workxyzfac = 1._dl
		vpi%workrescalefac = 1._dl
		return
		! end disable.

	end subroutine void_makeworkxyz

	function minvalnonzero(arr)
		implicit none
		real(dl), intent(in) :: arr(:)
		real(dl) :: minvalnonzero
		integer :: n
		integer :: i
		real(dl) :: mymin

		n = size(arr)
		mymin=1.d30
		do i = 1, n
			if(arr(i)<mymin .and. arr(i) /= 0) mymin = arr(i)
		end do
		minvalnonzero = mymin
	end function minvalnonzero

	function GetOmegaCase(xyz,ain)
		implicit none
		real(dl), intent(in) :: xyz(:)
		real(dl), optional, intent(in) :: ain
		complex(dl) :: luvw(size(xyz))
		character(len=size(xyz)) :: thisresult
		character(len=1) :: GetOmegaCase(size(xyz))
		real(dl) :: lxyz(size(xyz)), lepsilon(size(xyz)), debug_eps(size(xyz))
		integer :: i, orderoflxyz(size(xyz)), n, j, ii(size(xyz)), largest_i
		character(len=1) :: largeorsmall(size(xyz))
		logical :: sortdone
		integer, parameter :: my_Epsilon(3,3)=Reshape((/1,2,3,2,3,1,3,1,2/),shape(my_Epsilon))
		integer, parameter :: localfeedback = 0

		n=size(xyz)
		do i = 1, n
			orderoflxyz(i)=i
		end do

		! only compare absolute values
		lxyz=abs(xyz)
		! only compare relative sizes:
		lxyz=abs(xyz)/(sum(abs(xyz)))


                GetOmegaCase="L"
                   
		do i = 1, n
                   if(lxyz(i)<tiny)then
                      GetOmegaCase(i)="Z"
                      lepsilon(i)=0._dl
                   end if
		end do
		if(localfeedback>0)write(*,*)'First step GetOmegaCase:',GetOmegaCase

		! In case we only want L's and Z's since we don't know a:
		if(.not.present(ain))then
                   ! necessary to roughly estimate this
                   ! because otherwise in void_uvw
                   ! we might encounter Q=0 due to overflows of the floating points.
                   ! if one of the parameters is tiny wrt the others,
                   ! Q might be nonzero and O(10^-7) only in the difference between to 
                   ! very large numbers...
                   lxyz=abs(xyz)/(xyz(1)+xyz(2)+xyz(3))
                   do i = 1, n
                      if(lxyz(i)**2<accuracy(dl) * rel_ACCEPTANCE(i))GetOmegaCase(i)='S'
                   end do
                   
					if(localfeedback>0)write(*,*)'Shortcut GetOmegaCase:',GetOmegaCase
                   return
		end if



		! otherwise:

                ! Ooops.. this is only valid for i=3
!		do i = 0, n-1
		do i = n-1, n-1
                   if(GetOmegaCase(i+1)=="Z")cycle
                   ! even permutions: 123, 231, 312
                   !			do j = 0, n-1
                   !				ii(j+1)=modulo(i+j,3)+1
                   !			end do
                   if ( ( xyz(my_Epsilon(i+1,2)) * ain**(xyz_apower(my_Epsilon(i+1,2))) + &
                        xyz(my_Epsilon(i+1,3)) * ain**(xyz_apower(my_Epsilon(i+1,3)))  ) /=0.d0)then
                     lepsilon(i+1) = xyz( my_Epsilon(i+1,1) ) * ain**(xyz_apower(my_Epsilon(i+1,1))) / &
                          ( xyz(my_Epsilon(i+1,2)) * ain**(xyz_apower(my_Epsilon(i+1,2))) + &
                          xyz(my_Epsilon(i+1,3)) * ain**(xyz_apower(my_Epsilon(i+1,3)))  )
                  else
                    lepsilon(i+1) =1.d30
                  end if
		end do
                lepsilon(1)=1e30 ! because omega_m should always dominate the big bang.
                lepsilon(2)=dabs(xyz(2)/xyz(1)*ain)
                
		lxyz = abs(lepsilon)

		! This should be faster:

		largest_i = 1
		do i = 2, n
			if (lxyz(i) > lxyz(largest_i))largest_i=i
		end do
		largeorsmall="L"
		do i = 1, n
			if(i==largest_i)cycle
			if(  lxyz(i)**2 < accuracy(dl) * rel_ACCEPTANCE(i) )then
				if(GetOmegaCase(i)/="Z")largeorsmall(i)="S"
				if( lxyz(i) == 0._dl  ) largeorsmall(i)="Z"
			end if
		end do

		GetOmegaCase = largeorsmall
	
		if(localfeedback>0)write(*,*)'Final GetOmegaCase:',GetOmegaCase


	contains

		subroutine mswapp(x1,x2)
			implicit none
			real(dl) :: x1, x2, x3
			x3=x1
			x1=x2
			x2=x3
		end subroutine mswapp

		subroutine mswapp_i(x1,x2)
			implicit none
			integer :: x1, x2, x3
			x3=x1
			x1=x2
			x2=x3
		end subroutine mswapp_i

	end function GetOmegaCase

#ifndef FMLIB
        subroutine MyQuadUVW(xyz,u,v,w) 
          use lltb_apprec
          implicit none
          real(dl), intent (in) :: xyz(3)
          complex(dl), intent(out) :: u,v,w
          !
          integer, parameter :: mykind=ql
          !		real(dl), intent(in) :: xyz(3)
          !		complex(dl), intent(out) :: uvw(3)
          real(dl):: smallestposrealroot
          complex(mykind) ::  Q13, Q23 , Q, thisu, thisv, thisw
          ! I've given up: just use quad precision for Q.
          real(mykind) :: x, y, z, Q_aprec(2)
          real(mykind), parameter :: f2p23=2._mykind**(2._mykind/3._mykind), f2p13=2._mykind**(1._mykind/3._mykind), f2p43=2._mykind**(4._mykind/3._mykind)
          real(mykind), parameter :: f3p23=3._mykind**(2._mykind/3._mykind), f3p13=3._mykind**(1._mykind/3._mykind), f3p32=3._mykind**(1.5_mykind)
          real(mykind), parameter :: f23p13 = (2._mykind/3._mykind)**(1._mykind/3._mykind)
          real(mykind), parameter :: sq3=3._mykind**0.5_mykind, third=1._mykind/3._mykind
          real(mykind) :: myscale

          myscale=1._mykind!
!          myscale=sum(xyz)
! This myscale business seems to fuck it up sometimes.
          x=xyz(1)/myscale
          y=xyz(2)/myscale
          z=xyz(3)/myscale


          Q = ((dcmplx(1._mykind,0._mykind) * ( 81._mykind * z**4 + 12._mykind * y**3 * z**3 / x**2 ))**0.5_mykind - dcmplx(1._mykind,0._mykind) *9._mykind*z**2 )*x ! x = always real and positive
          ! I have the impression that my own arbitrary precision implementation of the above seems to help.
			! but it is fucking slow...
#ifdef NOQUAD
          call APQ(Q_aprec,(/x,y,z/))
          Q=(dcmplx(1.d0,0.d0)*Q_aprec(1))**0.5_mykind-dcmplx(1.d0,0.d0)*Q_aprec(2)
#endif

          Q13 = Q**third
          Q23 = Q**(2._dl/3._dl)
          
          thisu = Q13/z/f2p13/f3p23 - f23p13 * y / Q13
          thisv = dcmplx(1._dl,sq3) * y / Q13 / f2p23 / f3p13 - dcmplx(1._dl,-sq3) * Q13 / z / f2p43 / f3p23
          if(all((/real(thisv),aimag(thisv)/)<accuracy(dl)*1e5_dl))	thisv = dcmplx(1._dl,sq3) * y / (Q *12._dl)**(third) - dcmplx(1._dl,-sq3) / z * (Q /144._dl)**third
          thisw = dcmplx(1._dl,-sq3) * y / Q13 / f2p23 / f3p13 - dcmplx(1._dl,sq3) * Q13 / z / f2p43 / f3p23
          if(all((/real(thisw),aimag(thisw)/)<accuracy(dl)*1e5_dl))  thisw = dcmplx(1._dl,-sq3) * y / (Q *12._dl)**(third) - dcmplx(1._dl,sq3) / z * (Q /144._dl)**third
          ! Original equation:
          ! Q = dcmplx(81._dl * x**2 * z**4 + 12._dl * y**3 * z**3,0._dl)**0.5_dl -dcmplx(9._dl*x*z**2 ,0._dl)

          if(lfeedback>1)write(*,*)'Using quad precision Q:',Q*myscale**3
			
			u=thisu
			v=thisv
			w=thisw
			! Check imaginary parts:
			call MakeRealIfReal(u,accuracy(dl) *acceptance)
			call MakeRealIfReal(v,accuracy(dl) *acceptance)
			call MakeRealIfReal(w,accuracy(dl) *acceptance)

          return
        end subroutine MyQuadUVW
#else
        subroutine MyQuadUVW(xyz,u,v,w) 
			USE FMZM

			implicit none
			real(dl), intent (in) :: xyz(3)
			complex(dl), intent(out) :: u,v,w
          !
			logical, save :: firstcall = .true.
		  ! 
			TYPE (ZM), SAVE :: xmp, ymp, zmp, ump, vmp, wmp
			TYPE (ZM), SAVE :: Q, Q13, Q23
!          real(mykind), parameter :: f2p23=2._mykind**(2._mykind/3._mykind), f2p13=2._mykind**(1._mykind/3._mykind), f2p43=2._mykind**(4._mykind/3._mykind)
!          real(mykind), parameter :: f3p23=3._mykind**(2._mykind/3._mykind), f3p13=3._mykind**(1._mykind/3._mykind), f3p32=3._mykind**(1.5_mykind)
!          real(mykind), parameter :: f23p13 = (2._mykind/3._mykind)**(1._mykind/3._mykind)
!          real(mykind), parameter :: sq3=3._mykind**0.5_mykind, third=1._mykind/3._mykind
          TYPE (ZM), SAVE :: f2p23, f2p13, f2p43
          TYPE (ZM), SAVE :: f3p23, f3p13, f3p32
          TYPE (ZM), SAVE :: f23p13
          TYPE (ZM), SAVE :: sq3, oneplusisq3, oneminisq3, third, twothird, twelve, twelvep2, three, one, two, nine, ninep2, half
			complex(dl) :: thisu, thisv, thisw

			
			if(firstcall)then
				firstcall = .false.
				! initialization..
				call FM_Set(128)
				twelve = TO_ZM(12)
				twelvep2 = TO_ZM(144)
				three = TO_ZM(3)
				one = TO_ZM(1)
				two = TO_ZM(2)
				half = one/two
				nine = TO_ZM(9)
				ninep2 = TO_ZM(81)
				third = one/three
				twothird = two/three
				sq3 = three**(one/two)
				oneplusisq3 = sq3 * TO_ZM('0 + i') + TO_ZM(1)
				oneminisq3 = sq3 * TO_ZM('0 - i') + TO_ZM(1)
				f2p23 = two**twothird
				f2p13 = two**third
				f2p43 = two**(2*twothird)
				f3p23 = three**twothird
				f3p13 = three**third
				f3p32 = three**(three/two)
				f23p13 = twothird**third
			end if

			xmp = TO_ZM(xyz(1))
			ymp = TO_ZM(xyz(2))
			zmp = TO_ZM(xyz(3))
			
			!Q = ((dcmplx(1._mykind,0._mykind) * ( 81._mykind * z**4 + 12._mykind * y**3 * z**3 / x**2 ))**0.5_mykind - dcmplx(1._mykind,0._mykind) *9._mykind*z**2 )*x ! x = always real and positive
			Q = (ninep2 * zmp**4 * xmp**2 + twelve * ymp**3 * zmp**3)**half - nine * zmp**2 * xmp
			Q13 = Q**third
			Q23 = Q**twothird

			ump = Q13/zmp/f2p13/f3p23 - f23p13 * ymp / Q13
			vmp = oneplusisq3 * ymp / Q13 / f2p23 / f3p13 - oneminisq3 * Q13 / zmp / f2p43 / f3p23
			wmp = oneminisq3 * ymp / Q13 / f2p23 / f3p13 - oneplusisq3 * Q13 / zmp / f2p43 / f3p23

			thisu = TO_DPZ(ump)
			thisv = TO_DPZ(vmp)
			thisw = TO_DPZ(wmp)

			if(all((/real(thisv),aimag(thisv)/)<accuracy(dl)*1e5_dl))	vmp = oneplusisq3 * ymp / (Q *twelve)**(third) - oneminisq3 / zmp * (Q /twelvep2)**third
			if(all((/real(thisw),aimag(thisw)/)<accuracy(dl)*1e5_dl))  wmp = oneminisq3 * ymp / (Q *twelve)**(third) - oneplusisq3 / zmp * (Q /twelvep2)**third

			thisv = TO_DPZ(vmp)
			thisw = TO_DPZ(wmp)

			u=thisu
			v=thisv
			w=thisw

			if(lfeedback>1)write(*,*)'Using multiple precision Q:',TO_DPZ(Q)
			if(lfeedback>1)then
				if(dabs(mreal(to_dpz(q)))<1.d-100)then
					write(*,*)xyz
					stop
				end if
			end	if

        end subroutine MyQuadUVW

#endif

  ! new subroutine 2012, double check that this real root is actually correct.
  ! some rare cases overestimate this real root.
  subroutine lltb_doublecheckposrealroot(xyz,realroot)
    use lltb_background_eqs
    implicit none
    real(dl), intent(in) :: xyz(3)
    real(dl), intent(inout) :: realroot

    real(dl) :: thisH2, toproot, botroot, newroot
    integer :: i

    thisH2 = HoverH0_squared(xyz,realroot,nocheck=.true.)
    newroot=realroot
    if(thisH2<0.)then
      !stop 'ladies and gentlemen, we got him.'
      i=0
      toproot=realroot
      botroot=toproot
      do while (thisH2<0.)
        i=i+1
        botroot=botroot*0.95
        thisH2 = HoverH0_squared(xyz,botroot,nocheck=.true.)
          !write(*,*)i,thisH2,botroot
        if(i>100)then
          ! write(*,*)'cannot even find botroot. thisH2, botroot:',thisH2,botroot
          realroot=1.d-3*realroot
          return
        end if
      end do
      i=0
      do while ((abs(toproot/botroot-1.d0)>1.d-12 .or. thisH2<0) .and. (thisH2 /= 0.d0))
        i=i+1
        newroot=0.5d0*(toproot+botroot)
        thisH2 = HoverH0_squared(xyz,newroot,nocheck=.true.)
        if(thisH2>0)then
          botroot=newroot
        else
          toproot=newroot
        end if
          !write(*,*)i,thisH2,botroot,toproot
        if(i>500)then
          ! write(*,*)'cannot converge. thisH2, botroot:',thisH2,botroot
          realroot=botroot !1.d-3*realroot
          return          
        end if
      end do
    end if
    
    !write(*,*)'oldroot, newroot:',realroot,newroot
    realroot=newroot 
    
  end subroutine lltb_doublecheckposrealroot

end module lltb_prefuncs
