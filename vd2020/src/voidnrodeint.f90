

MODULE ode_path
        USE nrtype
        INTEGER(I4B) :: nok,nbad,kount, odeint_feedback = 0
        LOGICAL(LGT), SAVE :: save_steps=.false.
        REAL(DP) :: dxsav
        REAL(DP), DIMENSION(:), POINTER :: xp
        REAL(DP), DIMENSION(:,:), POINTER :: yp
        REAL(DP), DIMENSION(:,:), POINTER :: yperr
        logical :: nrsuccess
        logical, pointer :: OdeVoidError
END MODULE ode_path
!!$
module nrutil
  USE nrtype
  IMPLICIT NONE
  INTERFACE reallocate
     MODULE PROCEDURE reallocate_rv,reallocate_rm,&
          reallocate_iv,reallocate_im,reallocate_hv
  END INTERFACE
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE
  
interface isnan
    module procedure isnan_sc, isnan_vec
end interface isnan

contains


    function isnan_sc(var)
        implicit none
        real(kind(1.d0)) var
        logical isnan_sc

        isnan_sc = var /= var
    end function isnan_sc

    function isnan_vec(var)
        implicit none
        real(kind(1.d0)) var(:)
        logical isnan_vec(size(var))

        isnan_vec = var /= var
    end function isnan_vec

!BL
        FUNCTION reallocate_rv(p,n)
        REAL(DP), DIMENSION(:), POINTER :: p, reallocate_rv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_rv(n),stat=ierr)
        reallocate_rv=0_dp
        if (ierr /= 0) call &
                nrerror('reallocate_rv: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p)
        reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
        END FUNCTION reallocate_rv
!BL
        FUNCTION reallocate_iv(p,n)
        INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_iv(n),stat=ierr)
        reallocate_iv=0
        if (ierr /= 0) call &
                nrerror('reallocate_iv: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p)
        reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
        END FUNCTION reallocate_iv
!BL
        FUNCTION reallocate_hv(p,n)
        CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_hv(n),stat=ierr)
        reallocate_hv=""       
        if (ierr /= 0) call &
                nrerror('reallocate_hv: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p)
        reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
        END FUNCTION reallocate_hv
!BL
        FUNCTION reallocate_rm(p,n,m)
        REAL(DP), DIMENSION(:,:), POINTER :: p, reallocate_rm
        INTEGER(I4B), INTENT(IN) :: n,m
        INTEGER(I4B) :: nold,mold,ierr
        allocate(reallocate_rm(n,m),stat=ierr)
        reallocate_rm=0_dp
        if (ierr /= 0) call &
                nrerror('reallocate_rm: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p,1)
        mold=size(p,2)
        reallocate_rm(1:min(nold,n),1:min(mold,m))=&
                p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
        END FUNCTION reallocate_rm
!BL
        FUNCTION reallocate_im(p,n,m)
        INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
        INTEGER(I4B), INTENT(IN) :: n,m
        INTEGER(I4B) :: nold,mold,ierr
        allocate(reallocate_im(n,m),stat=ierr)
        reallocate_im=0
        if (ierr /= 0) call &
                nrerror('reallocate_im: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p,1)
        mold=size(p,2)
        reallocate_im(1:min(nold,n),1:min(mold,m))=&
                p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
        END FUNCTION reallocate_im
!BL
!BL
        SUBROUTINE nrerror(string)
          use ode_path
          CHARACTER(LEN=*), INTENT(IN) :: string
          if(odeint_feedback.gt.0)write (*,*) 'nrerror: ',string
          !STOP 'program terminated by nrerror'
          nrsuccess=.false.
        END SUBROUTINE nrerror

!BL
        SUBROUTINE assert1(n1,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1
        if (.not. n1) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert1'
        end if
        END SUBROUTINE assert1
!BL
        SUBROUTINE assert2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2
        if (.not. (n1 .and. n2)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert2'
        end if
        END SUBROUTINE assert2
!BL
        SUBROUTINE assert3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3
        if (.not. (n1 .and. n2 .and. n3)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert3'
        end if
        END SUBROUTINE assert3
!BL
        SUBROUTINE assert4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3,n4
        if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert4'
        end if
        END SUBROUTINE assert4
!BL
        SUBROUTINE assert_v(n,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, DIMENSION(:), INTENT(IN) :: n
        if (.not. all(n)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert_v'
        end if
        END SUBROUTINE assert_v
!BL
        FUNCTION assert_eq2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2
        INTEGER :: assert_eq2
        if (n1 == n2) then
                assert_eq2=n1
        else
                write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eq2'
        end if
        END FUNCTION assert_eq2
!BL
        FUNCTION assert_eq3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2,n3
        INTEGER :: assert_eq3
        if (n1 == n2 .and. n2 == n3) then
                assert_eq3=n1
        else
                write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eq3'
        end if
        END FUNCTION assert_eq3
!BL
        FUNCTION assert_eq4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2,n3,n4
        INTEGER :: assert_eq4
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
                assert_eq4=n1
        else
                write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eq4'
        end if
        END FUNCTION assert_eq4
!BL
        FUNCTION assert_eqn(nn,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, DIMENSION(:), INTENT(IN) :: nn
        INTEGER :: assert_eqn
        if (all(nn(2:) == nn(1))) then
                assert_eqn=nn(1)
        else
                write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                        string
                STOP 'program terminated by assert_eqn'
        end if
        END FUNCTION assert_eqn
!BL

      end module nrutil

Module NROdeInt
private
public odeint, rkqs, isnan

interface isnan
    module procedure isnan_sc, isnan_vec
end interface isnan

contains


    function isnan_sc(var)
        implicit none
        real(kind(1.d0)) var
        logical isnan_sc

        isnan_sc = var /= var
    end function isnan_sc

    function isnan_vec(var)
        implicit none
        real(kind(1.d0)) var(:)
        logical isnan_vec(size(var))

        isnan_vec = var /= var
    end function isnan_vec


!        SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
        SUBROUTINE odeint(ystart,intvarstart, GetExactList, stopvar,stopval,eps,h1,hmin,derivs,rkqs)
        USE nrtype; USE nrutil, ONLY : nrerror,reallocate
        USE ode_path
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: ystart
! WV NEW
!        REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
        REAL(DP), INTENT(IN) :: eps,h1,hmin
        real(dp) :: x1,x2, stopvarstart
        real(dp), Intent(in) :: GetExactList(:) 
        real(dp), intent(in) :: stopval, intvarstart
        integer, intent(in) :: stopvar
! END WV NEW
        INTERFACE
                SUBROUTINE derivs(x,y,dydx)
                USE nrtype
                IMPLICIT NONE
                REAL(DP), INTENT(IN) :: x
                REAL(DP), DIMENSION(:), INTENT(IN) :: y
                REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
                END SUBROUTINE derivs
!BL
                SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal, yerr,hdid,hnext,derivs)
                USE nrtype
                IMPLICIT NONE
                REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
                REAL(DP), DIMENSION(:), INTENT(OUT) :: yerr
                REAL(DP), DIMENSION(:), INTENT(IN) :: dydx,yscal
                REAL(DP), INTENT(INOUT) :: x
                REAL(DP), INTENT(IN) :: htry,eps
                REAL(DP), INTENT(OUT) :: hdid,hnext
                INTERFACE
                        SUBROUTINE derivs(x,y,dydx)
                        USE nrtype
                        IMPLICIT NONE
                        REAL(DP), INTENT(IN) :: x
                        REAL(DP), DIMENSION(:), INTENT(IN) :: y
                        REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
                        END SUBROUTINE derivs
                END INTERFACE
                END SUBROUTINE rkqs
        END INTERFACE
        REAL(DP), PARAMETER :: TINY=1.0e-30_DP
        INTEGER(I4B), PARAMETER :: MAXSTP=1000!0!500!10000
        INTEGER(I4B) :: nstp, mystat, i
        REAL(DP) :: h,hdid,hnext,x,xsav
        REAL(DP), DIMENSION(size(ystart)) :: dydx,y,yscal, yerr
! WV NEW
!        x=x1
        x1=intvarstart
        x = x1
        if (stopvar==0)then
           x2 = stopval
           stopvarstart = intvarstart
        else
           x2=1.d30
           stopvarstart = ystart(stopvar)
        end if

! END WV NEW
        h=sign(h1,x2-x1)
        nok=0
        nbad=0
        kount=0
        y(:)=ystart(:)
        if (save_steps) then
                xsav=x-2.0_DP*dxsav
                deallocate(xp,stat=mystat) ! prevent memory leak WV
                deallocate(yp,stat=mystat)
                deallocate(yperr,stat=mystat)
                nullify(xp,yp,yperr)
                allocate(xp(256))
                allocate(yp(size(ystart),size(xp)))
                allocate(yperr(size(ystart),size(xp)))
                xp=0_dp
                yp=0_dp
                yperr=0_dp
        end if
        do nstp=1,MAXSTP
                call derivs(x,y,dydx)
                yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY ! at each step y(:) has error eps
!                yscal(:)=abs(h*dydx(:))+TINY ! each stepsize dy(:) has error eps
!                yscal(:)=abs(h*y(:))+abs(h*dydx(:))+TINY ! each stepsize dy(:) has error eps
!write(*,'(5E14.5)')yscal, h
                if (save_steps .and. (abs(x-xsav) > abs(dxsav))) &
                        call save_a_step
 !if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
 !if(h.gt.1.d10)h=sign(h1,x2-x1)
                do i = 1, size(GetExactList)
                   if( ((x+h-GetExactList(i))*(x+h-x1) > 0.0) .and. ((x-GetExactList(i))*(x-x1) < 0.0) )h=GetExactList(i)-x
!                   write(*,'(A,3ES19.10)')'Setting dead-on:',x,GetExactlist(i),x+h
                end do
! Set maximum stepsize for testing purposes:
!if((GetExactList(1).gt.0.d0).and.(h.gt.GetExactList(1)/1.d3))h=GetExactList(1)/1.d3
                call rkqs(y,dydx,x,h,eps,yscal, yerr,hdid,hnext,derivs)

!!$                if(size(y).eq.1)then
!!$                   write(*,*)x,y
!!$                end if
                
                if(OdeVoidError)return

                if(.not.nrsuccess)return

                if (hdid == h) then
                        nok=nok+1
                else
                        nbad=nbad+1
                end if
 ! WV NEW
                if(stopvar == 0) then
                   if ((x-x2)*(x2-x1) >= 0.0) then
                      ystart(:)=y(:)
                      if (save_steps) call save_a_step
                      RETURN
                   end if
                else
                   
                   if ((y(stopvar)-stopval)*(stopval-stopvarstart) >= 0.0) then
                      ystart(:)=y(:)
                      if (save_steps) call save_a_step
                      RETURN
                   end if
                end if
! END WV NEW
                
                if (abs(hnext) < hmin)then
                   call nrerror('stepsize smaller than minimum in odeint')
                   return
                end if
                h=hnext

        end do
        call nrerror('too many steps in odeint')
        return
        CONTAINS
!BL
        SUBROUTINE save_a_step
        kount=kount+1
        if (kount > size(xp)) then
                xp=>reallocate(xp,2*size(xp))
                yp=>reallocate(yp,size(yp,1),size(xp))
                yperr=>reallocate(yperr,size(yperr,1),size(xp))
        end if
        xp(kount)=x
        yp(:,kount)=y(:)
        yperr(:,kount)=yerr(:)
        xsav=x
        END SUBROUTINE save_a_step
        END SUBROUTINE odeint

        SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal, yerr,hdid,hnext,derivs)
        USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
        use ode_path, only : OdeVoidError
!        USE nr, ONLY : rkck
        IMPLICIT NONE

        REAL(dp), DIMENSION(:), INTENT(INOUT) :: y
        REAL(dp), DIMENSION(:), INTENT(IN) :: dydx,yscal
        REAL(dp), DIMENSION(:), INTENT(OUT) :: yerr
        REAL(dp), INTENT(INOUT) :: x
        REAL(dp), INTENT(IN) :: htry,eps
        REAL(dp), INTENT(OUT) :: hdid,hnext
        INTERFACE
                SUBROUTINE derivs(x,y,dydx)
                USE nrtype
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: x
                REAL(dp), DIMENSION(:), INTENT(IN) :: y
                REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx

                END SUBROUTINE derivs
        END INTERFACE
        INTEGER(I4B) :: ndum

        REAL(dp) :: errmax,h,htemp,xnew
        REAL(dp), DIMENSION(size(y)) :: ytemp
!        REAL(dp), DIMENSION(size(y)) :: yerr,ytemp
        REAL(dp), PARAMETER :: SAFETY=0.9_dp,PGROW=-0.2_dp,PSHRNK=-0.25_dp,&

                ERRCON=1.89e-4
        ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
        h=htry
        do
                call rkck(y,dydx,x,h,ytemp,yerr,derivs)
              if(OdeVoidError)return
               errmax=maxval(abs(yerr(:)/yscal(:)))/eps
! WV NEW: force smaller stepsize if derivs subroutine cannot handle numbers.
                if(any(dydx.gt.1.d29).or.any(isnan(dydx)).or.any(isnan(y)))then
                   errmax=10.0
                   !write(*,*)'Forcing arrived.1'
                end if
! END WV NEW
                if (errmax <= 1.0) exit
                htemp=SAFETY*h*(errmax**PSHRNK)

                h=sign(max(abs(htemp),0.1_dp*abs(h)),h)
                xnew=x+h
                if (xnew == x) then
                   call nrerror('stepsize underflow in rkqs')
                   return
                end if
        end do
        if (errmax > ERRCON) then
                hnext=SAFETY*h*(errmax**PGROW)
        else

                hnext=5.0_dp*h

        end if
        hdid=h
        x=x+h
        y(:)=ytemp(:)
        END SUBROUTINE rkqs

        SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
        USE nrtype; USE nrutil, ONLY : assert_eq
        use ode_path, only : OdeVoidError
        IMPLICIT NONE

        REAL(dp), DIMENSION(:), INTENT(IN) :: y,dydx
        REAL(dp), INTENT(IN) :: x,h
        REAL(dp), DIMENSION(:), INTENT(OUT) :: yout,yerr
        INTERFACE
                SUBROUTINE derivs(x,y,dydx)
                USE nrtype
                IMPLICIT NONE

                REAL(dp), INTENT(IN) :: x
                REAL(dp), DIMENSION(:), INTENT(IN) :: y
                REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx

                END SUBROUTINE derivs
        END INTERFACE
        INTEGER(I4B) :: ndum

        REAL(dp), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
        REAL(dp), PARAMETER :: A2=0.2_dp,A3=0.3_dp,A4=0.6_dp,A5=1.0_dp,&
                A6=0.875_dp,B21=0.2_dp,B31=3.0_dp/40.0_dp,B32=9.0_dp/40.0_dp,&
                B41=0.3_dp,B42=-0.9_dp,B43=1.2_dp,B51=-11.0_dp/54.0_dp,&
                B52=2.5_dp,B53=-70.0_dp/27.0_dp,B54=35.0_dp/27.0_dp,&
                B61=1631.0_dp/55296.0_dp,B62=175.0_dp/512.0_dp,&
                B63=575.0_dp/13824.0_dp,B64=44275.0_dp/110592.0_dp,&
                B65=253.0_dp/4096.0_dp,C1=37.0_dp/378.0_dp,&
                C3=250.0_dp/621.0_dp,C4=125.0_dp/594.0_dp,&
                C6=512.0_dp/1771.0_dp,DC1=C1-2825.0_dp/27648.0_dp,&
                DC3=C3-18575.0_dp/48384.0_dp,DC4=C4-13525.0_dp/55296.0_dp,&
                DC5=-277.0_dp/14336.0_dp,DC6=C6-0.25_dp

        ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
        ytemp=y+B21*h*dydx
        call derivs(x+A2*h,ytemp,ak2)
        if(any(dabs(ak2).gt.1.d29))then
           yerr=1.d30
           !write(*,*)'Forcing arrived.2'
           return
        end if
      if(OdeVoidError)return
        ytemp=y+h*(B31*dydx+B32*ak2)
        call derivs(x+A3*h,ytemp,ak3)
        if(any(dabs(ak3).gt.1.d29))then
           yerr=1.d30
           !write(*,*)'Forcing arrived.3'
           return
        end if
      if(OdeVoidError)return
        ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
        call derivs(x+A4*h,ytemp,ak4)
        if(any(dabs(ak4).gt.1.d29))then
           yerr=1.d30
           !write(*,*)'Forcing arrived.4'
           return
        end if
      if(OdeVoidError)return
        ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
        call derivs(x+A5*h,ytemp,ak5)
        if(any(dabs(ak5).gt.1.d29))then
           yerr=1.d30
           !write(*,*)'Forcing arrived.5'
           return
        end if
      if(OdeVoidError)return
        ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
        call derivs(x+A6*h,ytemp,ak6)
      if(OdeVoidError)return
        yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
        yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
! WV NEW: force smaller stepsize if derivs subroutine cannot handle numbers.
        if(any(dabs(ak2).gt.1.d29).or.any(dabs(ak3).gt.1.d29).or.any(dabs(ak4).gt.1.d29).or.any(dabs(ak5).gt.1.d29).or.any(dabs(ak6).gt.1.d29))then
           yerr=1.d30
           !write(*,*)'Forcing arrived.6'
        end if
! END WV NEW
        END SUBROUTINE rkck
end Module NROdeInt


