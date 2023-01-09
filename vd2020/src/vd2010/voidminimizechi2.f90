	subroutine void_nr_powell(p,xi,ftol,iter,fret,func)
!	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
!	USE nr, ONLY : linmin
	IMPLICIT NONE
	REAL(dl), DIMENSION(:), INTENT(INOUT) :: p
	REAL(dl), DIMENSION(:,:), INTENT(INOUT) :: xi
	INTEGER, INTENT(OUT) :: iter
	REAL(dl), INTENT(IN) :: ftol
	REAL(dl), INTENT(OUT) :: fret
	INTERFACE
		FUNCTION func(p)
			use precision
		IMPLICIT NONE
		REAL(dl), DIMENSION(:), INTENT(IN) :: p
		REAL(dl) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=20
	REAL(dl), PARAMETER :: TINY=1.0e-25 
	INTEGER :: i,ibig,n
	REAL(dl) :: del,fp,fptt,t
	REAL(dl), DIMENSION(size(p)) :: pt,ptt,xit
!  EXTERNAL :: VoidErrorReport
	n=size(p)
  if(.not.all(n==(/size(xi,1),size(xi,2)/)))stop 'size do not match in void_nr_powell.'
	fret=func(p)
	pt(:)=p(:)
	iter=0
	do
		iter=iter+1
		fp=fret
		ibig=0
		del=0.0
		do i=1,n
			xit(:)=xi(:,i)
			fptt=fret
			call void_nr_linmin(p,xit,fret,func)
			if (fptt-fret > del) then
				del=fptt-fret
				ibig=i
			end if
		end do
    !write(*,'(4ES14.5)'),fp,fret,2.0 *(fp-fret), ftol*(abs(fp)+abs(fret))+TINY
		if (2.0 *(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) RETURN
		if (iter == ITMAX) then
!      write(0,*)'powell exceeding maximum iterations',2.0 *(fp-fret), ftol*(abs(fp)+abs(fret))+TINY,2.0 *(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY
!      call VoidErrorReport()
return
!      stop
    end if 
		ptt(:)=2.0 *p(:)-pt(:)
		xit(:)=p(:)-pt(:)
		pt(:)=p(:)
		fptt=func(ptt)
		if (fptt >= fp) cycle
		t=2.0 *(fp-2.0 *fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
		if (t >= 0.0) cycle
		call void_nr_linmin(p,xit,fret,func)
		xi(:,ibig)=xi(:,n)
		xi(:,n)=xit(:)
	end do
	END subroutine void_nr_powell



	subroutine void_nr_linmin(p,xi,fret,func)
!	USE nrtype; USE nrutil, ONLY : assert_eq
!	USE nr, ONLY : mnbrak,brent
!	USE f1dim_mod
	IMPLICIT NONE
    INTERFACE
      FUNCTION func(x)
			use precision
!      USE nrtype
      REAL(dl), DIMENSION(:), INTENT(IN) :: x
      REAL(dl) :: func
      END FUNCTION func
    END INTERFACE
	INTEGER :: ncom
	REAL(dl), DIMENSION(:), POINTER :: pcom,xicom
	REAL(dl), INTENT(OUT) :: fret
	REAL(dl), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
	REAL(dl), PARAMETER :: TOL=1.0e-4 
	REAL(dl) :: ax,bx,fa,fb,fx,xmin,xx
!	ncom=assert_eq(size(p),size(xi),'linmin')
	ncom=size(p)
  if(ncom/=size(xi))then
    write(0,*)'sizes do not match in void_nr_linmin'
    stop
  end if
	pcom=>p
	xicom=>xi
	ax=0.0
	xx=1.0
	call void_nr_mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
	fret=void_nr_brent(ax,xx,bx,f1dim,TOL,xmin)
	xi=xmin*xi
	p=p+xi
  CONTAINS
  !BL
    FUNCTION f1dim(x)
    IMPLICIT NONE
    REAL(dl), INTENT(IN) :: x
    REAL(dl) :: f1dim
    REAL(dl), DIMENSION(:), ALLOCATABLE :: xt
    allocate(xt(ncom))
    xt(:)=pcom(:)+x*xicom(:)
    f1dim=func(xt)
    deallocate(xt)
    END FUNCTION f1dim
	END subroutine void_nr_linmin

	subroutine void_nr_mnbrak(ax,bx,cx,fa,fb,fc,func)
!	USE nrtype; USE nrutil, ONLY : swap
	IMPLICIT NONE
	REAL(dl), INTENT(INOUT) :: ax,bx
	REAL(dl), INTENT(OUT) :: cx,fa,fb,fc
	INTERFACE
		FUNCTION func(x)
			use precision
!		USE nrtype
		IMPLICIT NONE
		REAL(dl), INTENT(IN) :: x
		REAL(dl) :: func
		END FUNCTION func
	END INTERFACE
	REAL(dl), PARAMETER :: GOLD=1.618034 ,GLIMIT=100.0 ,TINY=1.0e-20 
	REAL(dl) :: fu,q,r,u,ulim, xswap
	fa=func(ax)
	fb=func(bx)
	if (fb > fa) then
		!call void_nr_swap(ax,bx)
    !
    xswap = ax
    ax = bx
    bx = xswap
    xswap = fa
    fa = fb
    fb = xswap
!		call void_nr_swap(fa,fb)
	end if
	cx=bx+GOLD*(bx-ax)
	fc=func(cx)
	do
		if (fb < fc) RETURN
		r=(bx-ax)*(fb-fc)
		q=(bx-cx)*(fb-fa)
		u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0 *sign(max(abs(q-r),TINY),q-r))
		ulim=bx+GLIMIT*(cx-bx)
		if ((bx-u)*(u-cx) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				ax=bx
				fa=fb
				bx=u
				fb=fu
				RETURN
			else if (fu > fb) then
				cx=u
				fc=fu
				RETURN
			end if
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		else if ((cx-u)*(u-ulim) > 0.0) then
			fu=func(u)
			if (fu < fc) then
				bx=cx
				cx=u
				u=cx+GOLD*(cx-bx)
				call void_nr_shft(fb,fc,fu,func(u))
			end if
		else if ((u-ulim)*(ulim-cx) >= 0.0) then
			u=ulim
			fu=func(u)
		else
			u=cx+GOLD*(cx-bx)
			fu=func(u)
		end if
		call void_nr_shft(ax,bx,cx,u)
		call void_nr_shft(fa,fb,fc,fu)
	end do
	CONTAINS
!BL
	subroutine void_nr_shft(a,b,c,d)
	REAL(dl), INTENT(OUT) :: a
	REAL(dl), INTENT(INOUT) :: b,c
	REAL(dl), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END subroutine void_nr_shft
	END subroutine void_nr_mnbrak

	FUNCTION void_nr_brent(ax,bx,cx,func,tol,xmin) result(brent)
!	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(dl), INTENT(IN) :: ax,bx,cx,tol
	REAL(dl), INTENT(OUT) :: xmin
	REAL(dl) :: brent
	INTERFACE
		FUNCTION func(x)
			use precision
!		USE nrtype
		IMPLICIT NONE
		REAL(dl), INTENT(IN) :: x
		REAL(dl) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER, PARAMETER :: ITMAX=100
	REAL(dl), PARAMETER :: CGOLD=0.3819660 ,ZEPS=1.0e-3 *epsilon(ax)
	INTEGER :: iter
	REAL(dl) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.0
	fx=func(x)
	fv=fx
	fw=fx
	do iter=1,ITMAX
		xm=0.5 *(a+b)
		tol1=tol*abs(x)+ZEPS
		tol2=2.0 *tol1
		if (abs(x-xm) <= (tol2-0.5 *(b-a))) then
			xmin=x
			brent=fx
			RETURN
		end if
		if (abs(e) > tol1) then
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0 *(q-r)
			if (q > 0.0) p=-p
			q=abs(q)
			etemp=e
			e=d
			if (abs(p) >= abs(0.5 *q*etemp) .or. &
				p <= q*(a-x) .or. p >= q*(b-x)) then
				e=merge(a-x,b-x, x >= xm )
				d=CGOLD*e
			else
				d=p/q
				u=x+d
				if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
			end if
		else
			e=merge(a-x,b-x, x >= xm )
			d=CGOLD*e
		end if
		u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
		fu=func(u)
		if (fu <= fx) then
			if (u >= x) then
				a=x
			else
				b=x
			end if
			call shft(v,w,x,u)
			call shft(fv,fw,fx,fu)
		else
			if (u < x) then
				a=u
			else
				b=u
			end if
			if (fu <= fw .or. w == x) then
				v=w
				fv=fw
				w=u
				fw=fu
			else if (fu <= fv .or. v == x .or. v == w) then
				v=u
				fv=fu
			end if
		end if
	end do
	write(0,*)'void_nr_brent: exceed maximum iterations'
!  call VoirErrorReport()
  stop
	CONTAINS
!BL
	SUBROUTINE shft(a,b,c,d)
	REAL(dl), INTENT(OUT) :: a
	REAL(dl), INTENT(INOUT) :: b,c
	REAL(dl), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END FUNCTION void_nr_brent
