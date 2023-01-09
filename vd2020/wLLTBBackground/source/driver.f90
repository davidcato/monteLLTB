program wLLTBTestDrive
  use wlltb_types
  use wlltb_constants
  use wlltb_kofr
  use wlltb
  implicit none
  
  type(AllCosmology) :: ac
  real(dl) :: H0_inf, Lambda, w, Mtilde, &
						r, t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp, tturn
  
  real(dl) :: Rltbw0, Rpw0, Rpdw0, Rpddw0, Sw0, Sdw0, Sddw0, aw0, apw0, apdw0, addw0, apddw0, Hw0, Hpw0, tturnw0
  real(dl) :: omegak, omegal, tbar
  real(dl) :: Hlocal,delta0,t0, Rltbplusdr, dr

  real(dl) :: preva, prevt, prevadd

  integer, parameter :: tsteps = 200, rsteps = 200, wsteps=5, lsteps=5, ksteps=2
  integer :: ti, ri, rcount, wi, li, ki, wimax,limax,kimax, tistart, ristart, rimax, thetime(3)
  real(dl) :: wmin, wmax, lmin, lmax, kmin, kmax
  real(dl) :: FMbar, For3
  logical :: testlltb = .true.
  
  character(len=64) :: filename
  
  wmin = -2.0d0
  wmax = -0.4d0

  lmin = 0.1d0
  lmax = 0.9d0
  
  kmin = -10.d0
  kmax = 20.d0
  
  write(*,*) "Hello world"


  ac%kofr_att%a = 4.d0
  ac%kofr_att%b = 2.d0
  ac%kofr_att%kmax = -10.d0
  ac%kofr_att%L = 1.d-2 ! units H0^-1, that is Gyr -> Glyr -> Gpc / 3


  ac%kofr => fk3
  ac%dkdr => fk3prime

  H0_inf = 70.d0

  if(testlltb)then
    wimax = 1
    limax = 1
    kimax = 1
    tistart = tsteps
!    ristart = 10
!    rimax = 10
    ristart = 1
    rimax = rsteps
  else
    wimax = wsteps
    limax = lsteps
    kimax = ksteps
    tistart = tsteps
    ristart = 1
    rimax = rsteps
  end if
  
  do wi = 1, wimax
  do li = 1, limax
  do ki = 1, kimax
  
call system_clock(count=thetime(1),count_rate=thetime(3))

  Omegal = 0.7d0
  omegak = 0.d0
  w = -1.2d0
ac%kofr_att%kmax = kmax
  if(.not.testlltb)then
    omegal = dble(li-1)/dble(lsteps-1) * (lmax-lmin) + lmin
    w = dble(wi-1)/dble(wsteps-1) * (wmax-wmin) + wmin
    ac%kofr_att%kmax = dble(ki-1)/dble(ksteps-1) * (kmax-kmin) + kmin
  end if  


  Lambda = 3._dl * H0_inf**2 * OmegaL
  ac%kofr_att%kb = 4.d0 * pi / 3.d0 * ( Omegak / (1.d0 - Omegak - OmegaL) )
  a=0
  
  call lltb_normalization(H0_inf,ac%kofr_att%kb,ac%kofr_att%kmax,Lambda,w,Mtilde,Hlocal,delta0,t0,tbar)
!write(*,*)delta0
!stop
  if(.not.testlltb)then
    write(filename,'("output/testwltb_L",F3.1,"_w",F4.1,"_dm",F7.5,".tsv")')omegal,w,delta0
    write(*,'(A)')filename
  !  stop

    open(16,file=filename,status="REPLACE")
    write(16,'(A)')"# H0 t, H0 r, delta_m(r), H0 Y(r,t), Y'(r,t), w, H0/100.d0, OmegaLambda, F(r,t)/FM(r,tbar)"
  end if

  t=tbar
  dr=ac%kofr_att%L*1.d-4
  do ti = tistart, tsteps

      prevt = t
      t = tbar + dble(ti - 1)/dble(tsteps - 1) * (t0 - tbar)

      rcount = 0
!      do ri = 10,10
      do ri = ristart, rimax
        rcount = rcount + 1
        r = ac%kofr_att%L * dble(ri - 1)/dble(rsteps - 1) *1.1

  
!  Mtilde =  sqrt(3._dl * (H0_inf**2 - Lambda/3._dl)/8._dl/pi/(1._dl+0.75_dl * ac%kofr_att%kb/pi))
!  r=1.d0
!  t=t0
!  write(*,*)'T0:',t0
        preva = a*H
        prevadd = add
        call lltb_functions(H0_inf, Lambda, w, mykofr, mydkdr, mytofr, mytofr, Mtilde, &
						r, t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp)!, tturn, For3)

!        call lltb_functions(H0_inf, Lambda, w, mykofr, mydkdr, mytofr, mytofr, Mtilde, &
!						r+dr, t, &
!						Rltbplusdr)
        if(testlltb)then
          call lltb_functions(H0_inf, Lambda, -1.d0, mykofr, mydkdr, mytofr, mytofr, Mtilde, &
						r, t, &
						Rltbw0, Rpw0, Rpdw0, Rpddw0, Sw0, Sdw0, Sddw0, aw0, apw0, apdw0, addw0, apddw0, Hw0, Hpw0, tturnw0)
                              !  1,2,3          ,4        ,5          ,6            ,7      ,8       
          write(16,'(100ES14.5)')t,r,Rltb/Rltbw0, Rp/ Rpw0, Rpd/ Rpdw0, Rpdd/ Rpddw0, S/ Sw0, Sd/ Sdw0, &
              Sdd/ Sddw0 & !9   ,10     ,11       ,12         ,13         ,14           ,15     ,16       ,17   
                                ,a/ aw0 , ap/ apw0, apd/ apdw0, add/ addw0, apdd/ apddw0, H/ Hw0, Hp/ Hpw0, tturn/tturnw0
!          write(17,*)t,apdd,apddw0
        end if

        if(.not.testlltb)then
          call lltb_normalization(H0_inf,ac%kofr_att%kb,ac%kofr(r,ac%kofr_att),Lambda,w,Mtilde,Hlocal,delta0,t0,tbar)
          FMbar = delta0 + 1.d0
          write(16,'(100ES14.5)')H0_inf*t, H0_inf*r, delta0, H0_inf*Rltb, Rp, w, H0_inf/100.d0, omegal, For3 / H0_inf**2
          call lltb_normalization(H0_inf,ac%kofr_att%kb,ac%kofr_att%kmax,Lambda,w,Mtilde,Hlocal,delta0,t0,tbar)
        end if

!        if(ti>1)write(19,*)t,0.5*(prevadd + add),(a*H-preva)/(t-prevt)
        
!        write(21,*)t,(Rltbplusdr - Rltb)/dr, Rpw0
        
if(.false.)then
if(ri==1)then
  For3 = H
 call lltb_functions(H0_inf, Lambda, w, mykofr, mydkdr, mytofr, mytofr, Mtilde, &
						ac%kofr_att%L*1.1, t, &
						Rltb, Rp, Rpd, Rpdd, S, Sd, Sdd, a, ap, apd, add, apdd, H, Hp, tturn)
  write(*,*)For3,For3/H,For3/H0_inf, H, H0_inf, a       
end if
end if


      end do
  end do

      close(16)
  
call system_clock(count=thetime(2),count_rate=thetime(3))
write(*,'(A,F7.3,A)')'Spent ',dble(thetime(2)-thetime(1))/dble(thetime(3)),' seconds in one cosmology.'

  end do
  end do
  end do
  

contains

  function mykofr(r)
    implicit none
    real(dl), intent(in) :: r
    real(dl) :: mykofr
    
    mykofr = ac%kofr(r,ac%kofr_att)
    
  end function mykofr

  function mydkdr(r)
    implicit none
    real(dl), intent(in) :: r
    real(dl) :: mydkdr
    
    mydkdr = ac%dkdr(r,ac%kofr_att)

  end function mydkdr

  function mytofr(r)
    implicit none
    real(dl), intent(in) :: r
    real(dl) :: mytofr
    
    mytofr = 0.d0
  end function mytofr

end program wLLTBTestDrive