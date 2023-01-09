  ! New VDII:
  function VoidAgeGYR()
    implicit none
    real(dl) :: VoidAgeGYR
    
    VoidAgeGYR = VP%t0_gyr
    
  end function VoidAgeGYR


!  subroutine VoidGetThetaOut(CMB)
!    type(FakeCMBParams) :: CMB
!    CMB%VCP%theta = voidTh%out
!  end subroutine VoidGetThetaOut

!  subroutine VoidSaveTheta(CMB,theta)
!    type(FakeCMBParams) :: CMB
!    real, intent(in) :: theta
!    if(CMB%VCP%invoid)VoidTh%in = theta
!!    if(VP%zB.le.0 .and. .not.CMB%VCP%invoid)then ! When we are reset
!    if(dabs(VP%zB).lt.1.d-20 .and. .not.CMB%VCP%invoid)then ! When we are reset
!       VoidTh%out = theta
!    else if(.not.CMB%VCP%invoid)then
!       VoidTh%obs = theta
!    end if
!
!!!$   if(CMB%VCP%invoid)write(*,*)'theta-in set.',VoidTh%in
!!!$   if(VP%zB.le.0 .and. .not.CMB%VCP%invoid)then
!!!$      write(*,*)'theta-out set.',VoidTh%out
!!!$    else if(.not.CMB%VCP%invoid)then
!!!$       write(*,*)'theta-obs set.',VoidTh%obs
!!!$    end if
!
!  end subroutine VoidSaveTheta


  function VoidHlocal()
    real(dl) :: VoidHlocal

!    if(VP%zB > 0.d0)then
!    if(dabs(VP%zB) > 1.d-20)then
    if(.not.VP%is_pure_flrw)then
       VoidHlocal = VP%Hlocal_paper
    else ! we've been reset
       voidHlocal = -1.d0
    end if
       
    
  end function VoidHlocal

    function t0_effective()
    real(dl) :: t0_effective
    !real(dl) :: zboundary
    real(dl) :: r

  ! first get (r,t) for FLRW equivalent observer:
    !r = VP%r0 !0.d0
    r = VP%r0 !0.d0
    t0_effective = NT_flrw(r = r)

  end function t0_effective

  function VoidTVOverTF(tCOMP)
!    real(dl) :: t
    real(dl) :: VoidTVOverTF
    real(dl), optional :: tCOMP
    character(len=64) :: thisT
    character(len=64), save :: prevT

    real(dl) :: a0_f_obs, aL_f, aL_v, tl!, tf
    real(dl) :: daun, vr, tnr

!    if(voidfeedback.gt.0)write(*,*)'Getting TVoverTF.'
    if(VoidError.or..not.associated(VP%LTB_FLRW_r))then
       VoidTVOverTF=1.d0
       if(voidfeedback.gt.0)write(*,*)'Skipping TVoverTF.'
       return
    end if
!    if(VP%zB.eq.0.d0)then
!    if(dabs(VP%zB)<1.d-20)then
    if(VP%is_pure_flrw)then
       VoidTVOverTF=1.d0
    else
       
       ! Matching is at t=tl
       if(.not.present(tCOMP))then
        tl = nt(r=VP%LTB_FLRW_r)
       else
        tl = tCOMP
       end if
       aL_v = 1.d0 / (nz(t=tl) + 1.d0)
       aL_f = 1.d0 / (nz_FLRW(t=tl) + 1.d0)
       ! FLRW observer is always at r=0.
       ! This is fine because r = irrelevant in FLRW.
       ! See voidasymp.f90.
       a0_f_obs =  1.d0 / (nz_FLRW(0.d0) + 1.d0)
!!$       aL_v = 1.d0 / (nz(VP%LTB_FLRW_r) + 1.d0)
!!$       aL_f = 1.d0 / (nz_FLRW(VP%LTB_FLRW_r) + 1.d0)
!!$       a0_f_obs =  1.d0 / (nz_FLRW(0.d0) + 1.d0)
!       a0_f_obs =  1.d0 / (nz_FLRW(VP%r0) + 1.d0)
       VoidTVOverTF = a0_f_obs / aL_f * aL_v
       !(1.d0 + nz_FLRWt(t)) / (1.d0 + NZ(VP%LTB_FLRW_r))
    end if
    

    if(voidfeedback.gt.1)then
      thisT=""
      write(thisT,'(A,4E14.5)')' # VoidTVOverTF - 1:',VoidTVOverTF-1.d0
      if(thisT/=prevT)then
        write(*,*)trim(thisT)
        prevT = thisT

        if(voidfeedback>2)then
           write(*,*)'Double check DA to t(r=L). r''s should differ if r_obs/=0, DA''s should equal. L=',VP%L
daun = DA_unnorm(0.d0,(/tl,VP%LTB_FLRW_r/))
vr = voidR(t=tl,r=nr_FLRW(t=tl))
tnr = nr_FLRW(t=tl)
           write(*,'(3(A,ES24.15))')' DA:',daun,' r:',VP%LTB_FLRW_r,' obs: LTB'
           al_v = VP%L
           VP%L = 0.d0
           write(*,'(3(A,ES24.15))')' DA:',vr,' r:',tnr,' obs: FLRW'
           VP%L = al_v
        end if

      end if
    end if
!    if(voidfeedback.gt.0)write(*,*)'Done TVoverTF.'

  end function VoidTVOverTF

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!
!
!      Make DA
!
!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

!
!  subroutine VoidDzArray(VDI)
!!    use VoidDistSNChi2
!    use ode_path, only : nrsuccess
!!    use camb
!    implicit none
!    type(VoidDistInfo), pointer :: vdi
!    real(dl) :: Hrat, H0kms,zB
!!    real(dl) :: z, dum
!!    real(dl) :: dum
!    integer :: i!, mystat
!
!    if(voidfeedback.gt.0)write(*,*)'Making DA(SN)..'
!
!    zB=VP%zb
!    Hrat=VP%Hlocal/VP%H0
!    H0kms=VP%H0
!
!!call testmakedz
!
!!    if(zB.gt.0.d0)then
!    call MakeDz(vdi)
!
!    if(voidfeedback.gt.0)write(*,*)'Made DA(SN).'
!
!    if( .not. nrsuccess) VoidError=.true.
!
!  end subroutine VoidDzArray
!
!
!  subroutine testMakeDz
!!    use camb
!    use VoidDistSNChi2
!    implicit none
!    real(dl), allocatable :: zz(:), Dzdz(:)
!    real(dl) :: Thish,ThisHrat,ThiszB
!    integer :: i, imax
!    type(VoidDistInfo), target :: vdi
!    type(VoidDistInfo), pointer :: tvdi
!
!    !----------------------------------------!
!    ! In order to test da(flrw) uncomment:
!    ! call init_flrw
!    ! VP%L=0.d0
!    !----------------------------------------!
!
!    write(*,*) 'Testing MakeDz'
!    imax = 250
!    allocate(zz(imax))
!    allocate(Dzdz(imax))
!    thiszb = 2.!nz_flrw(r=VP%L)
!    do i=1,imax
!!       zz(i) = dble(i)/dble(imax)*VP%zB*1.25
!       zz(i) = dble(i)/dble(imax)*thiszb*1.25
!    end do
!
!    tvdi => vdi
!
!    call IniMyVDInfo(tvdi,imax)
!
!!!$    allocate(vdi%snz(imax))
!!!$    allocate(vdi%nda(imax))
!!!$    allocate(vdi%da(1,1))
!!!$    allocate(vdi%weights(1,1))
!
!
!    vdi%snz = zz
!
!    call MakeDz(tvdi)
!
!
!    ThisHrat = VP%Hlocal / VP%H0
!    Thish = VP%H0/100.d0
!!    ThiszB = VP%zB
!    open(123,file='testda_z.txt')
!    do i=1,imax
!       write(123,'(4E14.5)')zz(i),vdi%da(i,1)
!       write(*,'(4E14.5)')zz(i),vdi%da(i,1)
!    end do
!    close(123)
!    deallocate(zz)
!    deallocate(Dzdz)
!    deallocate(vdi%snz, vdi%nda, vdi%da)
!    thiszB=VoidTVOverTF()
!    write(*,*)'Hlocal, voidtvovertf:',VP%Hlocal_paper,thiszB
!    write(*,'(A,3E17.8)')'L, zB, r(zB):',VP%L,VP%zB,NR(VP%zB)
!    write(*,'(A,3E17.8)')'VP%kb_om:',VP%kb_om
!    write(*,'(A,3E17.8)')'VP%Mtilde:',VP%Mtilde
!    stop 'Tested makedz? See testda_z.txt'
!
!  end subroutine testMakeDz



!  subroutine MakeDz(vdi)
!!    use camb
!    use VoidDistSNChi2
!    implicit none
!    type(voiddistinfo), pointer :: vdi
!    real(dl), pointer :: thisda(:), thesedas(:,:), thisweight(:), theseweights(:,:)
!!    real(dl) :: deltachi2, thisdelta
!    integer :: i, znum, maxn, mystat
!
!    if(maxval(vdi%snz).gt.void_zmax)then
!       stop 'Increase void_maxz in voiddistances.f90 to match at least the largest z in your SN sample.'
!    end if
!
!    znum = size(vdi%snz)
!
!    if(VoidError)then
!
!       call IniMyVData(vdi,znum,1)
!       vdi%da = 0.d0
!       vdi%weights = 1.d0
!       vdi%mu = 0.d0
!       vdi%diffs = 1.d30
!       vdi%nda = 1
!       return
!    end if
!
!
!    allocate(thesedas(znum,10))
!    allocate(theseweights(znum,10))
!
!    thesedas=1.d0
!    theseweights=0.d0
!
!
!
!    allocate(thisda(1))
!    allocate(thisweight(1))
!    do i=1,znum
!       call DegenDA(vdi%snz(i),thisda,thisweight)
!       vdi%nda(i) = size(thisda)
!       thesedas(i,1:vdi%nda(i))=thisda
!       theseweights(i,1:vdi%nda(i))=thisweight
!
!    end do
!
!    if(VoidError)then
!
!       call IniMyVData(vdi,znum,1)
!       vdi%da = 0.d0
!       vdi%weights = 1.d0
!       vdi%mu = 0.d0
!       vdi%diffs = 1.d30
!       vdi%nda = 1
!       return
!    end if
!
!    maxn = maxval(vdi%nda)
!    call IniMyVData(vdi,znum,maxn)
!    vdi%da = thesedas(:,1:maxn)
!    vdi%weights = theseweights(:,1:maxn)
!
!    call VoidMakeMuDiff(vdi,maxn)
!
!
!
!    deallocate(thesedas, theseweights)
!    deallocate(thisDA,stat=mystat)
!    deallocate(thisweight,stat=mystat)
!
!  end subroutine MakeDz

!  subroutine VoidMakeMuDiff(vdi,maxn)
!    use VoidDistSNChi2
!    implicit none
!    type(VoidDistInfo), pointer :: vdi
!    integer :: maxn, i
!    ! wv new simplify marginalization:
!    real(dl) :: avv_diff
!
!    do i = 1, maxn
!       vdi%mu(:,i)=5.d0*log10((1.d0+vdi%snz(:))**2*vdi%DA(:,i))+25.d0 
!       vdi%diffs(:,i)= vdi%mu(:,i) - vdi%moduli(:)
!    end do
!
!    ! cheat our way into marginalization over H0: normalize to average diff,
!    ! such that best fit range {minH < H0 < maxH} is good enough not to be 
!    ! rejected a priori by our integrator.
!    avv_diff=0.d0
!    do i = 1, maxn
!       avv_diff=avv_diff + sum(vdi%diffs(:,i))
!    end do
!    avv_diff = avv_diff / dble(maxn) / dble(size(vdi%snz))
!    
!    vdi%diffs=vdi%diffs - avv_diff 
!    vdi%mu=vdi%mu-avv_diff
!    vdi%avvdiff = avv_diff
!
!!!$do i = 1, size(vdi%snz)
!!!$   write(1234,*)vdi%snz(i),vdi%mu(i,1)
!!!$end do
!!!$stop 'see fort.1234'
!  end subroutine VoidMakeMuDiff

  subroutine DegenDA(z,thisDA, DAweight)
    implicit none
    real(dl), pointer :: tr(:,:), thisDA(:), DAweight(:)
    real(dl) :: z
    integer :: mystat, i, n
    
    deallocate(thisDA,stat=mystat)
    deallocate(DAweight,stat=mystat)

    if(VoidError)then
       allocate(thisDA(1))
       allocate(DAweight(1))
       thisDA=1.d30
       DAweight = 0.d0
       return
    end if
    
    allocate(tr(2,1))
!    call masplint(VPNR%ztr(1,:),VPNR%ztr(2:3,:),z,tr,VPNR%ddztr(2:3,:),VPNR%kmaxz,VPNR%kminz)
    ! new 23/02/2011
    call MultiValuedInvert(VPNR,z,tr)

    if(Lbound(tr,1)>1)then
       call VoidErrorReport()
       stop 'It happened: lbound(tr:,1)>1. In Degenda.'
    end if
    n = size(tr(1,:))
    allocate(thisDA(n),DAweight(n))

    do i = 1, n
       thisDA(i) = DA(0.d0,tr(:,i))
       ! Define weighting of each value here:
       ! Flat weighting:
       DAweight(i) = 1.d0 / dble(n)
       !write(*,*)'degenda',thisda(i)
    end do

    deallocate(tr, stat=mystat)
    return
  end subroutine DegenDA

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!
!
!    Make CMBin, CMBout, and CMBobs
!
!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

  subroutine EffectiveParams(CMB)
    type(FakeCMBParams), intent(in) :: CMB
    
    real(dl) :: T0cmb, N_ur, T0cmb_eff
    real(dl) :: H0_out, H0_eff,h_eff
    real(dl) :: t0eff, tL, zL, zL_eff, magical_z
    real(dl) :: OmegaG_eff, OmegaUR_eff, OmegaK_eff, OmegaB_eff, OmegaDM_eff !,OmegaDE_eff lets CLASS takes care of OmegaDE
    real(dl) :: omegah2_b_eff, omegah2_dm_eff

    t0eff = t0_effective()

    tL = nt(r=VP%LTB_FLRW_r)  ! or tl, time at z = zb
    
    zL = nz(t=tL) ! zL is effectively direfferent to zb
    zL_eff = nz_FLRW(t=tL)
    magical_z = (1+zL) / (1+zL_eff)

    !Now real part to effective parameters
    H0_out = voidHaout(VP%t0bg)
    H0_eff = voidHaout(t0eff)
    h_eff = H0_eff /100.d0

    N_ur=3.046d0

    T0cmb=2.7255d0
    T0cmb_eff=magical_z*T0cmb
    OmegaG_eff=2.469e-5_dl*(T0cmb_eff/T0cmb)**4/h_eff**2
    OmegaUR_eff=(7.d0/8.d0)*N_ur*(4.d0/11.d0)**(4./3.)*OmegaG_eff
    OmegaK_eff=(H0_out/H0_eff)**2 * (voidaout(VP%t0bg)/voidaout(t0eff))**2 * CMB%omk
    OmegaB_eff=CMB%omegabh2*magical_z**3/h_eff**2
    OmegaDM_eff=CMB%omegadmh2*magical_z**3/h_eff**2

    omegah2_b_eff =CMB%omegabh2*magical_z**3
    omegah2_dm_eff =CMB%omegadmh2*magical_z**3

    
    write(*,*) 'EffectiveQuantities:'
    !write(*, '(100A20)') 't0_eff', 'zL', 'zL_eff', 'magical_z', 'H0_out', 'H0_eff', 'Tcmb_eff'
    !write(*, '(100E20.11)') t0eff, zL, zL_eff, magical_z, H0_out, H0_eff, T0cmb_eff
    write(*, '(100A20)') 't0_eff', 'zL', 'zL_eff', 'magical_z', 'H0_out'
    write(*, '(100E20.11)') t0eff, zL, zL_eff, magical_z, H0_out

    write(*,*) 'EffectiveParams:'
    !write(*, '(100A20)') 'OmegaG_eff', 'OmegaUR_eff', 'OmegaK_eff', 'OmegaB_eff', 'OmegaDM_eff', 'omegah2_b_eff', 'omegah2_dm_eff'
    !write(*, '(100E20.11)') OmegaG_eff, OmegaUR_eff, OmegaK_eff, OmegaB_eff, OmegaDM_eff, omegah2_b_eff, omegah2_dm_eff
    write(*, '(100A20)') 'H0_eff', 'Tcmb_eff', 'OmegaK_eff', 'OmegaB_eff', 'OmegaDM_eff', 'omegah2_b_eff', 'omegah2_dm_eff'
    write(*, '(100E20.11)') H0_eff, T0cmb_eff, OmegaK_eff, OmegaB_eff, OmegaDM_eff, omegah2_b_eff, omegah2_dm_eff

  end subroutine EffectiveParams

  subroutine VoidConvert(CMBout,CMBin, CMBobs,z)!, CambInside)
    type(FakeCMBParams), intent(in) :: CMBout
    type(FakeCMBParams), intent(out) :: CMBin, CMBobs
    real(dl), optional :: z
!    logical, optional :: CambInside
    
    ! ----------------- !
    ! To copy all other parameters:

    CMBin = CMBout
    CMBobs = CMBout

    if(voiderror)return
    ! ----------------- !

    ! These are only relevant after evrything is done, for output:
    CMBin%vcp%Theta = VoidTh%in
    CMBobs%vcp%Theta = VoidTh%obs

!    if(dabs(CMBout%vcp%zb) .lt. 1.d-20)then
    if(VP%is_pure_flrw)then
       return
    end if

    if(present(z))then
       call VoidMakeCMBin(CMBout,CMBin, z)
    else
       call VoidMakeCMBin(CMBout,CMBin)
    end if

    call VoidMakeCMBobs(CMBout,CMBobs)

    if(voidtesting .or. VoidFeedback>1)then
       write(*,*)'Omega_dm h2, out:',CMBout%omdmh2
       write(*,*)'Omega_b h2, out:',CMBout%ombh2
       write(*,*)'Omega_m, out:',(CMBout%omdmh2+CMBout%ombh2)/CMBout%h0**2*1.d4
       write(*,*)'Omega_k, out:',CMBout%omk
       write(*,*)'Omega_v, out:',CMBout%omv
    end if

!!$    call VoidTestCMB(CMBobs)
!!$    if(present(CambInside))then
!!$       if(CambInside)then
!!$          call VoidTestCMB(CMBin)
!!$       end if
!!$    end if
!!$
!!$    call VoidBigAlert('Exchanging inner and outer CMB!!.')
!!$    CMBobs = CMBin

    ! set parameter priors for cosmology
    call VoidCMBParamsPriors(CMBout,CMBin,CMBobs)
    if(voiderror)return

  end subroutine VoidConvert

  
  ! Apply all possible omega-priors here:
  ! e.g. OmegaM > 0
  !      OmegaL > -1
  !      OmegaK > -1
  ! etc
  subroutine VoidCMBParamsPriors(CMBout,CMBin,CMBobs)
    implicit none
    type(FakeCMBParams), intent(in) :: CMBout, CMBin, CMBobs

    if(any((/CMBout%omdmh2,CMBin%omdmh2,CMBobs%omdmh2/)<1.d-20))VoidError=.true.
    
    if(CMBobs%omv < -1.d0)VoidError=.true.

  end subroutine VoidCMBParamsPriors

  subroutine VoidMakeCMBobs(CMBout,CMBobs)
    type(FakeCMBParams), intent(in) :: CMBout
    type(FakeCMBParams), intent(inout) :: CMBobs
    
    real(dl) :: r,t, a0L, aF
    real(dl) :: H2obs, Omm_out, omk_out, omdm_out, omb_out
    real(dl) :: omm_obs, omk_obs, omdm_obs, omb_obs, omdmh2_obs , ombh2_obs 

    ! but:
    CMBobs%vcp%invoid = .false.
    CMBobs%vcp%fromMCMC = .false.
    CMBobs%vcp%effobs = .true.
    CMBobs%vcp%outside = .false.
    
    if(VoidSplitISW)then
       CMBobs%vcp%splitISW = .true.
       CMBobs%vcp%withISW = .false.
    end if

    CMBobs%vcp%zb = VP%zb
    ! first get (r,t) for FLRW equivalent observer:
    r = VP%r0 !0.d0
    t = NT_flrw(r = r)

    ! get scale factors:
    a0L = voida(r=VP%LTB_FLRW_r, t = VP%t0)
    aF = voida(r=VP%LTB_FLRW_r,t=t)


    omm_out = (CMBout%omdmh2 + CMBout%ombh2) / (CMBout%H0/100.d0)**2
    omk_out = CMBout%omk

    omdm_out = (CMBout%omdmh2) / (CMBout%H0/100.d0)**2
    omb_out = (CMBout%ombh2) / (CMBout%H0/100.d0)**2

    H2obs = CMBout%H0**2 * (  omm_out * (a0L/aF)**3 +  omk_out * (a0L/aF)**2  + CMBout%omv*(a0L/aF)**(3*(1+VP%w)) )
    
    omm_obs =   CMBout%H0**2 / H2obs  * (  omm_out * (a0L/aF)**3 )
    omk_obs =   CMBout%H0**2 / H2obs  * (  omk_out * (a0L/aF)**2 )

    omdm_obs = omdm_out / omm_out * omm_obs
    omb_obs = omb_out / omm_out * omm_obs

    omdmh2_obs = omdm_obs * H2obs / 1.d4
    ombh2_obs = omb_obs  * H2obs / 1.d4 

    CMBobs%H0 = dsqrt(H2obs)
    CMBobs%vcp%hlocal = VP%Hlocal_paper
    CMBobs%vcp%TVoverTF = VoidTVoverTF()
    CMBobs%omdmh2 = omdmh2_obs
    CMBobs%ombh2 = ombh2_obs
    CMBobs%omk = omk_obs
!write(*,*)'test this line (make obs):'
CMBobs%omv = 1._dl - omm_obs - omk_obs
if(voidtesting)then
CMBobs%omv = 1._dl - omm_obs - omk_obs
write(*,*)CMBobs%omv,CMBout%H0**2 / H2obs * CMBout%omv* (  (a0L/aF)**(3*(1+VP%w)) ), CMBout%omv * CMBout%H0**2 / H2obs, CMBout%omv
!stop 'make cmbobs'
end if

    if(voidtesting .or. VoidFeedback>1)then
       write(*,*)'Omega_m, obs:',omm_obs
       write(*,*)'Omega_k, obs:',omk_obs
       write(*,*)'Omega_v, obs:',CMBobs%omv
       write(*,*)'H(tf), obs:',dsqrt(H2obs)
       write(*,*)'Omega_m, obs + Omega_k, obs + Omega_v, obs:',omm_obs  + omk_obs + CMBobs%omv
        if(voidtesting .or. VoidFeedback>2)then
           write(*,*)'Omega_dm h2, obs:',omdmh2_obs
           write(*,*)'Omega_b h2, obs:',ombh2_obs
           write(*,*)'H(t0):',CMBout%H0
           write(*,*)'H(tf), obs / H(t0):',dsqrt(H2obs)/ CMBout%H0
           write(*,*)'Omega_m, obs + Omega_k, obs:',omm_obs  + omk_obs
        end if
    end if
!!$     
!!$    CMBobs%H0 = dsqrt(H2obs)
!!$    CMBobs%vcp%hlocal = VP%Hlocal_paper
!!$    CMBobs%vcp%TVoverTF = VoidTVoverTF()
!!$    CMBobs%omdmh2 = omdmh2_obs
!!$    CMBobs%ombh2 = ombh2_obs
!!$    CMBobs%omk = omk_obs
!!$write(*,*)'test this line (make obs):'
!!$CMBobs%omv = 1._dl - omm_obs - omk_obs
    
    CMBobs%vcp%t0 = VP%t0bg * 977.8d0

  end subroutine VoidMakeCMBobs



  subroutine VoidMakeCMBin(CMBout,CMBin, z)
    type(FakeCMBParams), intent(in) :: CMBout
    type(FakeCMBParams), intent(inout) :: CMBin!, CMBobs
    real(dl) :: Omdmh2_out, Ombh2_out, H0_out
    real(dl), optional :: z
    real(dl) :: rbao, hbao2
    real(dl) :: thisa1, thisa2, thisa3
    real(dl) :: thish

    ! but:
    CMBin%vcp%invoid = .true.
    CMBin%vcp%fromMCMC = .false.
    CMBin%vcp%effobs = .false.
    CMBin%vcp%outside = .false.
  
    if(VoidSplitISW)then
       CMBin%vcp%splitISW = .true.
       CMBin%vcp%withISW = .true.
    end if

    CMBin%vcp%zb = VP%zb

    Omdmh2_out = CMBout%omdmh2
    Ombh2_out = CMBout%ombh2
    H0_out = CMBout%H0

    if(.not.present(z))then
       
       CMBin%H0 = VP%Hlocal_paper

       CMBin%vcp%Hlocal = VP%Hlocal_paper
       !    CMBin%omdmh2 =  Omdmh2_out / voida(0.d0,VP%t0)**3*voida(VP%L,VP%t0)**3
       !    CMBin%ombh2 =  Ombh2_out / voida(0.d0,VP%t0)**3*voida(VP%L,VP%t0)**3
!       CMBin%omdmh2 =  Omdmh2_out / voida(0.d0,VP%t0)**3*voida(VP%LTB_FLRW_r,VP%t0)**3
!       CMBin%ombh2 =  Ombh2_out / voida(0.d0,VP%t0)**3*voida(VP%LTB_FLRW_r,VP%t0)**3
       
!       CMBin%omk = 1.d0 - (CMBin%omdmh2+CMBin%ombh2) / CMBin%H0**2 * 1.d4
       
       rbao = VP%r0 ! 0.d0
       
    else ! this can only be done after the void has been initialized

       rbao = nr(z)

       ! eq. II.5
       hbao2 = 8.d0 * pi * VP%Mtilde_paper**2 / 3.d0 * (1.d0 / voida(rbao,VP%t0)**3*voida(VP%LTB_FLRW_r,VP%t0)**3&
                 + 3.d0 * voidkofr(rbao) / 4.d0 / pi / voida(rbao,VP%t0)**2*voida(VP%LTB_FLRW_r,VP%t0)**2 ) + VP%Lambda/3._dl

       CMBin%H0 = Hbao2**0.5d0
! or:
       CMBin%H0 = voidH(rbao,VP%t0)

       CMBin%vcp%Hlocal = voidH(VP%r0,VP%t0) !Hbao2**0.5d0 


    end if

! Headache:
!!$    CMBin%omdmh2 =  Omdmh2_out / voida(rbao,VP%t0)**3*voida(VP%LTB_FLRW_r,VP%t0)**3 * CMBin%H0**2 / H0_out**2
!!$    CMBin%ombh2 =  Ombh2_out / voida(rbao,VP%t0)**3*voida(VP%LTB_FLRW_r,VP%t0)**3 * CMBin%H0**2 / H0_out**2

    CMBin%omdmh2 =  Omdmh2_out / voida(rbao,VP%t0)**3*voida(VP%LTB_FLRW_r,VP%t0)**3 
    CMBin%ombh2 =  Ombh2_out / voida(rbao,VP%t0)**3*voida(VP%LTB_FLRW_r,VP%t0)**3 
    

!!$thish=voidH(r=0._dl,t=VP%t0)
!!$    write(*,'(A,ES14.5,A,ES14.5)')'H0,in Compare:',VP%Hlocal_paper,' to:', thish
!!$thish=voidH(r=VP%LTB_FLRW_r,t=VP%t0)
!!$    write(*,'(A,ES14.5,A,ES14.5)')'H0,out Compare:',H0_out,' to:', thish 
!!$writE(*,'(A,ES14.5,A,ES14.5)')'OM,out Compare:',(CMBout%omdmh2 + CMBout%ombh2) / CMBout%H0**2 * 1.d4,' to:',8.d0 * pi *VP%Mtilde_paper**2 / 3.d0 / CMBout%H0**2
!!$writE(*,'(A,ES14.5,A,ES14.5)')'OK,out Compare:',CMBout%omk,' to:',voidkofr(VP%LTB_FLRW_r)*2._dl*VP%Mtilde_paper**2 /  CMBout%H0**2
!    CMBin%omk = 1.d0 - (CMBin%omdmh2+CMBin%ombh2) / CMBin%H0**2 * 1.d4

    CMBin%omk = voidkofr(rbao)*2._dl*VP%Mtilde_paper**2 /  CMBin%H0**2*(voida(VP%LTB_FLRW_r,VP%t0)/voida(rbao,VP%t0))**2
    
    CMBin%vcp%TVoverTF = 1.d0

    CMBin%vcp%t0 = VP%t0bg * 977.8d0

!write(*,*)'test this line (makein):'
!CMBin%omv = 1._dl - CMBin%omk - ((CMBin%omdmh2 + CMBin%ombh2) / CMBin%H0**2 * 1.d4)
CMBin%omv = VP%Lambda/3._dl / CMBin%H0**2 ! this line is safer against imprecise subtraction:
                                          ! O(1) - O(1) = O(1e-8) \neq 0, for single precision
                                          ! while Lambda == omv / h0^2 * 3, which stays zero always
                                          ! once omv=0 anywhere.
                                          ! This is relevant for omv=0 cases, obviously.
!writE(*,'(A,ES14.5,A,ES14.5)')'OL Compare:',CMBin%omv,' to:',VP%Lambda/3._dl / CMBin%H0**2
!writE(*,'(A,ES14.5,A,ES14.5)')'OK Compare:',CMBin%omk,' to:',voidkofr(rbao)*2._dl*VP%Mtilde_paper**2 /  CMBin%H0**2*(voida(VP%LTB_FLRW_r,VP%t0)/voida(rbao,VP%t0))**2
!writE(*,'(A,ES14.5,A,ES14.5)')'OM Compare:',(CMBin%omdmh2 + CMBin%ombh2) / CMBin%H0**2 * 1.d4,' to:',8.d0 * pi *VP%Mtilde_paper**2 / 3.d0 / CMBin%H0**2*(voida(VP%LTB_FLRW_r,VP%t0)/voida(rbao,VP%t0))**3
    if( voidtesting .or. VoidFeedback>2)then
       thisa1=voida(0.d0,VP%t0)
       thisa2=voida(VP%LTB_FLRW_r,VP%t0)
       thisa3=voida(rbao,VP%t0)
        write(*,*)'rbao/L:',rbao/VP%L
       if(VoidFeedback>2)then ! Some definitions are wrong, some are right.
          write(*,*)'H_out^2 / H_in^2:',(VP%H0/VP%Hlocal_paper)**2
          write(*,*)'Om,in / Om,out:',(1.d0 + 3.d0*VP%kb_om / 4.d0 / pi) * (1.d0 - VP%Mtilde_paper**2  * (VP%kmax + VP%kb_om)*2.d0/ thisa1**2*thisa2**2/VP%Hlocal_paper**2)
          write(*,*)'Om,in / Om,out:',(1.d0 + 3.d0*VP%kb_om / 4.d0 / pi)/(1.d0 + 3.d0 * (VP%kmax + VP%kb_om)/4.d0 / pi)
          write(*,*)'1+Delta_0:',1.d0+VP%delta0
          write(*,*)'Omega_m,inside:',8.d0 * pi *VP%Mtilde_paper**2 / 3.d0 / VP%Hlocal_paper**2
          write(*,*)'Omega_m,outside:',1.d0 - VP%Omk
          write(*,*)'Omega_m,inside /Omega_m,outside :',8.d0 * pi *VP%Mtilde_paper**2/ 3.d0 / VP%Hlocal_paper**2/(1.d0 - VP%Omk)
          write(*,*)'Omega_k,in:',VP%Mtilde_paper**2  * (VP%kmax + VP%kb_om)*2.d0/ thisa1**2*thisa2**2/VP%Hlocal_paper**2
          write(*,*)'Omega_m,in:',8.d0 * pi *VP%Mtilde_paper**2/ 3.d0 / VP%Hlocal_paper**2/ thisa1**3*thisa2**3
          write(*,*)'Omega_k,in+Omega_m,in:',VP%Mtilde_paper**2  * (VP%kmax + VP%kb_om)*2.d0/ thisa1**2*thisa2**2/VP%Hlocal_paper**2 + 8.d0 * pi *VP%Mtilde_paper**2/ 3.d0 / VP%Hlocal_paper**2/ thisa1**3*thisa2**3
          write(*,*)'a(L,t0):',thisa2
          write(*,*)'Alessio''s omk_in:',3.d0 * (VP%kmax + VP%kb_om) / 4.d0 / pi * H0_out**2 / CMBin%vcp%Hlocal**2 / thisa1**3*thisa2**3 
          write(*,*)'My omk_in:',3.d0 * (VP%kmax + VP%kb_om) / 4.d0 / pi *(1-CMBin%omk) * thisa1/thisa2 
          write(*,*)'eq III.17 with k_tot:',3.d0 * (VP%kmax + VP%kb_om) / 4.d0 / pi  / (1 + 3.d0 * (VP%kmax + VP%kb_om) / 4.d0 / pi ) 
       end if
       write(*,*)'omk_in for cosmomc:',CMBin%omk
       write(*,*)'omv_in for cosmomc:',CMBin%omv
       write(*,*)'omm_in for cosmomc:',(CMBin%omdmh2 + CMBin%ombh2) / CMBin%H0**2 * 1.d4
       write(*,*)'compare omm_in:', 8.*pi/3._dl * VP%Mtilde_paper**2 / CMBin%H0**2 / thisa3**3*thisa2**3
       write(*,*)'omk_in + omm_in:',CMBin%omk+(CMBin%omdmh2 + CMBin%ombh2) / CMBin%H0**2 * 1.d4
       write(*,*)'omk_in + omm_in + omv_in:',CMBin%omk+(CMBin%omdmh2 + CMBin%ombh2) / CMBin%H0**2 * 1.d4 + CMBin%omv
       write(*,*)'compare omk_in + omm_in + omv_in:',CMBin%omk+ 8.*pi/3._dl * VP%Mtilde_paper**2 / CMBin%H0**2 / thisa3**3*thisa2**3+ CMBin%omv
       write(*,*)'omdmh2_in for cosmomc:',CMBin%omdmh2
       write(*,*)'ombh2_in for cosmomc:',CMBin%ombh2

       write(*,*)'H_in for cosmomc:',CMBin%H0, CMBin%vcp%Hlocal
write(*,*)'VP%LTB_FLRW_r/VP%L:',VP%LTB_FLRW_r/VP%L
       !       stop
    end if

    
    
  end subroutine VoidMakeCMBin



  subroutine VoidTestCMB(CMBin)
    implicit none
    type(FakeCMBParams), intent(in) :: CMBin
!    type(FakeCMBParams) :: myCMB
!    type(VoidParameters) :: VPmem
!    real(dl) :: LTB_maxr, LTB_FLRW_r
    logical :: thiserror
    real(dl), parameter :: cambminq = 1.d-6, cambmaxq = 0.3d0,  cambc = 2.99792458e8_dl,  cambMpc = 3.085678e22_dl
    real(dl) :: cambcurv, cambr, cambtau0, cambchi0
    logical :: cambisclosed
    real(dl), external :: rombint
    real(dl) :: argument
    integer :: cambllmax

    thiserror = .false.

    if(voidfeedback>0)write(*,*)'Testing inner FLRW.'

!!           CP%curv=-CP%omegak/((c/1000)/CP%h0)**2
!!           CP%Ksign =sign(1._dl,CP%curv)
!!           CP%r=1._dl/sqrt(abs(CP%curv))
    cambcurv=-CMBin%omk/((cambc/1000)/CMBin%h0)**2
    cambr=1._dl/sqrt(abs(cambcurv))
    cambisclosed = CMBin%omk < 0.d0
    if(cambisclosed)return ! no problems expected

    cambtau0 = voidrombint(0.d0, 1.d0, 1.d-1)

!    write(*,*)cambr, cambtau0,'------'
    !stop

    cambchi0 = dsinh(cambtau0/cambr)
    cambllmax = nint(cambmaxq * cambr*cambchi0)
!write(*,*)cambllmax, cambchi0, cambtau0/cambr
!write(*,*)nint(2.d9)
!stop

!    argument = cambmaxq * cambr *dsinh(cambtau0/cambr + 6*pi / (cambminq * cambr))
    argument = (cambtau0/cambr + 6*pi / (cambminq * cambr))
!write(*,*)argument, dasinh(max_double/cambmaxq / cambr)
!stop'.'
!    if(argument > Max_Double)thiserror = .true.
    

!    thiserror  = cambr < cambtau0
    if(.not.VoidError)VoidError = thiserror
    
    if(voidfeedback>0 .and. voiderror)write(*,*)'Inner FLRW failed. Rejecting.'


  contains
    
    function f(a)
      implicit none
      real(dl) f
      real(dl), intent(in) :: a

      f = cambc / 1.d3 / CMBin%h0 / dsqrt( CMBin%omk * a**2 + (1.d0 - CMBin%omk)*a  +1.d-10)
      
    end function f
    
	include 'vd2010/voidrombint.f90'
	
  end subroutine VoidTestCMB
    

  function VoidH0obs(z)
    implicit none
    real(dl), optional, intent(in) :: z
    real(dl) :: VoidH0obs
    real(dl) :: myz, thisda, thisv, x

    if(voiderror)then
      VoidH0obs = 1.d30
      return
    end if

    if(present(z))then
       myz=z
    else
       myz=0.04
    end if
    
!    if(VP%is_pure_flrw)then
!       VoidH0obs = VP%H0
!    else
       
       thisda = VoidDA(myz)
       x = (1._dl+myz)
       thisv = (x**2 - 1._dl) / (x**2 + 1._dl) * VPc/1.d3
       
       VoidH0obs = thisv / thisda

!        write(*,'(2(A,E14.5))')'Hobs:',VoidH0obs, 'VP%H0:',VP%H0
!    end if
    
    
    
  end function VoidH0obs

  
!   subroutine VoidCalcDerivedPars(zb,Lltb,Lga,Tspot,spotsize,H0obs,Mga,deltam,deltaw,dens_nsigma_lnlike,r3p346,Copernican,defit_w,defit_wa,defit_omk,defit_omm,winvt)
!     implicit none
!     real(dl), intent(out) :: zb,Lltb,Lga,Tspot,spotsize,H0obs,Mga,deltam,deltaw,dens_nsigma_lnlike,r3p346,Copernican,defit_w,defit_wa,defit_omk,defit_omm,winvt(:)
!     real(dl) :: thisda,ddw,ddm
!     real(dl), parameter :: thisz = 1100.d0, third = 1.d0/3.d0
  
!       zb = VP%zb
!       Lltb = VoidComovDistMPC(0.d0,VP%L,VP%t0,VP%r0)
!       Lga = VoidComovDistMPC(0.d0,VoidGARadius(),VP%t0,VP%r0)
! !      Tspot = VoidColdNewtGauge(VP%t1100,r=0.d0) * third
!       Tspot = VoidColdGeodesic(VP%t1100,r=0.d0)
!       thisda = VoidDA(thisz)
!       spotsize = 2*Lltb / (thisda*(1+thisz)) / pi * 180.d0
!       H0obs = VoidH0obs(z=0.04_dl)
!       Mga = VoidLogGAMassM0()
!       if(isnan(Mga))Mga=0.d0
!       call lltb_functions(h0_inf=vp%h0, Lambda=VP%Lambda,w=VP%w,kofr=voidkofr,dkdr=voiddkdr,tbbofr=voidtbbofr, &
!         dtbbdr=voiddtbbdr,mtilde=VP%Mtilde,r=0.d0,t=VP%t0,deltarhow=ddw,deltam=ddm)
!       deltaw = ddw
!       deltam = ddm
!       dens_nsigma_lnlike = VP%dens_nsigma_lnlike
!       r3p346 = VP%r3p346
!       Copernican = VP%Copernican
!       defit_w = VP%defit_w
!       defit_wa = VP%defit_wa
!       defit_omk = VP%defit_omk
!       defit_omm = VP%defit_omm
!       call VoidWofZDispersion(winvt)
      
!   end subroutine VoidCalcDerivedPars


!   subroutine VoidCalcDerivedParsFLRW(zb,Lltb,Lga,Tspot,spotsize,H0obs,Mga,deltam,deltaw,dens_nsigma_lnlike,r3p346,Copernican,defit_w,defit_wa,defit_omk,defit_omm,winvt)
!     implicit none
!     real(dl), intent(out) :: zb,Lltb,Lga,Tspot,spotsize,H0obs,Mga,deltam,deltaw,dens_nsigma_lnlike,r3p346,Copernican,defit_w,defit_wa,defit_omk,defit_omm,winvt(:)
!     real(dl) :: thisda,ddw,ddm
!     real(dl), parameter :: thisz = 1100.d0, third = 1.d0/3.d0
  
!       zb = 0.d0
!       Lltb = 0.d0
!       Lga = 0.d0
! !      Tspot = VoidColdNewtGauge(VP%t1100,r=0.d0) * third
!       Tspot = 0.d0
!       spotsize = 0.d0
!       H0obs = VoidH0obs(z=0.04_dl)
!       Mga = 0.d0
!       if(isnan(Mga))Mga=0.d0

!       deltaw = 0.d0
!       deltam = 0.d0
!       dens_nsigma_lnlike = 0.d0
!       r3p346 = 0.d0
!       Copernican = 0.d0
!       defit_w=0.d0
!       defit_wa=0.d0
!       defit_omk=0.d0
!       defit_omm=0.d0
!       winvt = VP%w
      
!   end subroutine VoidCalcDerivedParsFLRW
