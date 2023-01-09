
!--------------------------------------------------------------------------!
!
!  Module VoidDistances, calculates the angular diameter distance
!  for a given void profile (refs...), to be used in combination
!  with VoidDistSNChi2.
!
!  Feedback is set by the standard cosmomc feedback parameter:
!  feedback = 2  -  no output from this module.
!  feedback = 3  -  quite some output
!             4  -  a lot of output
!             5  -  an enormous amount of completely irrelevant output.
!
!  There are three relevant public subroutines:
!
!  - InitVoid(CMB) which should be called as soon as the CMB parameters
!                  are known in cosmomc. CMB is the structure of 
!                  type(FakeCMBParams), as it appears throughout cosmomc.
!                  From this subroutine the actual magic happens, and 
!                  all relevant information is stored in arrays.
!                  CMB%VCP%Hlocal is set in this subroutine.
!
!  - VoidMakeCMBin(CMBout,CMBin,z(optional)) which returns CMBin as if it is inside the void.
!                  CMBout and CMBin are structures of type(FakeCMBParams), 
!                  as it appears throughout cosmomc.
!                  CMBout is intent(in) and should contain the FLRW parameters.
!                  CMBin is intent(out) and contains the 'fake' inside
!                  FLRW parameters. z is used if you want CMB at a certain z(bao).
!
!  - VoidMakeCMBin(CMBout,CMBin,CMBobs) which returns CMBin as if it is inside the void.
!                  and CMBobs for the FLRW observer that sees the CMB identically.
!                  CMBout and CMBin are structures of type(FakeCMBParams), 
!                  as it appears throughout cosmomc.
!                  CMBout is intent(in) and should contain the FLRW parameters.
!                  CMBin is intent(out) and contains the 'fake' inside
!                  FLRW parameters. z is used if you want CMB at a certain z(bao).
!
!  - VoidTVoverTF(), a double precision function which needs no arguments
!                    it can be called after InitVoid has been called, 
!                    and returns the CMB temperature inside the void
!                    today, divided by the CMB temperature outside 
!                    the void today -> T_CMB,void / T_CMB,FLRW.
!
!  - VoidHlocal(), a double precision function which needs no arguments
!                    it can be called after InitVoid has been called, 
!                    and returns the H0 inside the void
!                    today.
!
!  - VoidDzArray(vdi), with the argument vdi of type VoidDistInfo as
!                      defined in VoidSNChi2. This structure has to 
!                      be initialized with the proper SN data, and then
!                      is returned containing lists of 
!                      mu(theory) - mu(observed), and the weights of 
!                      different points, if more than one mu corresponds
!                      to one z_observed.
!
!  - VoidDA(z),   returns DA(z)
!
!
!  There are three public logical variables:
!
!  - VoidTesting: if .true., some testing routines are done, that
!                 will give some output. You should use this variable
!                 whenever you want to test and debug, and you want 
!                 to be able to switch between testing and serious
!                 output from the params.ini you run with. Make sure
!                 this variable is set properly in driver.F90, in which
!                 all parameters are read from params.ini.
!
!  - VoidTestingIntegrationOutput: if .true. and VoidTesting is .true., 
!                                  the code will output the actual 
!                                  t(r), z(r), S(r), Sdot(r) to a file.
!
!  - VoidTweakProfile: if .true., you can use for example a temporarily
!                      different profile, and use this logical parameter
!                      at each line that you rewrite. You may want to
!                      use 'call VoidBigAlert("some message here") to
!                      warn yourself about each line you change. 
!                      VoidBigAlert has a memory of messages already 
!                      displayed, and only displays messages that 
!                      haven't been shown yet.
!
!  - VoidKofRFunction: A character string which is read from params.ini
!                      which contains k(r) in Fortran format. Parameters
!                      to use in k(r) are a (for alpha), b (for beta),
!                      kmax, kb, L and r.
!
!  - VoiddKdRFunction: A character string which is read from params.ini
!                      which contains k'(r) in Fortran format. Parameters
!                      to use in k'(r) are a (for alpha), b (for beta),
!                      kmax, kb, L and r.
!
!
!
! This is version 10^12345 - kmax, zB and H_outside as true free parameters
!
! Cosmomc gives us H_outside (named H0), we return H_local and H_outside
!
!--------------------------------------------------------------------------!

! komt dit aan? en dit?

Module VoidDistances
  use precision
  use vdii_vdipatch_FakeCMBParams
!  use VoidDistSNChi2
!  use mpk, only : Do_bestfit_output
! NEw LLTB Library:
  use wlltb
  use NROdeInt, only: isnan

  implicit none

  private



  logical :: Do_bestfit_output = .false.

  ! Save theta:
  type Thetas
     real(dl) :: in,  obs, out
  end type Thetas
  type(thetas) :: VoidTh
  ! -


  integer, parameter :: feedback = 4

  logical :: VoidTesting, VoidTestingIntegrationOutput , voidwanterrors = .false., VoidTweakProfile=.false., VoidDisallowBlueshift=.false., VoidForceCollapse=.false., VoidLSSInside=.false., VoidSplitISW = .false., VoidDeltaMatterOnly=.true.
  real(dl) :: VoidRObs
  integer :: VoidZMaxType, VoidDelta0Type = 0
  integer :: voidfeedback = 0
  integer, parameter :: maxsteps = 10000
  integer, parameter :: Lmaxsteps = 500
  real(dl), parameter :: VoidDistPrec= 1.d-5!  3
  real(dl), parameter :: VoidMinPrec = 1.d-12 ! minimum for voiddistprec: boosted in case of bad VoidSetL
  real(dl), parameter :: VPc = 2.99792458d8 ! voiddistprec => nreps now.
  real(dl), parameter :: u_precision= 1.d-9 !12!7
  real(dl), parameter :: numerical_derivatives_precision=0.d0!1.d-7!30
  real(dl) :: void_zmax=50.d0
  real(dl), parameter :: void_zmax_fixed=2100.d0 ! 50.d0 ! Changed to 1100 to autimatically check that we can reach last scattering without the metric breaking up (.e.g S going imaginary, curvature to large). 
  ! changed zmax to 2100 to be tolerant for extreme z_drag from bao. If zdrag > zmax, reject model.
  real(dl), parameter :: VoidZBPrecision = 1.d-5
  real(dl), parameter :: kmax_max = 1.d3
  real(dl), parameter :: max_double = 1.d307
  real(dl), target :: VoidInfinity = 1.d30, VoidZero = 0.d0
  real(dl) :: SetLBoost
  logical, target :: VoidError, VoidIsReset=.true.
  real(dl), target :: VoidLTB_FLRW_r
  real(dl) :: VoidLSSMeanRedshift
  logical ::  VoidDensSigmaIsSet=.false.
  logical :: VoidTestGaussianField=.true.
  logical :: VoidIncludekSZ=.false.
  logical :: VoidIncludeDipoleRadius=.false.
!  logical :: VoidIncludeFLRWDEFit=.false.
  integer :: VoidIncludeFLRWDEFit=0 
  integer :: VoidWofZObsParams=0
  real(dl) :: VoidIncludeFLRWDEFitzmax=2.d0
  integer :: VoidIncludeFLRWDEFitWeighting=0
  real(dl) :: Ombin, Omcin, hh
  real(dl) :: YHe_pass, zreio_pass

  character(len=512), allocatable :: allalerts(:), oncealerts(:)
  !  character(len=512) :: VoidKofRFunction,VoiddKdRFunction, VoidPatchSize
  !  character(len=512) :: VoidModuleName, VoidProfileModule

  !  logical :: VoidSaveModule

  integer :: VoidProfile

  logical, save :: void_first_time = .true.

  type DerivsErrors
     real(dl) :: R,Rp,Rpdot,S,Sdot, dy(4)
  end type DerivsErrors
  type(DerivsErrors) :: DER

  type vpnr_type
     real(dl), pointer :: ztr(:,:), ddztr(:,:), ztrerr(:,:)
     ! new 23/02/2011
     real(dl), pointer :: rtz(:,:),ddrtz(:,:)
     logical :: monotonic_in_z
     integer :: nmaxz, nminz
     integer, pointer :: kmaxz(:),kminz(:)
  end type vpnr_type

  type VoidParameters
     real(dl) :: L
     real(dl), pointer :: LTB_maxr, LTB_FLRW_r
     real(dl) :: L_paper, kmax, zB, Hlocal, H0, Omk, rtopini, kb_om, delta0
     real(dl) :: t0, t0bg, Mtilde, Mtilde_paper, Hlocal_paper, t0_Gyr
     real(dl) :: alpha, beta
     real(dl) :: r0, tf
     real(dl) :: L2_1, d0_2, kmax_2
     logical :: ini, ini_flrw
     real(dl) :: omdmh2 ! only for error messages
     real(dl) :: deadon ! during integration: don't miss this r.
     ! New: Omega_Lambda:
     real(dl) :: omv, lambda, tturn, bgtturn, H0_paper, Lambda_paper
     logical :: is_pure_flrw
     ! New: coldspot
     real(dl) :: t1100, tdir, coldspot, coldspotmk
     ! New numerical wLLTB
     real(dl) :: tbar, w
     ! new Variance of gaussian field for matter density:
     real(dl) :: dens_nsigma_lnlike
     ! new kSZ
     real(dl) :: kSZzmax, kSZ3000
     ! new Compton y
     real(dl) :: ComptonY
     logical :: ComptonYIsSet, kSZIsSet
     ! new dipole
     real(dl) :: r3p346, Copernican
     ! new fitting flrw curve:
     real(dl) :: defit_w, defit_wa, defit_omk, defit_omm
     logical :: defit_set
      ! new w inversion
      type(FakeCMBParams) :: CMBParams
  end type VoidParameters

  type VoidLikesMem
     character(len=64) :: name
     real :: like
  end type VoidLikesMem
  
  type(VoidParameters), save, target :: VP

  type(vpnr_type), save, target :: VPNR
  type(vpnr_type), save, target :: VPNR_FLRW

  ! Physics and cosmomc:
!  public VoidDzArray, InitVoid, VoidTVOverTF, VoidError, VoidConvert, VoidBaoTofZ, VoidDistFuncs,VoidMakeCMBin, VoidHlocal, VoidDA
    public InitVoid
    public EffectiveParams, NT_FLRW

  ! Best fit output:
  public VoidReadParams !voidsnoutput, 
  public do_bestfit_output, VoidIniBestFitOutput, VoidFinishBestFitOutput, VoidBestfitOutputCMB,VoidBestfitOutputSN,VoidBestfitOutputHST

  ! Logicals from params.ini:
  public VoidTesting, VoidTestingIntegrationOutput, VoidTweakProfile, VoidDisallowBlueshift, VoidForceCollapse, VoidRObs, VoidH0obs, VoidSplitISW, VoidDelta0Type, VoidZMaxType

  ! String containing k(r):
  !  public VoidKofRFunction, VoiddKdRFunction, VoidPatchSize
  !  public VoidModuleName, VoidProfileModule, VoidSaveModule
  public VoidProfile, VoidLSSMeanRedshift

  ! Temporarily public:
  public Init_FLRW, ResetVoid

!  public VoidSaveTheta, VoidGetThetaOut, VoidLogGAMassM0, VoidComovLTBSize

  ! public VoidCalcDerivedPars, VoidCalcDerivedParsFLRW

  public VoidAgeGYR, VoidFeedBack, Voidzb, vdii_H0avv
  
  public VoidIsReset, VoidErrorReport

  public VoidRedMessage
  
!  public VoidTopHatGaussianField, VoidTestGaussianField

  ! BAO:
  public vdi_bao_Dv, vdi_bao_proptocomov3D, vdi_bao_proptocomovFLRWDrag

  ! kSZ:
  ! public VoidIncludekSZ, vdi_kSZSetRedshifts, vdi_kSZAddkSZToClScalar, vdi_kSZ3000
  ! public dtaudrz, betaz, vdi_SetComptonY, vdi_SetComptonY_ready
  public dtaudrz, betaz
  
  ! Hofz:
  public vdi_voiddtauda, vdi_NTGyr, voidHL
  
  ! Compton-Y:
!  public vdi_ComptonY

  ! Nsigma gaussian density field:
  public vdi_HomogeneityLnLike

  public isnan, rombint

  public FakeCMBParams, VPNR, vpnr_type, VoidDA, kflrw

  public voidR, avoidRp, avoidRpdot, voida, voidaout, voidH, voidHaout, avoidSout, voidu

Contains



  subroutine ResetVoid()
    use ode_path, only : nrsuccess
    logical :: reset = .true.
    type(FakeCMBParams) :: dummy
    !write(*,*)'Reset called.'
    VP%zB = 0.d0
    VP%L = 0.d0
    VP%kmax = 0.d0
    VP%t0 = 0.d0
    VP%Hlocal = -1.d0
    VP%H0 = -1.d0
    VP%is_pure_flrw=.true.
    VP%tdir = 1._dl
    VP%dens_nsigma_lnlike=0.d0
    VoidDensSigmaIsSet=.false.
    VP%ComptonYIsSet = .true.
    VP%kSZIsSet = .true.
    VP%defit_set = .false.
    VP%defit_w = -1.d0
    VP%defit_wa = 0.d0
    if(associated(VP%LTB_FLRW_r))nullify(VP%LTB_FLRW_r)
    if(associated(VP%LTB_maxr))nullify(VP%LTB_maxr)
    call InitVoid(reset=reset)
!call lltb_feedback(1)
!call lltb_precisiongoal(1.d-16)

    VoidError = .false.
    nrsuccess = .true.
    VoidIsReset = .true.
    if(voidfeedback>0)write(*,*)'Void reset done.'
  end subroutine ResetVoid


  subroutine InitVoid(CMB, reset)
    use ode_path, only : nrsuccess,OdeVoidError, odeint_feedback
!    use camb, only : FeedBack
!    use AMLUtils, only: FeedBack
    implicit none
    logical, optional :: reset
    type(FakeCMBParams), optional :: CMB
    type(FakeCMBParams) :: CMBdummy, CMBdummy2
!    real(dl) H0
    real, save :: memzB, memkmax, memH0, memomk, memomv, memw ! same kinds as CMB !!
    logical, save :: PrevVoidError
    integer :: i,j, imax,jmax, imin, jmin,ii
    real(dl) :: thisval, thisval2, thisval3, thisval4(3), thiscmbkmax!, kmaxmem
    character(len=256) :: mystring
    character(len=10000) :: ThisErrorReport
!    namelist /CAMBPNAMELIST/ CAMBPars
!    logical, save :: warned = .false.
    ! New
    integer :: thetime(3), ndig
	character(len=32) :: timeformat
    ! Get info to our memory for possible error reports:
    if(present(reset))then
       if(reset)then
          VoidFeedBack= FeedBack - 2
          memzB = 0.d0
          memkmax = 0.d0
          memH0 = 0.d0
          memomk = 0.d0
          VP%r0 = 0.d0
          PrevVoidError = .false.
          reset = .true.
          !          if(feedback>1)write(*,*)'Void reset.',VoidProfile
          return
       else
          continue
       end if
    else
       continue
       VoidIsReset = .false.
    end if

    call VoidErrorReport(CMB)
    if(voidfeedback<-2)call VoidErrorReport()
    VP%CMBParams = CMB
    
    ! New for lltb module:
    ! Retrieve the stored error report:
    call VoidErrorReport(GetMsg=ThisErrorReport)
    ! Send to lltb module as attachement for its error reports:
    call lltb_attach_to_error(ThisErrorReport)

    if(feedback>1)then
       call system_clock(count=thetime(1),count_rate=thetime(3))
    end if



    ! Just to initialize:
    VP%LTB_FLRW_r => VoidZero

    VP%r0 = 0.d0
    VP%alpha = 1.d0
    VP%beta = 1.d0
    VP%L = 1.d0
    VP%kmax = 1.d0
    VP%kb_om = 1.d0
    VP%L2_1 = 1.d0
    VP%kmax_2 = 1.d0

  VoidProfile = CMB%VCP%profile

    if(abs(CMB%VCP%zB).gt.1.d-20)then
       ! Initialize LTB_maxr pointer:
       thisval = 1.d0
       call VoidkofrContainer(k=thisval, dk = thisval2, r = thisval)
    else
       VP%LTB_maxr => VoidZero       
    end if
    VP%alpha = 0.d0
    VP%beta = 0.d0
    VP%L = 0.d0
    VP%kmax = 0.d0
    VP%kb_om = 0.d0
    VP%L2_1 = 0.d0
    VP%kmax_2 = 0.d0


    VoidFeedBack= FeedBack - 2
    if(voidfeedback>2)call lltb_feedback(voidfeedback-2)
    odeint_feedback = voidfeedback-1

    if(CMB%omk.lt.-1.d0)then
       VoidError = .true.
       if(voidfeedback.gt.1)write(*,*)'Rejecting omk < -1.'
       return
    end if

    ! ---------------------------------------------------------- !
    if(voidtesting)then
!       open(124,file='gridding_deltaT.txt')
       write(*,'(A,7E16.7)')'# CMB%VCP%delta0,H0,zB, omk, CMB%VCP%L2_1, CMB%VCP%d0_2(cosmomc),CMB%omv:',CMB%VCP%delta0,CMB%H0,CMB%VCP%zB, CMB%omk, CMB%VCP%L2_1, CMB%VCP%d0_2,CMB%omv
!!$   imax=30
!!$   jmax = 50
!!$   imin = 1
!!$   jmin = 1
       imin = 1
       imax = 1
       jmin = 1
       jmax = 1
    else
       imax = 1
       imin = 1
       jmax = 1
       jmin = 1
    end if

    do i=imin,imax
       do j=jmin,jmax
          if(voidtesting)then
             !   CMB%VCP%zB=dble(i)/50.d0
             ! CMB%VCP%zB=dble(i-imin)/dble(imax-imin)*0.3d0 + 3.d-3
             !   CMB%VCP%delta0 = 10.d0**(dble(j)/12.5d0)
             !   CMB%VCP%delta0=-dble(j)/50.d0
             !   CMB%VCP%delta0=dble(j)/100.d0
             ! CMB%VCP%delta0 = dble(j-jmin)/dble(jmax-jmin)*50.d0!dble(j)!-dble(j)/26.d0
             !   CMB%omk = dble(i-imin)/dble(imax-imin)-0.5d0

             !   voidfeedback=0!1!0!4!0
             voidwanterrors =.false. ! .true.
             voiderror = .false.
!             CMB%VCP%delta0=0.d0
!             CMB%VCP%zB = 0.d0
!             CMB%omk = .5d0 -dble(j-1)/dble(jmax-1)
          end if
          ! ---------------------------------------------------------- !




          VP%alpha = CMB%VCP%alpha 
          VP%beta = CMB%VCP%beta



          if(CMB%VCP%zB==memzB .and. CMB%VCP%delta0==memkmax  .and. CMB%H0==memH0 .and. CMB%omk==memomk .and. CMB%omv==memomv .and. CMB%w==memw)then
             if(voidfeedback>0)write(*,*)'Already initialized everything for this cosmology.'!zBgoal,kmax,thisr
             if(PrevVoidError)then
                VoidError=PrevVoidError
             end if
             if(.not.voidtesting)return
          end if
          memzB = CMB%VCP%zB
          memkmax = CMB%VCP%delta0
          memH0 = CMB%H0
          memomk = CMB%omk
          memomv=CMB%omv
          memw=CMB%w
          PrevVoidError = .false.

          OdeVoidError => VoidError 


          VP%ini=.true.
          nrsuccess=.true.
          VoidError = .false.
          VP%is_pure_flrw=.true.
          VP%tdir = 1._dl

          ! For debugging purposes / error messages:
          VP%omdmh2 = CMB%omdmh2


          ! Just to initialize:
          VP%LTB_FLRW_r => VoidZero


          SetLBoost = 1.d0
          void_zmax = void_zmax_fixed

          if(abs(CMB%VCP%zB).gt.1.d-30)then
             ! WV new collapses:
             VP%is_pure_flrw=.false.

             ! end wv new

             VP%H0_paper = 0._dl + dble(CMB%H0)

             VP%Hlocal = 0.d0!CMB%VCP%Hlocal
             VP%zB = CMB%VCP%zB
             
             VP%omk = CMB%omk ! Same sign!! omk = 1 - om_dm - om_DE
             if(abs(VP%omk)<1.d-7)VP%omk=0.d0
             ! Trick to force use of lambda even in flat universE:
             if(dabs(vp%omk)<1.e-15_dl)VP%omk=0._dl
             VP%omv = CMB%omv
             VP%w = CMB%w

             VP%kb_om = 4.d0 * pi/3.d0 * ( VP%omk / (1.d0 - VP%omk - VP%omv) )
             
             VP%lambda_paper = 3._dl * CMB%H0**2 * CMB%omv

             ! New: delta_0 or kmax:
             ! VP%kmax = CMB%VCP%delta0
             ! VP%kmax =
             thiscmbkmax = CMB%VCP%delta0
!!$             call VoidKmaxofCMBkmax(thiscmbkmax)
!!$             call sett0bg_Mt

!!$             if(VoidProfile==3)then
!!$                kmaxmem = VP%kmax
!!$                VP%L2_1 = CMB%VCP%L2_1
!!$                VP%d0_2 = CMB%VCP%d0_2
!!$                call VoidKmaxofCMBkmax(VP%d0_2)
!!$                VP%kmax_2 = VP%kmax
!!$                VP%kmax = kmaxmem
!!$                if(voidfeedback>2)then
!!$                   write(*,*)'kmax1:', VP%kmax
!!$                   write(*,*)'kmax2:', VP%kmax_2
!!$                   write(*,*)'d01:',CMB%VCP%delta0
!!$                   write(*,*)'d02:',CMB%VCP%d0_2
!!$                   write(*,*)'zB1:',VP%zB
!!$                   write(*,*)'L2/L1:',VP%L2_1
!!$                end if
!!$             end if
             ! if(voidfeedback.gt.1)then I don't needed
             if(voidfeedback.gt.3)then
                write(*,'(A,E17.8)')'CMB omegabh2',CMB%ombh2
                write(*,'(A,E17.8)')'CMB omegadmh2',CMB%omdmh2
                write(*,'(A,E17.8)')'CMB omegaH0',CMB%H0
                write(*,'(A,E17.8)')'CMB omegak',CMB%omk
                write(*,'(A,E17.8)')'CMB zb',CMB%VCP%zb
                write(*,'(A,E17.8)')'CMB delta0',CMB%VCP%delta0
                write(*,'(A,E17.8)')'CMB omegav',CMB%omv

             end if

             if(voidfeedback.gt.0)write(*,*)'Initializing void for this cosmology.'

             call VoidKmaxofCMBkmax(CMB%VCP)
             if(.not.(voiderror) .or. voidtesting)call sett0bg_Mt

             if(voidfeedback.gt.1)then
                write(*,'(A,6E16.7)')' # CMB%VCP%delta0,H0,zB, omk(cosmomc), omv(cosmomc):',CMB%VCP%delta0,CMB%H0,CMB%VCP%zB, CMB%omk, CMB%omv
                write(*,'(A,8E17.8)')' # VP%kmax, zB, omk, H0_out, omdmh2, omv: ',VP%kmax,VP%zB,VP%omk,VP%H0_paper, CMB%omdmh2,VP%omv
!!$                write(*,'(A,E17.8)')'CMB omegadmh2',CMB%omdmh2
!!$                write(*,'(A,E17.8)')'CMB omegabh2',CMB%ombh2
!!$                write(*,'(A,E17.8)')'CMB omegaH0',CMB%H0
!!$                write(*,'(A,E17.8)')'CMB omegak',CMB%omk
!!$                write(*,'(A,E17.8)')'CMB zb',CMB%VCP%zb
!!$                write(*,'(A,E17.8)')'CMB delta0',CMB%VCP%delta0
             end if

             VP%L = 0.d0

!              if(.not.voiderror .or. (voidtesting))then
! !                call init_FLRW
! !                vp%t1100=nt_FLRW(z=1100._dl)
! !                write(*,*)'t(z=1100):',vp%t1100
!              end if


             if( .not. nrsuccess ) VoidError = .true.

!!$             if(VoidTesting)then
!!$                VP%L=3.d-3
!!$                VP%kmax=-10.d0
!!$                call PlotDensity()
!!$                stop
!!$             end if
 !           
             !       call voidtestanafunc() ! if you like..
             if(voiderror)call Dothemagic() ! JUST TO GET FIRST ALLOCATIONS
             ! if(.not.voiderror)then! Background cosmology is ill, e.g. too large curvature or so      
             if(.not.voiderror)then! Typical computation, not bugged cosmology      
                call VoidSetL ! the main reason to be slow
                if(feedback>2)then
                   call system_clock(count=thetime(2))
                   ndig = int(log10(dble(thetime(2)-thetime(1))/dble(thetime(3))))
                    if(ndig>0 .and. ndig<30)then
                       write(timeformat,*)ndig+6
                       timeformat='(A,F'//trim(adjustl(timeformat))//'.3,A)' !'(A,F6.3,A)'
                       write(*,timeformat)' Done VoidSetL in',dble(thetime(2)-thetime(1))/dble(thetime(3)),' seconds.'
                    end if
                end if
              if(voiderror)call Dothemagic() ! JUST TO GET FIRST ALLOCATIONS
             end if
             !if(VoidTesting .and. VP%L>0.d0)stop ' check output'

!                call init_FLRW
!                vp%t1100=nt_FLRW(z=1100._dl)
                vp%t1100=nt(z=1099.9_dl)
                if(vp%t1100<=0.d0 .and. .not. voiderror)then ! this bug is solved.
                  call VoidErrorReport()
                  write(*,*)'t1100:',VP%t1100
                  write(0,*)'t1100:',VP%t1100
                  stop
                  
                end if
             if(.not.voiderror)call VoidSetAsymp()
             if(.not.voiderror)call VoidSetkSZZmax() ! in voidasymp.f90


!!$             if(voiderror .and. voidtesting)then
!!$                write(*,*)'Done with setasymp.'
!!$                VoidTestingIntegrationOutput = .true.
!!$                vp%t0 = nt_flrw(r=0.d0)
!!$                voiderror = .false.
!!$                call init_flrw
!!$                write(*,*)nr_FLRW(t=nt(VP%LTB_FLRW_r))-VP%LTB_FLRW_r,nr_FLRW(t=nt(VP%LTB_FLRW_r)),VP%LTB_FLRW_r
!!$                write(*,*)voiderror
!!$                stop 'Done and crashed.'
!!$             end if


             CMB%VCP%Hlocal = VP%Hlocal_paper
             
             ! New gaussian density field:
            ! call VoidTopHatGaussianField(CAMBPars,VP%dens_nsigma_lnlike)

             VP%r3p346 = 0.d0
             VP%Copernican = 1.d0
!             if(VoidTesting)call VoidAsympDipole_vs_r()
             
             !write(0,*) VoidIncludeFLRWDEFit
             !if(VoidIncludeDipoleRadius.and..not.voiderror)call VoidAsympDipole()
              
             !if(VoidIncludeFLRWDEFit>0.and..not.voiderror) call vdi_GetIgnObsDEPars()
  
             if(voidfeedback.gt.0)write(*,*)'Void initialized.'

!              if(Do_bestfit_output)then
! !write(*,*)'Skipping plotdensity. This is a temporary hack.'
!                 if(VoidError)stop 'Error while plotting a best fit model?'
!                 call PlotDensity()
! !                iF(VP%delta0<0.)then
!                 !if(VoidRObs==0.d0)
!                 ! call VoidAsympDipole()
!                    ! crunches are rejected anyway. Should not crunch.
! !                else
! !                   write(*,*)'Skipping r_3.55mk calculation for overdensities: risk of big crunch at centre.'
! !                end iF
!              end if

          else
             if(voidfeedback.gt.1)write(*,*)'VoidDistances is initializing FLRW cosmology.'
             ! WV new collapses:
             if((.not.VoidTesting).and.(.not.Do_bestfit_output))then
                VP%is_pure_flrw=.true.
             else
                VP%is_pure_flrw=.false.
             end if
             ! end wv new

             VP%kmax = 0.d0
             VP%L = 0.d0
             VP%H0_paper = CMB%H0
             VP%Hlocal = 0.d0
             VP%zB = CMB%VCP%zB
             VP%omk = CMB%omk ! Same sign!! omk = 1 - om_dm
             VP%omv = CMB%omv
             VP%w = CMB%w
             VP%kb_om = 4.d0 * pi/3.d0 * ( VP%omk / (1.d0 - VP%omk - VP%omv) )
             VP%lambda_paper = 3._dl * CMB%H0**2 * CMB%omv
             thiscmbkmax = CMB%VCP%delta0
!             call VoidKmaxofCMBkmax(thiscmbkmax)
             call VoidKmaxofCMBkmax(CMB%VCP)
             call sett0bg_Mt
             CMB%VCP%Hlocal = CMB%H0
             call init_FLRW
          end if
          
          if(voidtesting)then ! uncomment whatever you want.
!              call vdi_GetIgnObsDEPars()
!             write(*,*)'While testing, boost now:',SetLBoost
!             SetLBoost = SetLBoost * 1.d-2
!             write(*,*)'Changed into:',SetLBoost
!             if(voiderror)stop 'You want to test, but we crashed with an error already.'
!             call testmakedz ! if you like..
!             stop
             ! call CompareWithMath()
             ! call VoidTestGA()
             ! call PlotKofr()
             ! call VoidTestVolumes()

!!$             voiderror=.false.
!!$             thisval = VP%r0
!!$             VP%r0=0.d0
!!$             void_zmax=void_zmax_fixed
!!$             call Dothemagic()
!!$             VP%r0 = thisval
            ! call testVoidTopHatGaussianField()
             ! if(present(CAMBPars))call VoidTopHatGaussianField(CAMBPars,VP%dens_nsigma_lnlike)
             ! if(present(CAMBPars))write(*,CAMBPNAMELIST)
            ! call PlotDensity()
!              write(*,*)nz(r=VP%L)
!              write(*,*)VoidDA(z=1000.d0)
!              stop 'that was z(L) and DA(1000)'
             ! call VoidMakeCritDens()
             ! call TestColdSpot() 
             ! call VoidGetColdSpot()
             ! call VoidAsymp()
             ! call VoidAsympDipole()
             ! call VoidWofZ()
             ! call testVoidDumpDAofz()
             ! call VoidAsympDipole_vs_r()
             ! call VoidErrorReport()
             ! call VoidDensAtCMB()
             !! VP%coldspotmk=VoidLogGAMassM0()
             ! call VoidSNMock()
             ! call wlltb_error('VoidTesting')
             continue
          end if

          VP%ini=.false.

          if( .not. nrsuccess ) VoidError = .true.


          if(voidfeedback>0)then
             if(VoidError.and.(dabs(VP%omk)+dabs(VP%kmax/20.d0).lt.1.d0))then ! very negative omk and large kmax both give difficulties that we are aware of.
                write(*,*)'Rejecting model due to numerical trouble.'
                writE(*,'(A,5E17.8)')' #kmax, zB, omk, H0_out, omdmh2',VP%kmax,VP%zB,VP%omk,VP%H0, CMB%omdmh2
                call EXIT(voidfeedback - voidfeedback)
               
             end if

              !if(VoidError.and.SetLBoost*VoidDistPrec.lt.VoidMinPrec)then ! is not possible reach a good precision for L
              if(VoidError.and.SetLBoost*VoidDistPrec.lt.VoidMinPrec)then ! should avoid all errors
                write(*,*)'Rejecting model due to numerical trouble.'
                writE(*,'(A,5E17.8)')' #kmax, zB, omk, H0_out, omdmh2',VP%kmax,VP%zB,VP%omk,VP%H0, CMB%omdmh2
                call EXIT(voidfeedback - voidfeedback)
               
             end if
          end if
          if(voiderror)PrevVoidError = VoidError


          !if(voidtesting)call EstimateErrors()

          ! ---------------------------------------------------------- !
          ! if(voidtesting)call InitVoidMakeOutput
       end do
       if(voidtesting)then
          write(*,*)' '
!          write(124,*)' '
       end if
    end do
    if(voidtesting)close(124)
!    if(voidtesting)call VoidConvert(CMB,CMBdummy, CMBdummy2)
    if(feedback>2)then
       call system_clock(count=thetime(2))
       ! write(*,*)'Done InitVoid in',dble(thetime(2)-thetime(1))/dble(thetime(3)),'seconds.'
 				   ndig = int(log10(dble(thetime(2)-thetime(1))/dble(thetime(3))))
				   write(timeformat,*)ndig+6
				   timeformat='(A,F'//trim(adjustl(timeformat))//'.3,A)' !'(A,F6.3,A)'
       write(*,timeformat)' Done VoidDistances in ',dble(thetime(2)-thetime(1))/dble(thetime(3)),' seconds.'
    end if
      !if(voidtesting)stop 'Done testing.'
      ! moved the stop call to CAMB_TransfersToPowers in vdii_cambcamb.f90'
    ! ---------------------------------------------------------- !

    ! call InitVoidTestzFLRW ! if you like..

    return
  contains

!     subroutine InitVoidMakeOutput()
!       implicit none
!       real(dl) :: kg2t02, krb2

!       thisval3 = VPNR%ztrerr(1, GetNearestIndex(VPNR%ztr(1,:),VP%zB)     )

!       thisval4 = VPNR_FLRW%ztrerr(1, GetNearestIndex(VPNR_FLRW%ztr(1,:),z_FLRW(VP%L))   )

!       !(1.d0 + z_FLRW(VP%L)) / (1.d0 + VP%zB)

!       thisval=VoidTVoverTF()

!       thisval2=AlessiosEstimate()

!       thisval3=dsqrt( ( thisval4(1)/(1.d0+VP%zB) )**2 + (  (1+z_FLRW(VP%L))/(1+VP%zB)**2 * thisval3  )**2  ) 

!       if(VP%kb_om == 0.d0) then
!          thisval4(1) = 1.d0/void_B21()
         
!          thisval4(2) = 1.d0/void_B12()
!       else
!          thisval4 = 0.d0
!       end if

!       call TirthosTinyNumbers(kg2t02,krb2)

!       ! With error estimates:
!       !       write(mystring,'(9F16.7, E14.5)')CMB%VCP%zB,VP%L*VP%Mtilde,CMB%VCP%delta0, VP%delta0,1.d0-thisval, thisval2,1.d0-thisval4(1),  (1.d0 - thisval)/thisval2 - 1.d0,  (1.d0-thisval)/(1.d0-thisval4(1))  - 1.d0, thisval3
!       ! With expansions parameters:
!       !       write(mystring,'(9F16.7, E14.5)')CMB%VCP%zB,VP%L*VP%Mtilde,CMB%VCP%delta0, VP%delta0,1.d0-thisval, thisval2,1.d0-thisval4(1),  kg2t02,krb2
!       !       write(mystring,'(12ES16.7)')CMB%VCP%zB,VP%L_paper*VP%Mtilde_paper,VP%kmax, VP%delta0,1.d0-thisval, thisval2,1.d0-thisval4(1),  kg2t02,krb2,CMB%omk
!       write(mystring,'(12ES16.7)')CMB%VCP%zB,VP%L*VP%Mtilde,VP%kmax, VP%delta0,1.d0-thisval, thisval2,1.d0-thisval4(1),1.d0-thisval4(2),  kg2t02,krb2,VP%kmax*(VP%L*VP%Mtilde)**2
! !      write(mystring,'(12ES16.7)')CMB%omk, DA(1000.d0)

! !      write(124,'(A)')trim(mystring)
!       write(*,'(A)')' '//trim(mystring)

!     end subroutine InitVoidMakeOutput


    ! subroutine InitVoidTestzFLRW()


    !   do i = 1, 25
    !      thisval = VP%L*dble(i)/25.d0
    !      thisval3 = z_FLRW(thisval)
    !      thisval2 = simpler(thisval3) * VP%t0/VP%t0bg 
    !      write(*,'(5F16.9)')thisval,thisval2, thisval/thisval2, thisval3, z_FLRW(thisval2)
    !   end do
    !   stop

    ! end subroutine InitVoidTestzFLRW


  end subroutine InitVoid

#include "vd2010/rombint.f90"

#include "vd2010/void_bao_dist.f90"

#include "vd2010/voidanafuncs.f90"

#include "vd2010/voidnumfuncs.f90"

#include "vd2010/voidasymp.f90"

! #include "vd2010/voidcoldspot.f90"

#include "vd2010/voidGA.f90"

! #include "vd2010/voidvolumes.f90"

!  include "void/voiderror.f90"

! #include "vd2010/voidestimates.f90"

#include "vd2010/voidintegrator_dvode.f90"

#include "vd2010/voidkofr.f90"

! #include "vd2010/voidminimizechi2.f90"

#include "vd2010/voidtbbofr.f90"

#include "vd2010/voidnsolve.f90"

#include "vd2010/voidprefuncs.f90"

#include "vd2010/voidprecomm.f90"

#include "vd2010/voidpostcomm.f90"

#include "vd2010/voidpostfuncs.f90"

! #include "vd2010/voidsnoutput.f90"

#include "vd2010/voidtestrout.f90"

#include "vd2010/voidutils.f90"

#include "vd2010/voidoutput.f90"

#include "vd2010/voidksz.f90"

! #include "vd2010/voidminchi2.f90"

! #include "vd2010/voidminimc.f90"

end Module VoidDistances
