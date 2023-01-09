
subroutine VoidReadParams()
    write(*,*)"VoidReadParams does nothing."
    return
!  use IniFile
!  implicit none
!  logical, save :: firstcall = .true.

!  if(.not.firstcall)return
!  firstcall = .false.

!  VoidTesting = Ini_Read_Logical('VoidTesting',.false.)
!!  VoidTesting = Ini_Read_Logical('VoidTesting',.true.)
!  Do_bestfit_output = Ini_Read_Logical('Do_bestfit_output',.false.)
!  VoidTestingIntegrationOutput = Ini_Read_Logical('VoidTestingIntegrationOutput',.false.)
!!  VoidTestingIntegrationOutput = Ini_Read_Logical('VoidTestingIntegrationOutput',.true.)
!  VoidDisallowBlueshift = Ini_Read_Logical('VoidDisallowBlueshift',.false.)
!  VoidForceCollapse = Ini_Read_Logical('VoidForceCollapse',.false.)
!  VoidProfile = Ini_Read_Int('VoidProfile',19)
!  VoidLSSMeanRedshift = dble(Ini_Read_Real('VoidLSSMeanRedshift',0.))
!  VoidLSSInside = Ini_Read_Logical('VoidLSSInside',.false.)
!  VoidRObs = dble(Ini_Read_Real('VoidRObs',0.))
!  VoidSplitISW = Ini_Read_Logical('VoidSplitISW',.false.)
!
!  VoidDelta0Type = Ini_Read_Int('VoidDelta0Type',0)
!  VoidDeltaMatterOnly = Ini_Read_Logical('VoidDeltaMatterOnly',.true.)
!!  VoidZMaxType = Ini_Read_Int('VoidZMaxType',0)
!  VoidZMaxType = Ini_Read_Int('VoidZMaxType',2)
!  VoidTestGaussianField = Ini_Read_Logical('VoidTestGaussianField',.true.)
!  ! kSZ:
!  VoidIncludekSZ = Ini_Read_Logical('VoidIncludekSZ',.true.)
!  ! Dipole:
!  VoidIncludeDipoleRadius = Ini_Read_Logical('VoidIncludeDipoleRadius',.false.)
!  ! flrw w fit:
!  VoidIncludeFLRWDEFit = Ini_Read_Int('VoidIncludeFLRWDEFit',0)
!  VoidIncludeFLRWDEFitzmax = Ini_Read_Real('VoidIncludeFLRWDEFitzmax',2.0)
!  VoidIncludeFLRWDEFitWeighting = Ini_Read_Int('VoidIncludeFLRWDEFitWeighting',0)
!  VoidWofZObsParams = Ini_Read_Int('VoidWofZObsParams',0)

end subroutine VoidReadParams
