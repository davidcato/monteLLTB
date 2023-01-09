program VoidDistances2020
    use VoidDistances
    use precision
    type(FakeCMBParams) CMB
    integer :: i
    character :: action
    logical :: do_basic, do_complete
    real(dl) :: zdrag_in, tdrag_BAO, YHe_in, zreio_in
    write(*,*) "Hello void world!"

    write(*,*) "Usage:"
    write(*,*) trim(adjustl(GetCommandName())), " action H0 omega_b_h2 omega_dm_h2 Omega_k Omega_DE w_DE delta0 z_Boundary z_drag"
    write(*,*) ""
    write(*,*) "where action is either:"
    write(*,*) "  b  (basic/fastest computation, not suitable for Y Compton & kSZ)"
    write(*,*) "  c  (complete/slowest computation, used to Y Compton & kSZ) you should also pass YHe z_reio"
    write(*,*) ""
    write(*,*) "and all other parameters are real values."
    write(*,*)

    write(*,*) "Will crash if you do not supply these."

    CMB%VCP%profile = 19
    CMB%VCP%alpha = 1.0
    CMB%VCP%beta = 1.0

    i = 0

    i = i + 1
    call get_command_argument(i, action)

    do_basic = action == "b"
    do_complete = action == "c"

    i = i + 1
    CMB%H0 = GetArgFloat(i)
    i = i + 1
    CMB%omegabh2 = GetArgFloat(i)
    i = i + 1
    CMB%omegadmh2 = GetArgFloat(i)
    i = i + 1
    CMB%omk = GetArgFloat(i)
    i = i + 1
    CMB%omv = GetArgFloat(i)
    i = i + 1
    CMB%w = GetArgFloat(i)
    i = i + 1
    CMB%VCP%delta0 = GetArgFloat(i)
    i = i + 1
    CMB%VCP%zB = GetArgFloat(i)    
    i = i + 1
    zdrag_in = GetArgFloat(i)
    if ( do_complete ) then
        i = i + 1
        YHe_in = GetArgFloat(i)
        i = i + 1
        zreio_in = GetArgFloat(i)
    else
        YHe_in = 0.d0
        zreio_in = 0.d0
    end if

    CMB%omdmh2 = -10000 ! unused.

    CMB%YHe = YHe_in
    CMB%zreio = zreio_in

    call ResetVoid()

    VoidTesting = .true.
    VoidTestingIntegrationOutput = .false.

    call InitVoid(CMB)

    call EffectiveParams(CMB)

    if ( do_basic ) then
        call DumpBasic(6)
    end if

    if ( do_complete ) then
        call DumpComplete(6)
    end if

contains

    function GetCommandName()
        character(len=256) GetCommandName

        call get_command_argument(0, GetCommandName)

    end function GetCommandName

    function GetArgFloat(i)
        integer, intent(in) :: i
        real(kind(1.d0)) GetArgFloat
        character(len=64) :: arg

        call get_command_argument(i, arg)

        read(arg, *) GetArgFloat

    end function GetArgFloat

    subroutine DumpBasic(outfn)

        integer :: outfn, i, imax

        ! imax = size(VPNR%ztr(1,:))
        imax = 180 ! zmax - 20. bigger than z_reio

        tdrag_BAO = NT_FLRW(zdrag_in)

        !write(outfn, '(100A15)') "z", "t", "r", "d_A", "HL",  "a", "a_out", "a0_r", "H", "H_out", "H0_r", "R", "R0_r", "Rz_tdrag", "Rp", "Rp0_r", "Rpz_tdrag", "Rpdot"
        write(outfn, '(100A15)') "z", "t", "r", "d_A", "HL",  "a", "a_out", "a0_r", "H", "H_out", "H0_r", "R", "R0_r", "Rz_tdrag", "Rp", "Rp0_r", "Rpz_tdrag", "Rpdot", "Rp0dot_r"
        do i = 1, imax

            ! write(outfn, '(100E20.11E3)') VPNR%ztr(1:3, i), VoidDA(VPNR%ztr(1, i)), voidHL(VPNR%ztr(1, i)), voida(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidaout(VPNR%ztr(2, i)), voida(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidH(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidHaout(VPNR%ztr(2, i)), voidH(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidR(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidR(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidR(VPNR%ztr(3, i),tdrag_BAO,voidu(VPNR%ztr(3, i),VPNR%ztr(2, i))), nvoidRp(VPNR%ztr(3, i),VPNR%ztr(2, i)), nvoidRp(VPNR%ztr(3, i),VPNR%ztr(2, 1)), nvoidRp(VPNR%ztr(3, i),tdrag_BAO,voidu(VPNR%ztr(3, i),VPNR%ztr(2, i))), nvoidRpdot(VPNR%ztr(3, i),VPNR%ztr(2, i))
            write(outfn, '(100E20.11E3)') VPNR%ztr(1:3, i), VoidDA(VPNR%ztr(1, i)), voidHL(VPNR%ztr(1, i)), voida(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidaout(VPNR%ztr(2, i)), voida(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidH(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidHaout(VPNR%ztr(2, i)), voidH(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidR(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidR(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidR(VPNR%ztr(3, i),tdrag_BAO), avoidRp(VPNR%ztr(3, i),VPNR%ztr(2, i)), avoidRp(VPNR%ztr(3, i),VPNR%ztr(2, 1)), avoidRp(VPNR%ztr(3, i),tdrag_BAO), avoidRpdot(VPNR%ztr(3, i),VPNR%ztr(2, i)), avoidRpdot(VPNR%ztr(3, i),VPNR%ztr(2, 1))

        end do

        write(*,*) imax

    end subroutine DumpBasic


    subroutine DumpComplete(outfn)

        integer :: outfn, i, imax
        real(dl) :: YHe_pass, zreio_pass

        YHe_pass = CMB%YHe * 1.d0
        zreio_pass = CMB%zreio * 1.d0

        ! imax = size(VPNR%ztr(1,:))
        imax = 180 ! zmax - 20. bigger than z_reio

        tdrag_BAO = NT_FLRW(zdrag_in)

        ! write(outfn, '(100A15)') "z", "t", "r", "d_A", "HL", "a", "a_out", "a0_r", "H", "H_out", "H0_r", "R", "Rz_tdrag", "Rp", "Rpz_tdrag", "Rpdot", "dtaudrz", "betaz", "k_back"
        write(outfn, '(100A15)') "z", "t", "r", "d_A", "HL", "a", "a_out", "a0_r", "H", "H_out", "H0_r", "R", "R0_r", "Rz_tdrag", "Rp", "Rp0_r", "Rpz_tdrag", "Rpdot", "Rp0dot_r", "dtaudrz", "betaz", "k_back"
        do i = 1, imax

            ! write(outfn, '(100E20.11E3)') VPNR%ztr(1:3, i), VoidDA(VPNR%ztr(1, i)), voidHL(VPNR%ztr(1, i)), voida(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidaout(VPNR%ztr(2, i)), voida(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidH(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidHaout(VPNR%ztr(2, i)), voidH(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidR(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidR(VPNR%ztr(3, i),tdrag_BAO,voidu(VPNR%ztr(3, i),VPNR%ztr(2, i))), nvoidRp(VPNR%ztr(3, i),VPNR%ztr(2, i)), nvoidRp(VPNR%ztr(3, i),tdrag_BAO,voidu(VPNR%ztr(3, i),VPNR%ztr(2, i))), nvoidRpdot(VPNR%ztr(3, i),VPNR%ztr(2, i)), dtaudrz(VPNR%ztr(1, i),YHe_pass,CMB%omegabh2,CMB%omegadmh2,CMB%H0), betaz(VPNR%ztr(1, i),zreio_pass), kflrw(VPNR%ztr(1, i))
            write(outfn, '(100E20.11E3)') VPNR%ztr(1:3, i), VoidDA(VPNR%ztr(1, i)), voidHL(VPNR%ztr(1, i)), voida(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidaout(VPNR%ztr(2, i)), voida(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidH(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidHaout(VPNR%ztr(2, i)), voidH(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidR(VPNR%ztr(3, i),VPNR%ztr(2, i)), voidR(VPNR%ztr(3, i),VPNR%ztr(2, 1)), voidR(VPNR%ztr(3, i),tdrag_BAO), avoidRp(VPNR%ztr(3, i),VPNR%ztr(2, i)), avoidRp(VPNR%ztr(3, i),VPNR%ztr(2, 1)), avoidRp(VPNR%ztr(3, i),tdrag_BAO), avoidRpdot(VPNR%ztr(3, i),VPNR%ztr(2, i)), avoidRpdot(VPNR%ztr(3, i),VPNR%ztr(2, 1)), dtaudrz(VPNR%ztr(1, i),YHe_pass,CMB%omegabh2,CMB%omegadmh2,CMB%H0), betaz(VPNR%ztr(1, i),zreio_pass), kflrw(VPNR%ztr(1, i))

        end do

        write(*,*) imax

    end subroutine DumpComplete

end program VoidDistances2020
