!  AdjointOceanModel.f90 
!
!  FUNCTIONS:
!  AdjointOceanModel - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: AdjointOceanModel
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
    program AdjointOceanModel
    use Const
    use Ocean_Model
    use Adjoint_Model
    
    ! read constants and parameters
    call Const_init
    call Ocean_Model_allocate
    call Adjoint_Model_allocate
    
    ! 如果不优化开边界、底摩擦，说明要运行正向模式
    if ( (.not. optOpenBoundary) .and. (.not. optFrictionData) ) then
        write(*, *) "Now running the 2-D forward simulation."
        call Ocean_Model_Run
        ! 调和分析
        call Harmonic
        ! 输出整个区域的调和常数
        call dump2D(simZeta, 'Zeta.dat')
        call dump2D(simSigma, 'Sigma.dat')
        ! 输出水位、流场变化
        call dump3D(H, 'H.dat')
        call dump3D(U, 'U.dat')
        call dump3D(V, 'V.dat')
    else
        write(*, *) "Now running the adjoint optimization."
        call Adjoint_Model_Run
        if ( optOpenBoundary ) then
            ! 注意这里优化出的A、B仅仅是开边界上面的
            call dump2D(Tide_A, 'A.dat')
            call dump2D(Tide_B, 'B.dat')
            ! 注意这里的振幅迟角是针对整个模拟区域的
            call dump2D(simZeta, 'Zeta.dat')
            call dump2D(simSigma, 'Sigma.dat')
            ! 注意这里的观测点与模拟值之差是只在观测点上的
            call dump2D(diffZeta, 'diffZeta.dat')
            call dump2D(diffSigma, 'diffSigma.dat')
        end if
        if ( optFrictionData ) then
            ! 注意这里输出的底摩擦则是针对整个流场的
            call dump2D(Friction, 'Friction.dat')
        end if
    endif
    
    end program AdjointOceanModel
