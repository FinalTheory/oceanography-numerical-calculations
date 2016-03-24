    module Const
    implicit none
    ! 定义全局变量精度
    integer, parameter              :: PREC = selected_real_kind(15)
    
    real(kind=PREC), parameter      :: PI = dacos(-1.d0), deg2rad = PI / 180.d0, &
                                        Omega = 2.d0 * PI / ( 24.d0 * 3600.d0 )
    real(kind=PREC)                 :: g, TideOmega, EarthRadius
    real(kind=PREC)                 :: Lon, Lat, Resolution
    
    logical                         :: Coriolis, optOpenBoundary, optFrictionData
    character*80                    :: case_dir, output_dir, water_depth_file, border_file, &
                                        observe_data_file, open_boundary_file, friction_file, config_name = 'config.nml'
    ! 时空步数，观测点个数，开边界点个数，最大迭代次数
    integer                         :: numX, numY, numSteps, numObserveData, numOpenBoundary, maxIteration, observeFormat
    
    ! 本模式暂时不支持优化涡动粘性系数
    real(kind=PREC)                 :: T, dt, dy, ddy, globalFriction, globalViscosity, &
                                        IntegrationEps, eps, SigmaEps, ZetaEps, alphaFriction, &
                                        alphaTideAB, alphaTideABdec, alphaTideABmin, alpha
    integer                         :: x, y, sumObserve, ExitType
    
    real(kind=PREC), allocatable, &
    dimension(:)                    :: dx, ddx, f
    
    ! 水深，观测点上的振幅、迟角（不变），底摩擦（可变）
    ! 开边界上当前的、优化后的调和常数(a, b)，优化出的振幅、迟角
    real(kind=PREC), allocatable, &
    dimension(:, :)                 :: WaterDepth, obsZeta, obsSigma, Friction, newFriction, Viscosity, &
                                        Tide_A, Tide_B, newTide_A, newTide_B, optZeta, optSigma
    
    ! 代表调和常数观测值导出的水位变化（即观测点处的水位变化）
    real(kind=PREC), allocatable, &
    dimension(:, :, :)              :: H0
    
    integer, allocatable, &
    dimension(:, :)                 :: maskU, maskV, maskH, maskBoundary, maskObserve
    
    
    namelist /CONF/                 numX, numY, numSteps, numObserveData, numOpenBoundary, maxIteration, &
                                    observe_data_file, open_boundary_file, friction_file, &
                                    optOpenBoundary, optFrictionData, observeFormat, &
                                    g, TideOmega, Lon, Lat, EarthRadius, Resolution, &
                                    water_depth_file, border_file, ExitType, &
                                    Coriolis, globalFriction, globalViscosity, &
                                    IntegrationEps, eps, SigmaEps, ZetaEps, alphaFriction, alpha, &
                                    alphaTideABdec, alphaTideABmin, alphaTideAB
                    
    contains
    
    subroutine Const_init
        implicit none
        real            :: rnd_zeta, rnd_sigma
        real(kind=PREC) :: CFL
        
        if ( iargc() < 1 ) then
            write(*, *) 'Error: please input directory name for input data!'
            call exit(1)
        end if
        
        call getarg(1, case_dir)
        if ( iargc() > 1 ) call getarg(2, config_name) 
#ifdef __WIN32__
        case_dir = trim(case_dir) // '\'
        output_dir = trim(case_dir) // 'output\'
#else
        case_dir = trim(case_dir) // '/'
        output_dir = trim(case_dir) // 'output/'
#endif
        open(10, file=trim(case_dir)//config_name, status='old')
        read(10, nml=CONF)
        close(10)
        
        ! 分配内存空间
        call allocate_vars
        
        ! 初始化科里奥利力相关数据
        call init_coriolis
        
        ! 初始化掩码数据
        call init_masks
        
        ! 初始化对应的数据
        if ( len_trim(observe_data_file) > 0 ) call init_observe
        if ( len_trim(friction_file) > 0 ) call init_friction
        if ( len_trim(open_boundary_file) > 0 ) then 
            call init_open_boundary
        else
            ! 对于边界条件，需要将其随机初始化
            ! 非常像神经网络呢，这里肯定可以有改进的地方
            call RANDOM_NUMBER(rnd_zeta)
            call RANDOM_NUMBER(rnd_sigma)
            do x = 1, numX
                do y = 1, numY
                    if ( maskBoundary(x, y) > 1 ) then
                        newTide_A(x, y) = rnd_zeta * dcos(rnd_sigma * 2 * PI)
                        newTide_B(x, y) = rnd_zeta * dsin(rnd_sigma * 2 * PI)
                    end if
                end do
            end do
        end if
        
        ! 检查CFL条件
        CFL = dx(numY) / dsqrt(2. * g * maxval(WaterDepth))
        write(*, "(1X, 'dt =', F10.3, '  CFL =', F10.3)") dt, CFL
    end subroutine
    
    subroutine allocate_vars
        implicit none
        integer :: t
        
        allocate(f(numY))
        allocate(dx(numY))
        allocate(ddx(numY))
        
        allocate(WaterDepth(numX, numY))
        allocate(obsZeta(numX, numY))
        allocate(obsSigma(numX, numY))
        allocate(Friction(numX, numY))
        allocate(newFriction(numX, numY))
        allocate(Viscosity(numX, numY))
        allocate(Tide_A(numX, numY))
        allocate(Tide_B(numX, numY))
        allocate(newTide_A(numX, numY))
        allocate(newTide_B(numX, numY))
        
        allocate(maskH(numX, numY))
        allocate(maskU(numX, numY))
        allocate(maskV(numX, numY))
        allocate(maskBoundary(numX, numY))
        allocate(maskObserve(numX, numY))
        
        allocate(H0(numSteps, numX, numY))
        
        ! 初始化所有分配的空间
        maskU = 0
        maskV = 0
        maskH = 0
        maskBoundary = 0
        maskObserve = 0
        
        f = 0.
        H0 = 0.
        dx = 0.
        ddx = 0.
        
        WaterDepth = 0.
        obsZeta = 0.
        obsSigma = 0.
        Tide_A = 0.
        Tide_B = 0.
        newTide_A = 0.
        newTide_B = 0.
        
    end subroutine
    
    subroutine init_coriolis
        implicit none
        real(kind=PREC) :: theta
        
        ! convert resolution to degree
        Resolution = Resolution / 60.
        
        T = 2 * PI / TideOmega
        dt = T / (numSteps - 1)
        dy = 2 * PI * EarthRadius / 360. * Resolution
        ddy = dt / dy
        
        if ( Coriolis ) then
            ! 初始化网格宽度、科氏参数（随纬度变化）
            theta = ( Lat + Resolution / 2. ) * deg2rad
            do y = 1, numY
                dx(y) = dy * dcos(theta)
                ddx(y) = dt / dx(y)
                f(y) = 2 * Omega * dsin(theta)
                theta = theta + Resolution * deg2rad
            end do
        else
            dx(:) = dy
            ddx(:) = ddy
            f(y) = 0.d0
        end if
        
    end subroutine
        
    subroutine init_masks
        implicit none
        
        ! 读入水深
        open(10, file=trim(case_dir)//water_depth_file, status='old')
        do y = numY, 1, -1
            read(10, *) (WaterDepth(x, y),x=1,numX)
        end do
        close(10)
        
        ! 读入边界掩码
        open(10, file=trim(case_dir)//border_file, status='old')
        do y = numY, 1, -1
            read(10, '(<numX>I1)') (maskBoundary(x, y),x=1,numX)
        end do
        close(10)
        
        ! 计算水深掩码
        where ( WaterDepth > eps )
            maskH = 1
        end where
        
        ! 检查水深掩码与边界掩码是否一致
        do y = 1, numY
            do x = 1, numX
                if ( maskH(x, y) * maskBoundary(x, y) == 0 &
                    .and. maskH(x, y) + maskBoundary(x, y) /= 0 ) then
                    write(*, *) "Error: water depth do not match with border mask."
                    call exit(2)
                end if
            end do
        end do
        
        ! 计算流速U, V的掩码
        do y = 1, numY
            do x = 1, numX - 1
                if ( WaterDepth(x, y) * WaterDepth(x + 1, y) > eps ) maskU(x, y) = 1
            end do
        end do
        
        do y = 1, numY - 1
            do x = 1, numX
                if ( WaterDepth(x, y) * WaterDepth(x, y + 1) > eps ) maskV(x, y) = 1
            end do
        end do
        
        ! 初始化底摩擦、涡动粘性系数
        where ( maskH == 1 )
            Friction = globalFriction
            newFriction = globalFriction
            Viscosity = globalViscosity
        else where
            Friction = 0.
            newFriction = 0.
            Viscosity = 0.
        end where
        
    end subroutine
    
    subroutine init_observe
        implicit none
        real(kind=PREC) :: xx, yy, cur_zeta, cur_sigma, cur_lon, cur_lat
        integer         :: i
        logical         :: found
        
        ! 读入观测值数据，并对应到最近的网格点上面
        open(10, file=trim(case_dir)//observe_data_file, status='old')
        ! 两种读入观测数据的模式
        if ( observeFormat == 0 ) then
            ! 1. 将观测点的经纬度对应到最近的网格上面
            do i = 1, numObserveData
                read(10, *) xx, yy, cur_sigma, cur_zeta
                found = .false.
                do x = 1, numX
                    do y = 1, numY
                        cur_lon = Lon + (x - 1) * Resolution
                        cur_lat = Lat + (y - 1) * Resolution
                        if ( dabs(cur_lon - xx) < Resolution / 2. &
                        .and. dabs(cur_lat - yy) < Resolution / 2. ) then
                            obsZeta(x, y) = cur_zeta / 100.
                            obsSigma(x, y) = cur_sigma
                            found = .true.
                            exit
                        end if
                    end do
                    if ( found ) exit
                end do
            end do
        else if ( observeFormat == 1 ) then
            ! 2. 数据中包含的是就网格索引
            do i = 1, numObserveData
                read(10, *) x, y, cur_sigma, cur_zeta
                obsZeta(x, y) = cur_zeta / 100.
                obsSigma(x, y) = cur_sigma
            end do
        endif
        close(10)
        
        ! 设置观测点位置
        where ( obsZeta > eps .and. maskBoundary == 1 .and. maskH == 1 )
            maskObserve = 1
        end where
        
        OPEN(10,file=trim(output_dir)//'observe_mask.dat', status='replace')
        DO y = numY, 1, -1
            write(10,'(<numX>I3)') (maskObserve(x,y),x=1,numX)
        ENDDO
        CLOSE(10)
        
        call dump2D(real(maskBoundary, kind=PREC), 'boundary_mask.dat')
        
        ! 输出观测点总数
        sumObserve = sum(maskObserve)
        write(*, *) "Number of observe point: ", sumObserve
        
        ! 计算整个周期上观测点处的水位变化情况
        do t = 1, numSteps
            do x = 1, numX
                do y = 1, numY
                    H0(t, x, y) = obsZeta(x, y) * dcos(obsSigma(x, y) * deg2rad) * dcos(TideOmega*(t-1)*dt) &
                                + obsZeta(x, y) * dsin(obsSigma(x, y) * deg2rad) * dsin(TideOmega*(t-1)*dt)
                end do
            end do
        end do
        
    end subroutine
    
    subroutine init_open_boundary
        implicit none
        integer             :: i
        real(kind=PREC)     :: a, b
        
        open(10, file=trim(case_dir)//open_boundary_file, status='old')
        do i = numOpenBoundary, 1, -1
            read(10, *) x, y, a, b
            if ( maskBoundary(x, y) > 1 ) then
                Tide_A(x, y) = a
                Tide_B(x, y) = b
                newTide_A(x, y) = a
                newTide_B(x, y) = b
            else
                write(*, *) "Error: harmonic constant is not on boundary."
                call exit(3)
            end if
        end do
        close(10)
    end subroutine
    
    subroutine init_friction
        implicit none
        ! 直接读入文件中的底摩擦系数
        open(10, file=trim(case_dir)//friction_file, status='old')
        do y = numY, 1, -1
            read(10, *) (Friction(x, y),x=1,numX)
        end do
        close(10)
        newFriction = Friction
    end subroutine
    
    subroutine dump1D(array, filename, start, finish, length)
        implicit none
        integer         :: start, finish, length, i
        real(kind=PREC) :: array(length)
        character(len=*):: filename
        open(10, file=trim(output_dir)//filename, status='replace')
        write(10, "(<finish-start+1>F15.7)") (array(i), i = start, finish)
        close(10)
    end subroutine
    
    subroutine dump2D(array, filename)
        implicit none
        real(kind=PREC) :: array(numX, numY)
        character(len=*):: filename
        open(10, file=trim(output_dir)//filename, status='replace')
        do y = numY, 1, -1
            write(10, "(<numX>F15.7)") (array(x, y), x = 1, numX)
        end do
        close(10)
    end subroutine
    
    subroutine dump3D(array, filename)
        implicit none
        real(kind=PREC) :: array(numSteps, numX, numY)
        character(len=*):: filename
        open(10, file=trim(output_dir)//filename, status='replace')
        do t = 1, numSteps
            do y = numY, 1, -1
                write(10, "(<numX>F15.7)") (array(t, x, y), x = 1, numX)
            end do
        end do
        close(10)
    end subroutine

    end module