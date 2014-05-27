! title: 二维潮波数值模拟程序
! brief: 本Fortran程序为命令行后端，用于实现数值模拟运算以及数据输出的功能；
!        如果单独运行本程序，需要以正确的参数按照一定顺序来调用，具体见下面的说明。
! Author: 黄龑
! date: 2014-05-21
! gfortran version: 4.5.2

! 变量命名约定:
! 1) 所有符号尽可能按照公式中的形式输入；
! 2) 英文字母以其原有形式输入；
! 3) 希腊字母等非ASCII字符以其LaTex缩写表示；
! 4) 以下划线后面的字符来表示下角标；
! 5) 其他有意义的变量遵循Javascript语言的命名规则。

! 程序调用参数顺序约定:
! 1.  水深数据文件名
! 2.  开边界调和常数数据文件名
! 3.  x方向节点数
! 4.  y方向节点数
! 5.  时间步数
! 6.  空间分辨率（单位:分）
! 7.  总周期时间（单位:秒）
! 8.  是否允许平流项（0表示否，1表示是）
! 9.  是否允许粘度项（0表示否，1表示是）
! 10. 是否允许底摩擦（0表示否，1表示是）
! 11. 是否使用科氏力参数（0表示否，1表示常数科氏力，2表示可变科氏力）
! 12. 是否使用随着纬度变化的x方向空间步长（0表示否，1表示是）
! 13. 流体粘度系数
! 14. 底摩擦系数
! 15. 半隐半显差分格式的系数α
! 16. 模拟区域的底部对应纬度
! 17. 输出文件的编号

module constants
    implicit none

    !--------------------常量设定部分--------------------

    ! 数值模拟的最大迭代次数
    integer, parameter :: maxLoops = 30
    ! 迭代精度控制变量、正负精度控制变量、常数π
    real(kind=8), parameter :: eps = 1d-2, epsZero = 1e-4, PI = 3.14159265358979323846d0, conv = PI / 180.d0
    ! 地球半径（米）、自转角速度（rad/s）、重力加速度(m/s^2)
    real(kind=8), parameter :: earthRadius = 6378.1d3, Omega = 2.d0 * PI / ( 24.d0 * 3600.d0 ), g = 9.8d0
    ! 时间步数为numSteps、开边界节点数量为numOpenNodes
    integer :: numX, numY, numSteps, numOpenNodes

    !--------------------输入参数声明--------------------

    ! 定义空间最小分辨率以及研究区域底边的纬度
    real(kind=8) :: resolution, latitude
    ! 定义时间步长、潮汐过程周期总长度、潮汐角速度等参数
    real(kind=8) :: delta_t, T, w
    ! 定义粘度系数、底摩擦系数、半隐半显格式的控制系数
    real(kind=8) :: A, K, alpha
    ! 定义开边界上的调和常数、水深数据h
    real(kind=8), allocatable :: Boundary(:, :), h(:, :)

    !--------------------运行参数声明--------------------

    ! 设定是否允许平流项、粘度项、底摩擦、科氏力
    integer :: enableAdvection, enableViscosity, enableFriction, enableCoriolisForce
    ! 设定是否使用自动变化的x方向空间步长；
    integer :: autoXStep

    !--------------------临时变量声明--------------------

    integer :: n, i, j, idx

    contains

    subroutine parseParameters( fileidx )
        implicit none
        integer :: ioStatus = 0, fileidx
        character(len=80) :: buf, name1, name2

        if ( iargc() /= 17 ) then
            write(*, *) 'Not enough input parameters!'
            call exit(-1)
        end if

        call getarg(1, buf)
        read(buf, *) name1
        call getarg(2, buf)
        read(buf, *) name2
        call getarg(3, buf)
        read(buf, *) numX
        call getarg(4, buf)
        read(buf, *) numY
        call getarg(5, buf)
        read(buf, *) numSteps
        call getarg(6, buf)
        read(buf, *) resolution
        ! 分辨率由"分"转化到"度"
        resolution = resolution / 60.d0
        call getarg(7, buf)
        read(buf, *) T
        ! 计算分潮角速度
        w = 2.d0 * PI / T
        call getarg(8, buf)
        read(buf, *) enableAdvection
        call getarg(9, buf)
        read(buf, *) enableViscosity
        call getarg(10, buf)
        read(buf, *) enableFriction
        call getarg(11, buf)
        read(buf, *) enableCoriolisForce
        call getarg(12, buf)
        read(buf, *) autoXStep
        call getarg(13, buf)
        read(buf, *) A
        call getarg(14, buf)
        read(buf, *) K
        call getarg(15, buf)
        read(buf, *) alpha
        call getarg(16, buf)
        read(buf, *) latitude
        call getarg(17, buf)
        read(buf, *) fileidx

        ! 在栈内存上分配数据空间
        allocate(h(numX, numY))
        allocate(Boundary(numX*numY, 4))

        delta_t = T / numSteps

        ! 读入水深数据和开边界数据
        open(11, file = name1, action = 'read')
        open(12, file = name2, action = 'read')
        read(11, *) ( h(:, j), j = numY, 1, -1 )
        i = 1
        do while ( .true. )
            read(12, *, iostat = ioStatus) Boundary(i, :)
            if ( ioStatus /= 0 ) exit
            i = i + 1
        end do
        numOpenNodes = i - 1
        close(11)
        close(12)

    end subroutine

end module

module calculate
    use constants
    implicit none

    !--------------------数据分配部分--------------------

    ! 定义海面波动，流速u、v
    real(kind=8), allocatable :: u(:, :, :), v(:, :, :), zeta(:, :, :)
    ! 定义随维度不同而发生变化的x方向空间步长、固定的y方向空间步长、科氏力参数
    real(kind=8), allocatable :: delta_x(:), delta_y, f(:)
    ! 定义开边界上的调和常数数据以及水位数据
    real(kind=8), allocatable :: boundaryWaterLevel(:, :, :)
    ! 定义计算区域掩码以及流速u、v的掩码
    integer, allocatable :: mask(:, :), ctrlU(:, :), ctrlV(:, :)
    ! 记录方程左侧系数以及右侧数值、平均速度等临时变量
    real(kind=8) :: equalLeft, equalRight, temp, u_aver, v_aver

    contains

    ! 计算在一个周期中边界条件上的水位情况
    subroutine calculateBoundary()
        implicit none
        ! 定义振动强迫当前的相位
        real(kind=8) :: curPhase

        boundaryWaterLevel = 0.d0

        do n = 1, numSteps + 1
            do idx = 1, numOpenNodes
                curPhase = w * ( n - 1 ) * delta_t
                i = int(Boundary(idx, 1) + epsZero)
                j = int(Boundary(idx, 2) + epsZero)
                boundaryWaterLevel(i, j, n) = Boundary(idx, 3) * dcos(curPhase) + Boundary(idx, 4) * dsin(curPhase)
                ! write(*, *) boundaryWaterLevel(i, j, n)
            end do
        end do

    end subroutine

    subroutine initData()
        implicit none
        real(kind=8) :: theta

        allocate(f(numY))
        allocate(delta_x(numY))
        allocate(mask(numX, numY))
        allocate(ctrlU(numX, numY))
        allocate(ctrlV(numX, numY))
        allocate(u(numX, numY, numSteps + 1))
        allocate(v(numX, numY, numSteps + 1))
        allocate(zeta(numX, numY, numSteps + 1))
        allocate(boundaryWaterLevel(numX, numY, numSteps + 2))

        ! 初始纬度选取为模拟区域的中位线
        theta = ( numY / 2 * resolution + latitude ) * conv
        delta_y = 2.d0 * PI * earthRadius * resolution / 360.d0
        delta_x(:) = dcos(theta) * delta_y
        if ( enableCoriolisForce == 0 ) then
            f(:) = 0.d0
        else if ( enableCoriolisForce == 1 ) then
            f(:) = 2.d0 * Omega * dsin(theta)
        end if

        if ( autoXStep /= 0 ) then
            theta = ( resolution / 2.d0 + latitude ) * conv
            do j = 1, numY
                delta_x(j) = dcos(theta) * delta_y
                theta = theta + resolution * conv
            end do
        end if

        if ( enableCoriolisForce == 2 ) then
            theta = ( resolution / 2.d0 + latitude ) * conv
            do j = 1, numY
                f(j) = 2.d0 * Omega * dsin(theta)
                theta = theta + resolution * conv
            end do
        end if

        ! 根据水深创建mask
        where ( h > epsZero )
            mask = 1
        elsewhere
            mask = 0
        endwhere

        ! 强制令开边界上掩码为1
        do idx = 1, numOpenNodes
            ! 浮点误差修正
            i = int(Boundary(idx, 1) + epsZero)
            j = int(Boundary(idx, 2) + epsZero)
            mask(i, j) = 5
        end do

        ! 计算流速掩码
        ctrlU = 1
        ctrlV = 1
        do j = 1, numY
            do i = 1, numX
                if ( mask(i, j) == 0 .or. mask(i + 1, j) == 0 ) ctrlU(i, j) = 0
                if ( mask(i, j) == 0 .or. mask(i, j + 1) == 0 ) ctrlV(i, j) = 0
            end do
        end do

!        write(*, *) f
!        write(*, *) delta_x
!        write(*, *) delta_t
!        write(*, *) resolution, latitude
!        write(*, *) numX, numY, numSteps, numOpenNodes
!        write(*, *) Boundary(:, :)

        call calculateBoundary()

        ! 将水位、流速初始化为零
        ! 并将初始时刻的水位初始化为开边界水位
        u = 0.d0
        v = 0.d0
        zeta = 0.d0
        zeta(:, :, numSteps + 1) = boundaryWaterLevel(:, :, 1)

    end subroutine

    ! 根据当前步迭代下一步的水位
    subroutine calculateZeta(curStep)
        implicit none
        integer :: curStep
        real(kind=8) :: h_r, h_l, h_u, h_d

        do j = 1, numY
            do i = 1, numX
                ! 判断这个位置是否需要计算水位，以及是否为开边界
                if ( mask(i, j) == 0 ) then
                    zeta(i, j, curStep + 1) = 0.d0
                    cycle
                else if ( mask(i, j) == 5 ) then
                    zeta(i, j, curStep + 1) = boundaryWaterLevel(i, j, curStep + 1)
                    cycle
                end if

                h_r = ( h(i, j) + h(i + 1, j) + zeta(i, j, curStep) + zeta(i + 1, j, curStep) ) * 0.5d0
                h_l = ( h(i, j) + h(i - 1, j) + zeta(i, j, curStep) + zeta(i - 1, j, curStep) ) * 0.5d0
                h_u = ( h(i, j) + h(i, j + 1) + zeta(i, j, curStep) + zeta(i, j + 1, curStep) ) * 0.5d0
                h_d = ( h(i, j) + h(i, j - 1) + zeta(i, j, curStep) + zeta(i, j - 1, curStep) ) * 0.5d0

                zeta(i, j, curStep + 1) = zeta(i, j, curStep) - delta_t * &
                ( ( h_r * u(i, j, curStep) - h_l * u(i-1, j, curStep) ) / delta_x(j) + &
                ( h_u * v(i, j, curStep) - h_d * v(i, j-1, curStep) ) / delta_y )
            end do
        end do
    end subroutine

    ! 根据参数计算下一时间步的u
    ! 两个参数分别代表当前时间序号（以推出下一步）、计算时所依赖的v的时间序号
    subroutine calculateU(curStep, relyOnStep)
        implicit none
        integer :: curStep, relyOnStep

        do j = 1, numY
            do i = 1, numX
                if ( ctrlU(i, j) == 0 ) then
                    u(i, j, curStep + 1) = 0.d0
                    cycle
                end if

                v_aver = ( v(i, j, relyOnStep) + v(i + 1, j, relyOnStep) + &
                v(i, j - 1, relyOnStep) + v(i + 1, j - 1, relyOnStep) ) / &
                max(dble(ctrlV(i, j) + ctrlV(i + 1, j) + ctrlV(i, j - 1) + ctrlV(i + 1, j - 1)), 1.d0)

                ! 初始化方程两端
                equalLeft = 1.d0 / delta_t
                equalRight = u(i, j, curStep) / delta_t - g * ( zeta(i + 1, j, curStep + 1) - zeta(i, j, curStep + 1) ) / delta_x(j)

                ! 加入科氏力
                if ( enableCoriolisForce /= 0 ) then
                    equalRight = equalRight + f(j) * v_aver
                end if

                ! 加入黏性项
                if ( enableViscosity /= 0 ) then
                    equalRight = equalRight + A * ( ( u(i + 1, j, curStep) - 2 * u(i, j, curStep) + u(i - 1, j, curStep) ) &
                    / delta_x(j)**2 + ( u(i, j + 1, curStep) - 2 * u(i, j, curStep) + u(i, j - 1, curStep) ) / delta_y**2 )
                end if

                ! 加入平流项
                if ( enableAdvection /= 0 ) then
                    equalRight = equalRight - u(i, j, curStep) * ( u(i + 1, j, curStep) - u(i - 1, j, curStep) ) &
                    / delta_x(j) / 2.d0 - v_aver * ( u(i, j + 1, curStep) - u(i, j - 1, curStep) ) / delta_y / 2.d0
                end if

                ! 加入底摩擦
                if ( enableFriction /= 0 ) then
                    temp = 2.d0 * K * dsqrt( u(i, j, curStep)**2 + v_aver**2 ) / ( h(i, j) + h(i + 1, j) + &
                    zeta(i, j, curStep) + zeta(i + 1, j, curStep) )
                    equalLeft = equalLeft + ( 1.d0 - alpha ) * temp
                    equalRight = equalRight - temp * alpha * u(i, j, curStep)
                end if
                u(i, j, curStep + 1) = equalRight / equalLeft
            end do
        end do
    end subroutine

    subroutine calculateV(curStep, relyOnStep)
        implicit none
        integer :: curStep, relyOnStep

        do j = 1, numY
            do i = 1, numX
                if ( ctrlV(i, j) == 0 ) then
                    v(i, j, curStep + 1) = 0.d0
                    cycle
                end if

                u_aver = ( u(i, j, relyOnStep) + u(i, j + 1, relyOnStep) + &
                u(i - 1, j, relyOnStep) + u(i - 1, j + 1, relyOnStep) ) / &
                max(dble(ctrlU(i, j) + ctrlU(i, j + 1) + ctrlU(i - 1, j) + ctrlU(i - 1, j + 1)), 1.d0)

                ! 初始化方程两端
                equalLeft = 1.d0 / delta_t
                equalRight = v(i, j, curStep) / delta_t - g * ( zeta(i, j + 1, curStep + 1) - zeta(i, j, curStep + 1) ) / delta_y

                ! 加入科氏力
                if ( enableCoriolisForce /= 0 ) then
                    equalRight = equalRight - f(j) * u_aver
                end if

                ! 加入黏性项
                if ( enableViscosity /= 0 ) then
                    equalRight = equalRight + A * ( ( v(i + 1, j, curStep) - 2 * v(i, j, curStep) + v(i - 1, j, curStep) ) &
                    / delta_x(j)**2 + ( v(i, j + 1, curStep) - 2 * v(i, j, curStep) + v(i, j - 1, curStep) ) / delta_y**2 )
                end if

                ! 加入平流项
                if ( enableAdvection /= 0 ) then
                    equalRight = equalRight - u_aver * ( v(i + 1, j, curStep) - v(i - 1, j, curStep) ) / delta_x(j) / 2.d0 &
                    - v(i, j, curStep) * ( v(i, j + 1, curStep) - v(i, j - 1, curStep) ) / delta_y / 2.d0
                end if

                ! 加入底摩擦
                if ( enableFriction /= 0 ) then
                    temp = 2.d0 * K * dsqrt( v(i, j, curStep)**2 + u_aver**2 ) / ( h(i, j) + h(i, j + 1) + &
                    zeta(i, j, curStep) + zeta(i, j + 1, curStep) )
                    equalLeft = equalLeft + ( 1.d0 - alpha ) * temp
                    equalRight = equalRight - temp * alpha * v(i, j, curStep)
                end if
                v(i, j, curStep + 1) = equalRight / equalLeft
            end do
        end do
    end subroutine

    ! 计算调和常数并输出
    subroutine outputHarmonic()
        implicit none
        character(len = 80) :: str
        real(kind=8), dimension(numX, numY) :: harmonicA, harmonicB, amplitude, arg

        write(str, '(I0)') numX

        harmonicA = 0.d0
        harmonicB = 0.d0

        do n = 1, numSteps
            do j = 1, numY
                do i = 1, numX
                    harmonicA(i, j) = harmonicA(i, j) + zeta(i, j, n) * dcos(w*(n-1)*delta_t)
                    harmonicB(i, j) = harmonicB(i, j) + zeta(i, j, n) * dsin(w*(n-1)*delta_t)
                end do
            end do
        end do

        harmonicA = harmonicA / dble(numSteps) * 2.d0
        harmonicB = harmonicB / dble(numSteps) * 2.d0

        do j = 1, numY
            do i = 1, numX
                amplitude(i, j) = dsqrt(harmonicA(i, j)**2 + harmonicB(i, j)**2)
                ! 判断是否是无潮点
                if ( amplitude(i, j) < epsZero ) then
                    arg(i, j) = -1.1d5
                ! 判断是否是X轴正方向
                else if ( dabs(harmonicB(i, j)) < epsZero .and. harmonicA(i, j) >  amplitude(i, j) - epsZero ) then
                    arg(i, j) = 0.d0
                ! 判断是否是X轴负方向
                else if ( dabs(harmonicB(i, j)) < epsZero .and. harmonicA(i, j) < -amplitude(i, j) + epsZero ) then
                    arg(i, j) = PI
                ! 判断是否是Y轴正方向
                else if ( dabs(harmonicA(i, j)) < epsZero .and. harmonicB(i, j) <  amplitude(i, j) - epsZero ) then
                    arg(i, j) = PI / 2.d0
                ! 判断是否是Y轴负方向
                else if ( dabs(harmonicA(i, j)) < epsZero .and. harmonicB(i, j) < -amplitude(i, j) + epsZero ) then
                    arg(i, j) = PI / 2.d0 * 3.d0
                ! 判断是否是第一象限
                else if ( harmonicB(i, j) > 0.d0 .and. harmonicA(i, j) > 0.d0 ) then
                    arg(i, j) = datan(harmonicB(i, j) / harmonicA(i, j))
                ! 判断是否是第二象限
                else if ( harmonicB(i, j) > 0.d0 .and. harmonicA(i, j) < 0.d0 ) then
                    arg(i, j) = datan(harmonicB(i, j) / harmonicA(i, j)) + PI
                ! 判断是否是第三象限
                else if ( harmonicB(i, j) < 0.d0 .and. harmonicA(i, j) < 0.d0 ) then
                    arg(i, j) = datan(harmonicB(i, j) / harmonicA(i, j)) + PI
                ! 判断是否是第四象限
                else if ( harmonicB(i, j) < 0.d0 .and. harmonicA(i, j) > 0.d0 ) then
                    arg(i, j) = datan(harmonicB(i, j) / harmonicA(i, j))
                end if
            end do
        end do

        arg = arg / PI * 180.d0

        do j = 1, numY
            do i = 1, numX
                if ( mask(i, j) == 0 ) then
                    amplitude(i, j) = -1.1d5
                    arg(i, j) = -1.1d5
                end if
            end do
        end do

        open(14, file = 'amplitude.txt', status = 'replace', action = 'write')
        do j = numY, 1, -1
            write(14, '(' // str // '(1X, F15.6))') ( amplitude(i, j), i = 1, numX )
        end do
        close(14)
        open(14, file = 'arg.txt', status = 'replace', action = 'write')
        do j = numY, 1, -1
            write(14, '(' // str // '(1X, F15.6))') ( arg(i, j), i = 1, numX )
        end do
        close(14)

    end subroutine

end module

subroutine outputData( filename, matrix, numX, numY, numSteps )
    implicit none
    integer :: i, j, n, numX, numY, numSteps
    real(kind=8) :: matrix(numX, numY, numSteps)
    character(len = 80) :: filename, str
    write(str, '(I0)') numX

    open(13, file = filename, status = 'replace', action = 'write')

    do n = 1, numSteps
        do j = numY, 1, -1
            write(13, '(' // str // '(1X, F15.6))') ( matrix(i, j, n), i = 1, numX )
        end do
    end do

    close(13)

end subroutine

program TideModeling
    use calculate
    implicit none
    integer :: loopCounter, fileidx
    ! 定义迭代精度控制变量
    real(kind=8) :: accuracy
    character(len=80) :: filename, str

    call parseParameters(fileidx)
    call initData()

    do loopCounter = 1, maxLoops

        ! 将上一个周期的结束状态赋值给当前周期的起始状态
        u(:, :, 1) = u(:, :, numSteps + 1)
        v(:, :, 1) = v(:, :, numSteps + 1)
        zeta(:, :, 1) = zeta(:, :, numSteps + 1)

        ! 开始本周期的模拟
        do n = 1, numSteps - 1, 2
            call calculateZeta(n)
            call calculateV(n, n)
            call calculateU(n, n + 1)
            ! 交错步迭代
            call calculateZeta(n + 1)
            call calculateU(n + 1, n + 1)
            call calculateV(n + 1, n + 2)
        end do

        accuracy = maxval(dabs(zeta(:, :, numSteps + 1) - zeta(:, :, 1)))
        write(*, '("Current error: ", F15.6)') accuracy
        if ( accuracy < eps ) exit

    end do

    do j = 1, numY
        do i = 1, numX
            if ( ctrlU(i, j) == 0 ) u(i, j, :) = -1.1d5
            if ( ctrlV(i, j) == 0 ) v(i, j, :) = -1.1d5
            if ( mask(i, j) == 0 ) zeta(i, j, :) = -1.1d5
        end do
    end do

    write(str, '(I3.3)') fileidx
    filename = 'output_U_' // trim(str) // '.txt'
    call outputData(filename, u, numX, numY, numSteps)
    filename = 'output_V_' // trim(str) // '.txt'
    call outputData(filename, v, numX, numY, numSteps)
    filename = 'output_Zeta_' // trim(str) // '.txt'
    call outputData(filename, Zeta, numX, numY, numSteps)

    ! 计算调和常数并输出
    call outputHarmonic()

end program
