subroutine analyze( P, zeta, sigma, V, f, kappa, alpha, corr, H, g, S_0, M, Q )
    implicit none
    !数据总数、主分潮数，随从分潮数、中间值编号
    integer         :: M, P, Q
    !潮位数据，主、随从分潮角速度，天文初相角，交点因子，差比关系、
    real(kind=8)    :: zeta(M), sigma(P + Q), V(P + Q), f(P + Q), kappa(Q), alpha(Q)
    !随从分潮与主分潮的对应关系
    integer         :: corr(Q)
    !主分潮和随从分潮的中间变量η和ξ
    real(kind=8)    :: eta(P + Q), xi(P + Q)
    !方程的系数矩阵
    real(kind=8)    :: A(2*P + 1, 2*P + 1), b(2*P + 1), x(2*P + 1)
    !求解出的a和b，以及对前一次解的备份
    real(kind=8)    :: res_a(P + Q), res_b(P + Q), res_bak_a(P), res_bak_b(P), bak_S_0
    !调和常数H和g，以及平均海面S_0
    real(kind=8)    :: H(P + Q), g(P + Q), S_0
    !程序中的计数器等
    integer         :: i, j, n, mid
    !各种临时变量
    real(kind=8)    :: A_ij, B_ij, A_0i, A_0j, F_0, F_i, G_i, m_max
    !固定常数赋值
    real(kind=8)    :: pi = 3.14159265358979323846d0, eps = 1d-9

!下面是编译指导语句，用于给F2PY读取
!f2py intent(in,copy) P, zeta, sigma, V, f, kappa, alpha, corr
!f2py intent(out) H, g, S_0
!f2py integer,intent(hide),depend(zeta) :: M=len(zeta)
!f2py integer,intent(hide),depend(corr) :: Q=len(corr)

    !对η和ξ进行赋值
    do i = 1, P + Q
        eta(i) = f(i) * dsin(V(i))
         xi(i) = f(i) * dcos(V(i))
    end do

    !给随从分潮的a、b赋初始值，这里取零
    res_a(P + 1:P + Q) = 0.d0
    res_b(P + 1:P + Q) = 0.d0

    !记录前一次的解，用于确定主分潮的调和常数改变了多少；初始赋值为零
    bak_S_0 = 0.d0
    res_bak_a = 0.d0
    res_bak_b = 0.d0
    do while ( .true. )
        !对方程组左侧系数进行赋值
        !约定：前P个未知数为a_i，后P个未知数为b_i，最后一个为平均海面
        do i = 1, P
            A_0i = dsin(M*sigma(i)/2.d0) / dsin(sigma(i)/2.d0)
            A(i, 2*P + 1) = xi(i) * A_0i
            A(i + P, 2*P + 1) = eta(i) * A_0i
            do j = 1, P
                if ( i /= j ) then
                    A_ij = (dsin((sigma(i)-sigma(j))*M/2.d0) / dsin((sigma(i)-sigma(j))/2.d0) + &
                            dsin((sigma(i)+sigma(j))*M/2.d0) / dsin((sigma(i)+sigma(j))/2.d0)) / 2.d0
                    B_ij = (dsin((sigma(i)-sigma(j))*M/2.d0) / dsin((sigma(i)-sigma(j))/2.d0) - &
                            dsin((sigma(i)+sigma(j))*M/2.d0) / dsin((sigma(i)+sigma(j))/2.d0)) / 2.d0
                else
                    A_ij = ( M + dsin(M*sigma(i)) / dsin(sigma(i)) ) / 2.d0
                    B_ij = ( M - dsin(M*sigma(i)) / dsin(sigma(i)) ) / 2.d0
                end if
                A(i, j)         = xi(i) * xi(j) * A_ij + eta(i) * eta(j) * B_ij
                A(i, j + P)     = xi(i) * eta(j) * A_ij - eta(i) * xi(j) * B_ij
                A(i + P, j)     = eta(i) * xi(j) * A_ij - xi(i) * eta(j) * B_ij
                A(i + P, j + P) = eta(i) * eta(j) * A_ij + xi(i) * xi(j) * B_ij
            end do
        end do

        do j = 1, P
            A_0j = dsin(M*sigma(j)/2.d0) / dsin(sigma(j)/2.d0)
            A(2*P + 1, j) = xi(j) * A_0j
            A(2*P + 1, j + P) = eta(j) * A_0j
        end do
        A(2*P + 1, 2*P + 1) = M

        !计算数据的中间值，后面用于修正偏移
        mid = ( M - 1 ) / 2 + 1

        !对方程组右侧进行赋值
        F_0 = 0.d0
        do i = 1, M
            F_0 = F_0 + zeta(i)
        end do
        do i = 1, P
            F_i = 0.d0
            G_i = 0.d0
            do n = 1, M
                F_i = F_i + zeta(n) * dcos((n - mid)*sigma(i))
                G_i = G_i + zeta(n) * dsin((n - mid)*sigma(i))
            end do
            b(i) = xi(i) * F_i - eta(i) * G_i
            b(i + P) = eta(i) * F_i + xi(i) * G_i
            do j = P + 1, P + Q
                A_ij = (dsin((sigma(i)-sigma(j))*M/2.d0) / dsin((sigma(i)-sigma(j))/2.d0) + &
                        dsin((sigma(i)+sigma(j))*M/2.d0) / dsin((sigma(i)+sigma(j))/2.d0)) / 2.d0
                B_ij = (dsin((sigma(i)-sigma(j))*M/2.d0) / dsin((sigma(i)-sigma(j))/2.d0) - &
                        dsin((sigma(i)+sigma(j))*M/2.d0) / dsin((sigma(i)+sigma(j))/2.d0)) / 2.d0
                b(i)     = b(i)     - res_a(j) * ( xi(i) * xi(j) * A_ij + eta(i) * eta(j) * B_ij ) &
                                    - res_b(j) * ( xi(i) * eta(j) * A_ij - eta(i) * xi(j) * B_ij )
                b(i + P) = b(i + P) - res_a(j) * ( eta(i) * xi(j) * A_ij - xi(i) * eta(j) * B_ij ) &
                                    - res_b(j) * ( eta(i) * eta(j) * A_ij + xi(i) * xi(j) * B_ij )
            end do
        end do
        b(2*P + 1) = F_0
        do j = P + 1, P + Q
            A_0j = dsin(M*sigma(j)/2.d0) / dsin(sigma(j)/2.d0)
            b(2*P + 1) = b(2*P + 1) - res_a(j) * xi(j) * A_0j - res_b(j) * eta(j) * A_0j
        end do

        !调用高斯消元法模块求解线性方程组，方程中前P项赋值给a，后P项赋值给b，最后一项赋值给S_0
        call solve(A, b, x, 2*P + 1)
        res_a(1:P) = x(1:P)
        res_b(1:P) = x(P + 1:2 * P)
        S_0 = x(2*P + 1)

        m_max = dabs(S_0 - bak_S_0)
        do i = 1, P
            m_max = max(m_max, dabs(res_a(i) - res_bak_a(i)))
            m_max = max(m_max, dabs(res_b(i) - res_bak_b(i)))
        end do
        !如果迭代基本收敛，则退出循环，否则保存当前的解
        write(*,'("Current relative error: ", G17.9)') m_max
        if ( m_max < eps ) exit

        bak_S_0 = S_0
        res_bak_a = res_a(1:P)
        res_bak_b = res_b(1:P)

        !利用差比关系更新随从分潮的a、b
        do i = 1, Q
            res_a(i + P) = kappa(i) * ( res_a(corr(i)) * dcos(alpha(i)) - &
                                        res_b(corr(i)) * dsin(alpha(i)) )
            res_b(i + P) = kappa(i) * ( res_a(corr(i)) * dsin(alpha(i)) + &
                                        res_b(corr(i)) * dcos(alpha(i)) )
        end do
    end do

    !计算调和常数
    do i = 1, P + Q
        H(i) = dsqrt(res_a(i)**2 + res_b(i)**2)
        g(i) = datan(res_b(i) / res_a(i))
        !这一步是修正区时专用迟角g的取值范围，因为b和sin(g)显然应该同号
        if ( res_b(i) * dsin(g(i)) < 0.d0 ) g(i) = g(i) + pi
        if ( g(i) < 0.d0 ) g(i) = g(i) + 2 * pi
    end do

    return
end subroutine
