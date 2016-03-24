!**
! @file main.f90
! @brief  仅供调试时使用
! @author 黄龑
! @version 0.1
! @date 2014-04-06
!**
program main
    implicit none
    integer, parameter :: P = 11, Q = 2, M = 743
    real(kind=8)    :: sigma(P + Q), V(P + Q), f(P + Q), kappa(Q), alpha(Q)
    real(kind=8)    :: H(P + Q), g(P + Q), S_0
    integer         :: corr(Q), mask(6), i
    real(kind=8)    :: zeta(M), zeta_new(742*60 + 1)
    open(11, file = 'data_real.txt', action = 'read')
    read(11, *) ( zeta(i), i = 1, M )
    open(12, file = 'sigma.txt', action = 'read')
    read(12, *) ( sigma(i), i = 1, P + Q )
    open(13, file = 'V.txt', action = 'read')
    read(13, *) ( V(i), i = 1, P + Q )

    open(14, file = 'f.txt', action = 'read')
    read(14, *) ( f(i), i = 1, P + Q )

    open(15, file = 'kappa.txt', action = 'read')
    read(15, *) ( kappa(i), i = 1, Q )
    open(16, file = 'alpha.txt', action = 'read')
    read(16, *) ( alpha(i), i = 1, Q )
    open(17, file = 'corr.txt', action = 'read')
    read(17, *) ( corr(i), i = 1, Q )
    open(18, file = 'mask.txt', action = 'read')
    read(18, *) ( mask(i), i = 1, 6 )

    ! call Find_Abnormal(dat, 62*12)
    call analyze( P, zeta, sigma, V, f, kappa, alpha, corr, H, g, S_0, M, Q )
    write(*, *) S_0

!    open(19, file = 'H.txt', action = 'write')
!    write(19, *) ( H(i), i = 1, P + Q )
!    open(20, file = 'g.txt', action = 'write')
!    write(20, *) ( g(i), i = 1, P + Q )
!    open(21, file = 'S_0.txt', action = 'write')
!    write(21, *) S_0

    call predict( -371 * 60, 371 * 60, 13, zeta_new, sigma, V, f, H, g, S_0, mask, 6, 1.d0 / 60.d0 )
    open(22, file = 'data_predict.txt', action = 'write')
    write(22, *) ( zeta_new(i), i = 1, 742*60 + 1 )
end
