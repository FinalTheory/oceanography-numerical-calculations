subroutine predict( head, tail, num, zeta, sigma, V, f, H, g, S_0, mask, mask_num, delta )
    implicit none
    integer         :: head, tail, num, mask_num, i, j, idx
    real(kind=8)    :: zeta(tail - head + 1), sigma(num), V(num), f(num), H(num), g(num), S_0
    integer         :: mask(mask_num)
    real(kind=8)    :: t, delta

!f2py intent(in,copy) head, tail, sigma, V, f, H, g, mask, delta, S_0
!f2py intent(out) zeta
!f2py integer,intent(hide),depend(sigma) :: num=len(sigma)
!f2py integer,intent(hide),depend(mask) :: mask_num=len(mask)

    t = head * delta
    do i = 1, tail - head + 1
        zeta(i) = S_0
        do j = 1, mask_num
            idx = mask(j)
            zeta(i) = zeta(i) + f(idx) * H(idx) * dcos(sigma(idx) * t + V(idx) - g(idx))
        end do
        t = t + delta
    end do
    return
end subroutine
