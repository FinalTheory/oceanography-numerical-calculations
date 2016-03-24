!一个简单的浮点数交换
subroutine swap( a, b )
    implicit none
    real(kind=8) :: a, b, tmp
    tmp = a
    a = b
    b = tmp
    return
end subroutine swap

!高斯消元法解线性方程组
subroutine solve( A, b, x, n )
    implicit none
    integer         :: idx, i, j, k, n
    real(kind=8)    :: A(n, n), b(n), x(n)
    real(kind=8)    :: ssum, mmax, m
    real(kind=8)    :: eps = 1e-8
!f2py intent(in) A, b
!f2py intent(out) x
!f2py integer,intent(hide),depend(b) :: n=len(b)
    do k = 1, n - 1
        ! find the max
        mmax = dabs(A(k, k))
        idx = k
        do i = k + 1, n
            if ( dabs(A(i, k)) > mmax ) then
                mmax= dabs(A(i, k))
                idx = i
            end if
        end do
        ! handle the fail condition
        if ( mmax < eps ) then
            return
        end if
        ! swap
        if ( idx /= k ) then
            do j = 1, n
                call swap(A(k, j), A(idx, j))
            end do
            call swap(b(k), b(idx))
        end if
        ! elimination
        do i = k + 1, n
            if ( dabs(A(k, k)) < eps ) then
                return
            end if
            m = A(i, k) / A(k, k)
            do j = k, n
                A(i, j) = A(i, j) - A(k, j) * m
            end do
            b(i) = b(i) - b(k) * m
        end do
    end do
    ! back substitution
    do i = n, 1, -1
        ssum = 0.d0
        do j = i + 1, n
            ssum = ssum + A(i, j) * b(j)
        end do
        b(i) = ( b(i) - ssum ) / A(i, i)
    end do
    x = b
    return
end subroutine solve
