!资料统计部分，求平均值、标准差
subroutine Statistics(dat, n, aver, deviation)
    implicit none
    integer :: n, i
    real    :: dat(n)
    real    :: aver
    real    :: deviation
!f2py intent(in) dat
!f2py integer, intent(hide),depend(dat) :: n=len(dat)
!f2py intent(out) aver, deviation
    aver = sum(dat)
    aver = aver / real(n)
    deviation = 0.
    do i = 1, n
        deviation = deviation + ( dat(i) - aver ) ** 2
    end do
    deviation = deviation / real(n)
    deviation = sqrt(deviation)
end subroutine

!异常数据查找和处理
subroutine Find_Abnormal(dat, n)
    implicit none
    integer :: n
    real    :: dat(n)
    real    :: diff(n-1)
    integer :: flag(n-1)
    real    :: max_diff = 0
    integer :: cur_flag = 0
    integer :: start_pos = 1
    integer :: cur_len
    integer :: min_len
    integer :: found_pos_1, found_pos_2
    integer :: i
!f2py intent(in, out, inplace) dat
!f2py integer,intent(hide),depend(dat) :: n=len(dat)
    cur_len = n
    min_len = n

    diff = dat(2:n) - dat(1:n-1)
    where(diff < 0.)
        flag = -1
    elsewhere
        flag = 1
    end where
    diff = abs(diff)
    do i = 1, n - 1
        if ( diff(i) > max_diff ) then
            max_diff = diff(i)
            found_pos_1 = i
        end if
    end do
    do i = 1, n - 1
        if (cur_flag /= flag(i)) then
            if (cur_len < min_len) then
                found_pos_2 = start_pos + 1
                min_len = cur_len
            end if
            start_pos = i
            cur_flag = flag(i)
            cur_len = 1
        else
            cur_len = cur_len + 1
        end if
    end do
    if ( found_pos_1 == found_pos_2 ) then
        write(*, "('Index of abnormal position in data:', I5)") found_pos_1
        dat(found_pos_1) = ( dat(found_pos_1 - 1) + dat(found_pos_1 + 1) ) / 2.;
    end if
    return
end subroutine
