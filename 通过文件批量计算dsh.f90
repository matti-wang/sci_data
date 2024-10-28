!---------------------------------------------------------------
!...Created on Tue Jul 16 01:03:45 2024
!...@author: PMO
!...@author: Dr. Jianghui JI  (jijh@pmo.ac.cn)
!...Purple Mountain Observatory, Chinese Academy of Sciences
!...Function: calculate the DSH of two orbits of 
!...Codes: f90 version 
!---------------------------------------------------------------
!   精度问题
program calculate_asteroid_orbit
    implicit none

    ! character(len=128) :: title
    character(len=100) :: in_file_name, out_file_name
    ! character    :: t, e
    integer      :: in_file_id = 10,  out_file_id = 20,  ios, itera
    integer      :: row_width = 86, current_row_num = -1
    real(kind=8) :: a_A, a_B, e_A, e_B, q_A, q_B, i_A, i_B, N_A, N_B, w_A, w_B
    real(kind=8) :: D_SH

    in_file_name = 'input_test.txt'
    out_file_name = './testtest.txt'
    open(unit=in_file_id, file=in_file_name, iostat=ios, status="old", action="read")
    if ( ios /= 0 ) stop "Error opening in_file"
    open(unit=out_file_id, file=out_file_name, iostat=ios, status='replace', action='write')
    if ( ios /= 0 ) stop "Error opening out_file"

    do
        itera = 0
        current_row_num = current_row_num + 1
        ! if (current_row_num >= max_row_num) exit
        call fseek(10, current_row_num * row_width, 0)

        read(10, *, iostat=ios) a_A, e_A, i_A, N_A, w_A
        if (ios /= 0) exit

        do while (itera <= current_row_num)
            write(out_file_id,'(A, 4X)', advance='no') 'XXX'
            itera = itera + 1
        end do

        do
            read(10, *, iostat=ios) a_B, e_B, i_B, N_B, w_B
            if (ios /= 0) then
                write(out_file_id, *) ''
                ! write(*,*) '============='
                exit
            end if

            D_SH = sqrt(calculate_D_SH(a_A, a_B, e_A, e_B, i_A, i_B, N_A, N_B, w_A, w_B))
            if (D_SH >= 1.0) then
                ! write(out_file_id,'(A, 4X)', advance='no') 'NaN'
                write(out_file_id,'(F6.4, X)', advance='no') D_SH
            else
                ! write(out_file_id,'(1X, F6.4)', advance='no') D_SH
                write(out_file_id,'(F6.4, 1X)', advance='no') D_SH
            end if
        end do
        rewind(10)
    end do

    close(unit=10, iostat=ios)
        if ( ios /= 0 ) stop "Error closing in_file"
    close(unit=out_file_id, iostat=ios)
        if ( ios /= 0 ) stop "Error closing out_file"

    write(*,*)
    write(*,*) 'Calculation finished!', out_file_name
    write(*,*)
    ! 100 format(21X, F5.4, 3X, F6.4, 4X, F5.2, 1X, F18.15, 1X, F18.16, 1X, F18.16)
    ! 100 format()

contains

    function calculate_D_SH(a_A_, a_B_, e_A_, e_B_, i_A_, i_B_, N_A_, N_B_, w_A_, w_B_) &
             result(D_SH_squared)
     
        implicit none

        real(kind=8), intent(in)  :: a_A_, e_A_, i_A_, N_A_, w_A_, &
                                     a_B_, e_B_, i_B_, N_B_, w_B_

        real(kind=8), parameter   :: pi = acos(-1.0d0)
        real(kind=8), parameter   :: deg2rad = pi / 180.0                               
        real(kind=8)              :: q_A,  i_A,  N_A,  w_A,  &
                                     q_B,  i_B,  N_B,  w_B
        real(kind=8)              :: D_SH_squared
        real(kind=8)              :: I_BA           ! 轨道平面夹角
        real(kind=8)              :: pi_BA          ! 自轨道相交处起量的近日点经度之差
        real(kind=8)              :: delta_Omega    ! 近日点经度之差
        real(kind=8)              :: half_sum_i     ! 轨道倾角的平均值
        real(kind=8)              :: cos_I_BA       ! 轨道平面夹角的余弦值
        real(kind=8)              :: sin_half_I_BA  ! 轨道平面夹角的一半的正弦值
        real(kind=8)              :: sin_half_pi_BA ! 自轨道相交处起量的近日点经度之差的一半的正弦值
        real(kind=8)              :: mean_e         ! 偏心率的平均值
        real(kind=8)              :: arcsin_term    ! arcsin项的值

        q_A = a_A_ * (1.d0 - e_A_)
        q_B = a_B_ * (1.d0 - e_B_)
        i_A = i_A_ * deg2rad
        i_B = i_B_ * deg2rad
        w_A = w_A_ * deg2rad
        w_B = w_B_ * deg2rad

        delta_Omega = N_B_ - N_A_
        if ( delta_Omega > 180 ) then
            delta_Omega = (delta_Omega - 360.d0) * deg2rad
        elseif ( delta_Omega < -180 ) then
            delta_Omega = (delta_Omega + 360.d0) * deg2rad
        else
            delta_Omega = delta_Omega * deg2rad
        end if

        cos_I_BA = cos(i_B) * cos(i_A) + sin(i_B) * sin(i_A) * cos(delta_Omega)
        I_BA = acos(cos_I_BA)

        arcsin_term = asin(cos((i_B + i_A) / 2) * sin(delta_Omega / 2) / &
                      cos(I_BA / 2))

        if (abs(delta_Omega) > pi) then
            pi_BA = w_B - w_A - 2 * arcsin_term
        else
            pi_BA = w_B - w_A + 2 * arcsin_term
        end if

        D_SH_squared = (e_B_ - e_A_) ** 2 + (q_B - q_A) ** 2 + &
                       4.0 * (sin(I_BA / 2)) ** 2 + &
                       (e_B + e_A) ** 2 * (sin(pi_BA / 2)) ** 2
        ! write(*,*) a_A_, e_A_, i_A_, N_A_, w_A_
        ! write(*,*) a_B_, e_B_, i_B_, N_B_, w_B_
        ! write(*,*) 
    end function calculate_D_SH

end program calculate_asteroid_orbit