program main
    implicit none
    real,parameter ::  Ni = 0.0, Nf = -log(3001.0) ! e-folding numberの最初と最後を指定
    ! print *, Ni, Nf
    integer, parameter :: n = 2
    integer i
    real t, x1, x2, h, f

    t = 0.0
    x1 = 0.1
    x2 = 0.68
    h = 0.1

end program main

real function f(i, t, x1, x2)
    implicit none
    integer i
    real t, x1, x2
    if (i .eq. 1) then
        ! x1 := \Omega_r, x2 := \Omega_{\lambda} とする
        f = x1 * ( x1 - 3.0 * x2 - 1.0 ) ! \Omega_rに関する微分方程式
    else if (i .eq. 2) then
        f = x2 * ( x1 - 3.0 * x2 + 3.0 ) ! \Omega_{\lambda}に関する微分方程式
    end if
end function

subroutine runge4(t, x1, x2, n, h, f) ! Runge-Kutta methodを行うsubroutine
    implicit none 
    integer i, n
    real t, x1, x2, h, f
    real, allocatable, dimension(:)::k1, k2, k3, k4, work ! 任意次元の配列を用意
    allocate(k1(n), k2(n), k3(n), k4(n), work(n)) ! n次元を割り当てる

    do i = 1, n
        k1(i) = h * f(i, t, x(i))
        k2(i) = h * f(i, t + h/2.0, x(i) + k1(i)/2.0)
        k3(i) = h * f(i, t + h/2.0, x(i) + k2(i)/2.0)
        k4(i) = h * f(i, t + h, x(i) + k3(i))
        x(i) = x(i) + ( k1(i) + 2.0 * k2(i) + 2.0 * k3(i) + k4(i) )/6.0
    end do

    deallocate(k1, k2, k3, k4, work)
end subroutine runge4
