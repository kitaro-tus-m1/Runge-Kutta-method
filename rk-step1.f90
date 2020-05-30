! 独立変数1つ、従属変数1つの場合の1st-order ODE 初期値問題を解くためのRunge-Kutta methodを作る．
program main
    implicit none

    ! 解く方程式はy'(x) = x-yで，初期値はy(0) = 1
    ! x0 := x=0, y0 := y(x0) = y(0)
    real,parameter :: x0 = 0, y0 = 1

    ! 刻み幅h = 0.1とする
    real, parameter :: h = 0.1

    ! 独立変数x, 従属変数yを実数として定義し，初期値を代入しておく
    real :: x = x0, y = y0

    integer i
    real k1, k2, k3, k4
    ! real f

    ! (done)動作確認
    !print *, x, y

    ! 注意点
    ! 内部副プログラムで定義した函数を利用する際には，containsより前で呼び出さなければならない
    ! 内部副プログラムで函数を定義した場合，
    ! 14行目のような "real f" といった函数のデータ型指定は必要ない

    ! (done)動作確認
    !print *, f(10.0, 8.0)
    !print *, f(x, y)

    do i = 1, 100
        k1 = h * f(x, y)
        k2 = h * f(x+h/2.0, y+k1/2.0)
        k3 = h * f(x+h/2.0, y+k2/2.0)
        k4 = h* f(x+h, y+k3)
        x = x + h
        y = y + (k1 + 2*k2 + 2*k3 + k4)/6.0
        print *, x, y
    end do

    stop

! 以下，内部副プログラム
contains
    ! 1解微分の右辺を関数f(x,y)として定義する
    real function f(a, b)
        implicit none
        real a, b
        f = a - b
        return
    end function f

end program main