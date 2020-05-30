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
    !real f

    ! done 動作確認
    !print *, x, y

    stop

contains
    ! 1解微分の右辺を関数f(x,y)として定義する
real function f(a, b)
implicit none
real a, b
f = a - b
return
end function f

    print *, f(10.0 , 6.0)
    !print *, f(x, y)
end program main