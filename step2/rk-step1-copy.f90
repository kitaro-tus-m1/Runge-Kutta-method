! 独立変数2つ、従属変数1つの場合の連立1st-order ODE 初期値問題を解くためのRunge-Kutta methodを作る．
program main
    implicit none

    ! 解く方程式は x'(t) = x - 2xy =: f(x, y), y'(t) = -y + 3xy =: g(x, y)
    ! 初期値は x(0) = 1, y(0) = 2 とする
    ! t0 := t=0, x0 := x(0), y0 := y(0)
    ! しかし配列を使って，ベクトル表現で解くためには
    ! x(2) := (x, y), f(2) := (f, g), x0(2) := (x0, y0)とおくべき
    real,parameter :: t0 = 0.0, x0 = 1.0, y0 = 2.0

    ! 刻み幅h = 0.1とする
    real, parameter :: h = 0.1

    ! 独立変数tを定義し，初期値を代入しておく
    real :: t = t0
    ! 従属変数x, yを要素とする配列x(2)を定義する
    real x(2)
    ! 各要素に初期値を設定する
    data x /x0, y0/

    ! 繰り返し処理のためのiを定義する
    integer i
    ! 函数f1とf2を区別するためのjを定義する
    integer j
    ! Runge-Kutta method の各kを2要素配列として定義する
    real k1(2), k2(2), k3(2), k4(2)

    ! (done) 動作確認
    ! print *, x
    ! (done) 動作確認 配列xの第2成分を表示する
    ! print *, x(2)

    ! (done)動作確認
    !print *, f(10.0, 8.0)
    !print *, f(x, y)

    ! 初期値を画面に表示
    print *, t, x

    ! (done) 動作確認 函数fの初期値を画面に表示
    ! do i = 1, 2
    ! print *, f(i, t, x)
    ! end do

! 'rk-step2.dat'を開く
open(10, file='rk-step2.dat', status='replace')

! RUnge-Kutta methodの繰り返し処理を書く
do i = 1, 100
    ! f1, f2のそれぞれについて処理する
    do j = 1, 2
        k1(j) = h * f(j, t, x)
        k2(j) = h * f(j, t+h/2.0, x+k1(j)/2.0)
        k3(j) = h * f(j, t+h/2.0, x+k2(j)/2.0)
        k4(j) = h * f(j, t+h, x+k3(j))
        t = t + h
        x(j) = x(j) + (k1(j) + 2*k2(j) + 2*k3(j) + k4(j))/6.0
    end do
    ! 計算結果を画面に表示
    print *, t, x
    ! 計算結果を'rk-step2.dat'に出力
    write(10, *) t, x
end do
! 'rk-step2.dat'を閉じる
close(10)

stop

! 以下，内部副プログラム
contains
    ! 1解微分の右辺を関数f(j,t,x(1),x(2)) in { f1(t,x(1),x(2)), f2(t,x(1),x(2)) }として定義する
    real function f(j, t, x)
        implicit none
        real t, x(2)
        integer j
        if (j .eq. 1) then
            f = x(1) - 2 * x(1) * x(2)
        else if (j .eq. 2) then
            f = -x(2) + 3 * x(1) * x(2)
        end if
        return
    end function f


end program main