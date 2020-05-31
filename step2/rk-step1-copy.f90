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
    ! 函数を配列として定義する場合は不要
    ! integer j
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
do i = 1, 200
        k1 = h * f(t, x)
        k2 = h * f(t+h/2.0, x+k1/2.0)
        k3 = h * f(t+h/2.0, x+k2/2.0)
        k4 = h * f(t+h, x+k3)
        t = t + h
        x = x + (k1 + 2*k2 + 2*k3 + k4)/6.0
    ! 計算結果を画面に表示
    print '(3e17.9)', t, x
    ! 計算結果を'rk-step2.dat'に出力
    write(10, '(3e17.9)') t, x
end do
! 'rk-step2.dat'を閉じる
close(10)

! 'rk-step2.plt'を作成してgnuplotでAqauaTermによるプロットを表示させる
! 'rk-step2.plt'を開く
open(11, file='rk-step2.plt', status='replace')
    write (11, '(a)') 'plot "rk-step2.dat" using 1:2 w l'
    write (11, '(a)') 'replot "rk-step2.dat" using 1:3 w l'
! 'rk-step2.plt'を閉じる
close(11)

! gnuplotを起動し，'rk-step1.plt'を実行する
call execute_command_line('gnuplot "rk-step2.plt"')


! 'rk-step1-save.plt'を作成してgnuplotでpng画像としてプロットを保存する
! 'rk-step1-save.plt'を開く
open(12, file='rk-step2-save.plt', status='replace')
    ! 出力先をpngに設定
    write (12, '(a)') 'set terminal png'
    ! 出力ファイルを'rk-step1.png'に設定
    write (12, '(a)') 'set output "rk-step2.png"'
    ! 'rk-step2.dat'をプロット
    write (12, '(a)') 'plot "rk-step2.dat" using 1:2 w l'
    ! グラフを重ねるため再度，出力ファイルを'rk-step2.png'に設定
    write (12, '(a)') 'set output "rk-step2.png"'
    ! 1列目(t)と3列目(x(2))を使ってグラフを描く
    write (12, '(a)') 'replot "rk-step2.dat" using 1:3 w l'
! 'rk-step1-save.plt'を閉じる
close(12)

! gnuplotを起動し，'rk-step1-save.plt'を実行する
call execute_command_line('gnuplot "rk-step2-save.plt"')

stop

! 以下，内部副プログラム
contains
    ! 1解微分の右辺を関数f(j,t,x(1),x(2)) in { f1(t,x(1),x(2)), f2(t,x(1),x(2)) }として定義する
    ! 2つの函数を2成分をもつ1つの配列として定義できるかやってみる
    function f(t, x)
        implicit none
        real f(2)
        real t, x(2)
        !integer j
        !if (j .eq. 1) then
            f(1) = x(1) - 2 * x(1) * x(2)
        !else if (j .eq. 2) then
            f(2) = -x(2) + 3 * x(1) * x(2)
        !end if
        return
    end function f


end program main