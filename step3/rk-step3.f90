! 独立変数2つ、従属変数1つの場合の連立1st-order ODE 初期値問題を解くためのRunge-Kutta methodを作る．
program main
    implicit none
    ! 表記を簡単にするため，x := \Omega_r, y := \Omega_{\lambda}, t := N とする
    ! 解く方程式は x'(t) = x (x - 3y - 1) =: f(x, y), y'(t) = y (x - 3y + 3) =: g(x, y)
    ! 初期値は x(0) = unknown, y(0) = 0.68
    ! 境界条件は y(-Ln(3001)) = 1 - 2 x(-Ln(3001))
    ! x(0)がunknownなので，x(0)を色々と仮定して解いてみる
    ! その試行の中で，境界条件を満たすものが現れたら，それが解であるとみなす
    ! t = -Ln(3001) が計算結果の出力ファイルの何行目にあるかを把握しなければならない

    ! t0 := t=0, x0 := x(0), y0 := y(0)
    ! 配列を使って，ベクトル表現で解くためには
    ! x(2) := (x, y), f(2) := (f, g), x0(2) := (x0, y0)とおくべき
    ! 1+z = 0.1 から計算を始めるためにt0は次のようにした
    double precision,parameter :: t0 = -Log(0.10d0), x0 = 0.000000015d0, y0 = 0.99955d0
    ! x0, y0 を色々と調節した

    ! 刻み幅h = -0.1とする
    double precision, parameter :: h = -0.1d0

    ! 独立変数tを定義し，初期値を代入しておく
    double precision :: t = t0
    ! 従属変数x, yを要素とする配列x(2)を定義する
    double precision x(2)
    ! 各要素に初期値を設定する
    data x /x0, y0/

    ! 繰り返し処理のためのiを定義する
    integer i
    ! 函数f1とf2を区別するためのjを定義する
    ! 函数を配列として定義する場合は不要
    ! integer j
    ! Runge-Kutta method の各kを2要素配列として定義する
    double precision k1(2), k2(2), k3(2), k4(2)

    ! 初期値を画面に表示
    print *, Exp(-t), x, 1.0d0 - x(1) - x(2)

! 'rk-step3.dat'を開く
open(10, file='rk-step3.dat', status='replace')
    ! 初期値を'rk-step3.dat'に出力
    write(10, '(3e17.9)') Exp(-t), x, 1.0d0 - x(1) - x(2)

! Runge-Kutta methodの繰り返し処理を書く
do i = 1, 180
        k1 = h * f(t, x)
        k2 = h * f(t+h/2.0, x+k1/2.0)
        k3 = h * f(t+h/2.0, x+k2/2.0)
        k4 = h * f(t+h, x+k3)
        t = t + h
        x = x + (k1 + 2*k2 + 2*k3 + k4)/6.0
    ! 計算結果を画面に表示
    print '(4e17.9)', Exp(-t), x, 1.0d0 - x(1) - x(2)
    ! 計算結果を'rk-step3.dat'に出力
    write(10, '(4e17.9)') Exp(-t), x, 1.0d0 - x(1) - x(2)
end do
! 'rk-step3.dat'を閉じる
close(10)

! 'rk-step3.plt'を作成してgnuplotでAqauaTermによるプロットを表示させる
! 'rk-step3.plt'を開く
open(11, file='rk-step3.plt', status='replace')
    write (11, '(a)') 'set logscale x'
    write (11, '(a)') 'set format x "10^{%L}"'
    write (11, '(a)') 'set ytics 0.1'
    write (11, '(a)') 'set mytics 4'
    write (11, '(a)') 'set xrange [0.1:1000000]'
    write (11, '(a)') 'set yrange [-0.05:1.05]'
    write (11, '(a)') 'set xlabel "1+z (= Exp(-N))"'
    write (11, '(a)') 'set ylabel "Omega"'
    write (11, '(a)') 'plot "rk-step3.dat" using 1:2 w l'
    write (11, '(a)') 'replot "rk-step3.dat" using 1:3 w l'
    write (11, '(a)') 'replot "rk-step3.dat" using 1:4 w l'
    write (11, '(a)') 'set arrow from 0.1,0.68 to 1,0.68 nohead ds(10,5)'
    write (11, '(a)') 'set parametric'
    write (11, '(a)') 'set trange [-0.05:1.05]'
    write (11, '(a)') 'c1 = 1.0'
    write (11, '(a)') 'c2 = 3001.0'
    write (11, '(a)') 'replot c1, t'
    write (11, '(a)') 'replot c2, t'
    ! write (11, '(a)') 'unset parametric'
! 'rk-step3.plt'を閉じる
close(11)

! gnuplotを起動し，'rk-step1.plt'を実行する
call execute_command_line('gnuplot "rk-step3.plt"')


! 'rk-step3-save.plt'を作成してgnuplotでpng画像としてプロットを保存する
! 'rk-step3-save.plt'を開く
open(11, file='rk-step3-save.plt', status='replace')
    ! 出力先をpngに設定
    write (11, '(a)') 'set terminal postscript enhanced color'
    ! 出力ファイルを'rk-step1.png'に設定
    write (11, '(a)') 'set output "rk-step3.eps'
    ! 横軸を対数軸に設定
    write (11, '(a)') 'set logscale x'
    write (11, '(a)') 'set format x "10^{%L}"'
    write (11, '(a)') 'set ytics 0.1'
    write (11, '(a)') 'set mytics 4'
    write (11, '(a)') 'set xrange [0.1:1000000]'
    write (11, '(a)') 'set yrange [-0.05:1.05]'
    write (11, '(a)') 'set xlabel "1+z (= Exp(-N))"'
    write (11, '(a)') 'set ylabel "Omega"'
    ! 'rk-step3.dat'をプロット
    write (11, '(a)') 'plot "rk-step3.dat" using 1:2 w l'
    ! グラフを重ねるため再度，出力ファイルを'rk-step3.png'に設定
    write (11, '(a)') 'set output "rk-step3.eps"'
    ! 1列目(1+z)と3列目(x(2))を使ってグラフを描く
    write (11, '(a)') 'replot "rk-step3.dat" using 1:3 w l'
    ! グラフを重ねるため再度，出力ファイルを'rk-step3.png'に設定
    write (11, '(a)') 'set output "rk-step3.eps"'
    ! 1列目(1+z)と4列目(Omega_m)を使ってグラフを描く
    write (11, '(a)') 'replot "rk-step3.dat" using 1:4 w l'
    write (11, '(a)') 'set arrow from 0.1,0.68 to 1,0.68 nohead'
    write (11, '(a)') 'set parametric'
    write (11, '(a)') 'set trange [-0.05:1.05]'
    write (11, '(a)') 'c1 = 1'
    write (11, '(a)') 'c2 = 3001'
    ! グラフを重ねるため再度，出力ファイルを'rk-step3.png'に設定
    write (11, '(a)') 'set output "rk-step3.eps"'
    write (11, '(a)') 'replot c1, t'
    write (11, '(a)') 'set output "rk-step3.eps"'
    write (11, '(a)') 'replot c2, t'
! 'rk-step3-1save.plt'を閉じる
close(11)

! gnuplotを起動し，'rk-step3-save.plt'を実行する
call execute_command_line('gnuplot "rk-step3-save.plt"')

stop

! 以下，内部副プログラム
contains
    ! 1解微分の右辺を2成分配列としての関数f(t,x(1),x(2)) = { f1(t,x(1),x(2)), f2(t,x(1),x(2)) }として定義する
    function f(t, x)
        implicit none
        double precision f(2)
        double precision t, x(2)
            f(1) = x(1) * (x(1) - 3*x(2) - 1 )
            f(2) = x(2) * (x(1) - 3*x(2) + 3)
        return
    end function f


end program main