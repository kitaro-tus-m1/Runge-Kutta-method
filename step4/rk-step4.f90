! \Lambda-CDM modelにおける密度パラメータの時間変化をRunge-Kutta methodで求める
program main

! runge_kutta という subroutine を定義してリファクタリングする

    !-------------------------------------------------------------
    !----- 1. 変数定数 定義 ----------------------------------------
    !-------------------------------------------------------------
        implicit none
        ! 表記を簡単にするため，x := \Omega_r, y := \Omega_{\lambda}, t := N とする
        ! 解く方程式は x'(t) = x (x - 3y - 1) =: f(x, y), y'(t) = y (x - 3y + 3) =: g(x, y)
        ! 初期値は x(0) = unknown, y(0) = 0.68
        ! 境界条件は y(-Ln(3001)) = 1 - 2 x(-Ln(3001))
        ! x(0)がunknownなので，x(0)を色々と仮定して解いてみる
        ! その試行の中で，境界条件を満たすものが現れたら，それが解であるとみなす
    
        ! t0 := t=0, x0 := x(0), y0 := y(0)
        ! 配列を使って，ベクトル表現で解くためには
        ! x(2) := (x, y), f(2) := (f, g), x0(2) := (x0, y0)とおくべき
        ! 1+z = 0.1 から計算を始めるためにt0は次のようにした
        double precision,parameter :: t0 = -Log(0.10d0), x0 = 0.000000015d0, y0 = 0.99955d0
        ! x0, y0 を色々と調節した
        ! この調節をできれば自動化したかったが，いい方法が思いつかなかった
    
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
    
    
    !-------------------------------------------------------------
    !----- 2. Runge-Kutta method ---------------------------------
    !-------------------------------------------------------------
    
    open(10, file='rk-step4.dat', status='replace')
        ! 初期値を画面に表示
        print *, Exp(-t), x, 1.0d0 - x(1) - x(2)
        ! 初期値を'rk-step4.dat'に出力
        write(10, '(3e17.9)') Exp(-t), x, 1.0d0 - x(1) - x(2)

        ! Runge-Kutta の繰り返し
        do i = 1, 180
            call runge_kutta(t, x, h)
            ! 計算結果を画面に表示
            print '(4e17.9)', Exp(-t), x, 1.0d0 - x(1) - x(2)
            ! 計算結果を'rk-step4.dat'に出力
            write(10, '(4e17.9)') Exp(-t), x, 1.0d0 - x(1) - x(2)
        end do
    ! 'rk-step4.dat'を閉じる
    close(10)
    
    !-------------------------------------------------------------
    !----- 3. 結果を画面に表示 --------------------------------------
    !-------------------------------------------------------------
    
    ! 'rk-step4.plt'を作成してgnuplotでAqauaTermによるプロットを画面に表示させる
    ! 'rk-step4.plt'を開く
    open(11, file='rk-step4.plt', status='replace')
        ! 内部subroutine 'set_gnuplot_options' を呼び出しgnuplotの設定をする
        call set_gnuplot_options
        ! 内部subroutine 'gnuplot_plot_screen' を呼び出し，グラフを画面に表示する
        call gnuplot_plot_screen
    ! 'rk-step4.plt'を閉じる
    close(11)
    
    ! gnuplotを起動し，'rk-step4.plt'を実行する
    call execute_command_line('gnuplot "rk-step4.plt"')
    
    
    !-------------------------------------------------------------
    !----- 4. 結果をeps画像として保存 -------------------------------
    !-------------------------------------------------------------
    
    ! 'rk-step4-save.plt'を作成してgnuplotでeps画像としてプロットを保存する
    ! 'rk-step4-save.plt'を開く
    open(11, file='rk-step4-save.plt', status='replace')
        ! 出力先をeps(color)に設定
        write (11, '(a)') 'set terminal postscript enhanced color'
        ! 出力ファイルを'rk-step4.eps'に設定
        write (11, '(a)') 'set output "rk-step4.eps'
    
        ! gnuplotの設定
        call set_gnuplot_options
        ! 内部subroutine 'gnuplot_plot_output' を呼び出し，グラフをepsに出力保存する
        call gnuplot_plot_eps
    
    ! 'rk-step4-1save.plt'を閉じる
    close(11)
    
    ! gnuplotを起動し，'rk-step4-save.plt'を実行する
    call execute_command_line('gnuplot "rk-step4-save.plt"')
    
    stop
    
    !-------------------------------------------------------------
    !----- 5. 内部副プログラム --------------------------------------
    !-------------------------------------------------------------
    
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
    
        ! Runge-Kutta method を行う内部副サブルーチンを定義する
        subroutine runge_kutta(t, x, h)
            ! Runge-Kutta method の各kを2要素配列として定義する
            double precision k1(2), k2(2), k3(2), k4(2), t, x(2), h
            k1 = h * f(t, x)
            k2 = h * f(t+h/2.0, x+k1/2.0)
            k3 = h * f(t+h/2.0, x+k2/2.0)
            k4 = h * f(t+h, x+k3)
            t = t + h
            x = x + (k1 + 2*k2 + 2*k3 + k4)/6.0
        end subroutine

        ! gnuplotの設定を行う内部subroutine
        subroutine set_gnuplot_options
            write (11, '(a)') 'set logscale x'
            write (11, '(a)') 'set format x "10^{%L}"'
            write (11, '(a)') 'set ytics 0.1'
            write (11, '(a)') 'set mytics 5'
            write (11, '(a)') 'set xrange [0.1:1000000]'
            write (11, '(a)') 'set yrange [-0.05:1.05]'
            write (11, '(a)') 'set xlabel "1+z  ( = e^{-N} )"'
            write (11, '(a)') 'set ylabel "{/Symbol W}"'
            write (11, '(a)') 'set key outside invert reverse box Left width 2'
            write (11, '(a)') 'set size square'
        end subroutine
    
        ! plotを画面に表示するための内部subroutine
        subroutine gnuplot_plot_screen
            write (11, '(a)') 'set parametric'
            write (11, '(a)') 'set trange [-0.05:1.05]'
            write (11, '(a)') 'c1 = 1.0'
            write (11, '(a)') 'c2 = 3001.0'
            write (11, '(a)') 'c3 = 0.68'
            write (11, '(a)') 'plot c1, t lc "dark-gray" notitle'
            write (11, '(a)') 'replot c2, t lc "dark-gray" notitle'
            write (11, '(a)') 'replot t, c3 lc "dark-gray" notitle'
            write (11, '(a)') 'replot "rk-step4.dat" using 1:2 w l lt rgb "blue" title "{/Symbol W}_r"'
            write (11, '(a)') 'replot "rk-step4.dat" using 1:4 w l lt rgb "dark-green" title "{/Symbol W}_m"'
            write (11, '(a)') 'replot "rk-step4.dat" using 1:3 w l lt rgb "red" title "{/Symbol W}_{/Symbol L}"'
        end subroutine
    
        ! plotをepsに保存するための内部subroutine
        subroutine gnuplot_plot_eps
            write (11, '(a)') 'set parametric'
            write (11, '(a)') 'set trange [-0.05:1.05]'
            write (11, '(a)') 'c1 = 1.0'
            write (11, '(a)') 'c2 = 3001.0'
            write (11, '(a)') 'c3 = 0.68'
            write (11, '(a)') 'plot c1, t lc "dark-gray" notitle'
            ! グラフを重ねるため再度，出力ファイルを'rk-step4.eps'に設定
            write (11, '(a)') 'set output "rk-step4.eps"'
            write (11, '(a)') 'replot c2, t lc "dark-gray" notitle'
            ! グラフを重ねるため再度，出力ファイルを'rk-step4.eps'に設定
            write (11, '(a)') 'set output "rk-step4.eps"'
            write (11, '(a)') 'replot t, c3 lc "dark-gray" notitle'
            ! グラフを重ねるため再度，出力ファイルを'rk-step4.eps'に設定
            write (11, '(a)') 'set output "rk-step4.eps"'
            write (11, '(a)') 'replot "rk-step4.dat" using 1:2 w l lt rgb "blue" lw 2 title "{/Symbol W}_r"'
            ! グラフを重ねるため再度，出力ファイルを'rk-step4.eps'に設定
            write (11, '(a)') 'set output "rk-step4.eps"'
            write (11, '(a)') 'replot "rk-step4.dat" using 1:4 w l lt rgb "dark-green" lw 2 title "{/Symbol W}_m"'
            ! グラフを重ねるため再度，出力ファイルを'rk-step4.eps'に設定
            write (11, '(a)') 'set output "rk-step4.eps"'
            write (11, '(a)') 'replot "rk-step4.dat" using 1:3 w l lt rgb "red" lw 2 title "{/Symbol W}_{/Symbol L}"'
        end subroutine
    end program main