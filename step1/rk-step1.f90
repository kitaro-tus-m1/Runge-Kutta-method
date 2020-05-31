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

    ! 初期値を画面に表示
    print *, x, y

    ! 'rk-step1.dat'を開く
    open(10, file='rk-step1.dat', status='replace')
        ! 初期値を'rk-step1.dat'に出力
        write(10, *) x, y

    do i = 1, 100
        k1 = h * f(x, y)
        k2 = h * f(x+h/2.0, y+k1/2.0)
        k3 = h * f(x+h/2.0, y+k2/2.0)
        k4 = h* f(x+h, y+k3)
        x = x + h
        y = y + (k1 + 2*k2 + 2*k3 + k4)/6.0
        ! 計算結果を画面に表示
        print *, x, y

        ! 計算結果を'rk-step1.dat'に出力
            write(10, *) x, y
    end do

    ! 'rk-step1.dat'を閉じる
    close(10)

! 'rk-step1.plt'を作成してgnuplotでAqauaTermによるプロットを表示させる
    ! 'rk-step1.plt'を開くgnuplotのスクリプトファイル
    open(11, file='rk-step1.plt', status='replace')
        write (11, '(a)') 'plot "rk-step1.dat"'
        ! Returnを押すまで待機
        !write (11, '(a)') 'pause -1'
    ! 'rk-step1.plt'を閉じる
    close(11)

    ! gnuplotを起動し，'rk-step1.plt'を実行する
    call execute_command_line('gnuplot "rk-step1.plt"')

! 'rk-step1-save.plt'を作成してgnuplotでpng画像としてプロットを保存する
    ! 'rk-step1-save.plt'を開くgnuplotのスクリプトファイル
    open(12, file='rk-step1-save.plt', status='replace')
        ! 出力先をpngに設定
        write (12, '(a)') 'set terminal png'
        ! 出力ファイルを'rk-step1.png'に設定
        write (12, '(a)') 'set output "rk-step1.png"'
        ! 'rk-step1.dat'をプロット
        write (12, '(a)') 'plot "rk-step1.dat"'
    ! 'rk-step1-save.plt'を閉じる
    close(12)

    ! gnuplotを起動し，'rk-step1-save.plt'を実行する
    call execute_command_line('gnuplot "rk-step1-save.plt"')

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