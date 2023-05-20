program persamaan_diferensial_tugas_08
  implicit none
  integer, parameter :: dpr = kind(1.0D0)
  integer :: i, k, n
  real(dpr), allocatable :: X_analitik(:), Y_analitik(:), X_numerik(:), Y_numerik(:)
  real(dpr), allocatable :: sol_numerik(:,:)
  real(dpr) :: h, x0, xf, t0, y0, z0, m, delta_x

  !membaca input besar step h dan jumlah output
  write(*, *) 'masukkan jumlah titik data output yang diinginkan'
  write(*, *) 'jumlah titik data output ='
  read(*, *) k
  write(*, *) 'masukkan jarak h antar titik data output (x) yang diinginkan'
  write(*, *) 'jarak antar titik data (x) ='
  read(*, *) h

  !men-generate data x_analitik untuk plotting
  x0 = 0.0
  xf = h*k
  n = 10000
  m = n
  delta_x = (xf-x0)/(m-1)

  allocate(X_analitik(n), Y_analitik(n))

  do i = 1, n
     X_analitik(i) = x0 + (i - 1) * delta_x
  end do

  !men-generate solusi analitik
  Y_analitik = solusi_analitik(X_analitik, n)

  !men-generate solusi numerik

  t0 = 0.000001
  y0 = t0
  z0 = 1.0

  allocate(sol_numerik(k,2),X_numerik(k),Y_numerik(k))

  ! Call the function to solve the ODE
  sol_numerik = runge_kutta_4_solve_2nd_ode(t0, k, h, y0, z0)
  
  do i = 1,k 
     X_numerik(i) = sol_numerik(i,1)
     Y_numerik(i) = sol_numerik(i,2)
  end do

  !menyimpan file output
  open(unit=2, file='output_xa_pd2_tugas_08.txt', status='replace', action='write')
  do i = 1, n
    write(2, *) X_analitik(i)
  end do
  close(2)

  open(unit=3, file='output_ya_pd2_tugas_08.txt', status='replace', action='write')
  do i = 1, n
    write(3, *) Y_analitik(i)
  end do
  close(3)

  open(unit=4, file='output_xn_pd2_tugas_08.txt', status='replace', action='write')
  do i = 1, k
    write(4, *) X_numerik(i)
  end do
  close(4)

  open(unit=5, file='output_yn_pd2_tugas_08.txt', status='replace', action='write')
  do i = 1, k
    write(5, *) Y_numerik(i)
  end do
  close(5)

  deallocate(X_analitik, Y_analitik, X_numerik, Y_numerik)
  stop

contains

  function solusi_analitik(X, k) result(Y)
    implicit none
    integer, intent(in) :: k
    real(dpr), intent(in) :: X(k)
    real(dpr) :: Y(k)
    integer :: i

    !fungsi solusi analitik

    do i = 1, k
       Y(i) = fungsi(X(i))
    end do

  end function solusi_analitik

  function fungsi(x) result(y)
    implicit none
    real(dpr), intent(in) :: x
    real(dpr) :: y
    y = x * exp(-x)
  end function fungsi

  function runge_kutta_4_solve_2nd_ode(t0, k, h, y0, z0) result(t_y)
    implicit none
    real(dpr), intent(in) :: t0, h, y0, z0
    integer, intent(in) :: k
    integer :: i
    real(dpr) :: y(k), z(k), t(k), t_y(k,2)
    real(dpr) :: k1y, k2y, k3y, k4y
    real(dpr) :: k1z, k2z, k3z, k4z

    ! Set initial conditions

    do i = 1,k 
       y(i) = 0.0
       z(i) = 0.0
       t(i) = 0.0
    end do

    y(1) = y0
    z(1) = z0
    t(1) = t0

    ! Perform the iterative computation
    do i = 1, k-1
      ! Compute the intermediate values
      k1y = h * z(i)
      k1z = h * func(t(i), y(i), z(i))

      k2y = h * (z(i) + k1z / 2.0)
      k2z = h * func(t(i) + h / 2.0, y(i) + k1y / 2.0, z(i) + k1z / 2.0)

      k3y = h * (z(i) + k2z / 2.0)
      k3z = h * func(t(i) + h / 2.0, y(i) + k2y / 2.0, z(i) + k2z / 2.0)

      k4y = h * (z(i) + k3z)
      k4z = h * func(t(i) + h, y(i) + k3y, z(i) + k3z)

      ! Update the values
      y(i+1) = y(i) + (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0
      z(i+1) = z(i) + (k1z + 2.0 * k2z + 2.0 * k3z + k4z) / 6.0
      t(i+1) = t(i) + h
    end do

    ! simpan output ke t_y
    do i = 1,k 
       t_y(i,1) = t(i)
       t_y(i,2) = y(i)
    end do
    
  end function runge_kutta_4_solve_2nd_ode

  function func(t, y, z) result(d2y_dt2)
    implicit none
    real(dpr), intent(in) :: t, y, z
    real(dpr) :: d2y_dt2

    ! Define the derivative function
    d2y_dt2 = -y/t -z

  end function func

end program persamaan_diferensial_tugas_08
