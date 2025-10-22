program CSD_ML_v1
   use iso_fortran_env
   use omp_lib
   implicit none

   integer,  parameter :: dp    = real64
   real(dp), parameter :: pi    = 3.14159265358979323846264338327950288419_dp
   real(dp), parameter :: temp  = 0.325_dp
   real(dp), parameter :: beta  = 1.0_dp/temp
   real(dp), parameter :: dt    = 0.1_dp     
   integer,  parameter :: steps = 4000                
   integer,  parameter :: freq_steps = 5500
   integer,  parameter :: ntraj = 10000
   real(dp), parameter :: eta = 0.02_dp ! Artificial damping
   real(dp), parameter :: omega_0 = sqrt(2.0_dp*2.4040824893703618_dp)
   real(dp), parameter :: beta_eff = (exp(beta*omega_0)-1.0_dp)/omega_0

   real(dp)    :: q, p, f
   real(dp)    :: vacf(steps), vacf_sum(steps), time_lag(steps), vdos(freq_steps), omega(freq_steps)
   real(dp),    allocatable :: p_series(:)
   integer     :: i, j, k, traj

   call random_seed()

   vacf = 0.0_dp
   vacf_sum = 0.0_dp

   do i = 1, steps
      time_lag(i) = (i-1)*dt
   end do


   do traj = 1, ntraj

      allocate(p_series(steps))

      call init_config(q, p)
      p_series(1) = p
      vacf(1) = p_series(1)*p_series(1)

      do i = 2, steps
         call step_vv(q, p, f, dt)
         p_series(i) = p
         vacf(i) = p_series(i)*p_series(1) * exp(-eta*real(i-1,dp)*dt)
      end do

      vacf_sum = vacf_sum + vacf
      deallocate(p_series)
   end do

   vacf = vacf_sum * beta_eff/ (real(ntraj, dp)*beta)
   call write_two_col('CSPD_g43_T13_VACF.dat', time_lag, vacf, steps)
   call compute_vdos(vacf, vdos, omega, steps, dt)
   call write_two_col('CSPD_g43_T13_VDOS.dat', omega, vdos, freq_steps)


contains

   subroutine init_config(q, p)
      real(dp), intent(out) :: q, p
      p = sqrt(1.0_dp/beta_eff) * randn()
      q = sqrt(1.0_dp/(beta_eff*omega_0**2)) * randn()
   end subroutine init_config

   pure subroutine force(q, f)
      real(dp), intent(in)  :: q
      real(dp), intent(out) :: f
      real(dp), parameter :: c2 = 2.4040824893703618_dp
      real(dp), parameter :: c4 = 4.0179753918609205_dp
      real(dp), parameter :: c6 = 0.046363877903540011_dp

      f = -2.0_dp*c2*q - 4.0_dp*c4*q**3 - 6.0_dp*c6*q**5

   end subroutine force

   pure subroutine step_vv(q, p, f, dt)
      real(dp), intent(inout) :: q, p, f
      real(dp), intent(in)    :: dt
      real(dp) :: halfdt
      integer  :: i
      halfdt = 0.5_dp*dt
      call force(q, f)
      p = p + halfdt * f
      q = q + dt * p
      call force(q, f)
      p = p + halfdt * f
   end subroutine step_vv

   real(dp) function randn()
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1    = max(u1, 1.0e-12_dp)
      randn = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pi*u2)
   end function randn

   subroutine compute_vdos(vacf, vdos, omega, nvacf, dt_in)
      real(dp), intent(in) :: vacf(:)
      real(dp), intent(out) :: vdos(:), omega(:)
      integer, intent(in) :: nvacf
      real(dp), intent(in) :: dt_in
      integer :: i, j, nfreq
      real(dp), parameter :: d_omega = 0.001_dp

      nfreq = size(vdos)
      vdos = 0.0_dp
      omega = 0.0_dp

      do i = 1, nfreq
         omega(i) = real(i-1, dp) * d_omega
         vdos(i) = vdos(i) + 0.5_dp*vacf(1)  ! Trapezoid rule: half weight for first point
         do j = 2, nvacf
            vdos(i) = vdos(i) + vacf(j) * cos(omega(i) * real(j-1, dp) * dt_in)
         end do
      end do
      vdos = vdos * (2.0_dp * dt_in) / (pi*temp)

   end subroutine compute_vdos

   subroutine write_two_col(fname, x, y, n)
      character(*), intent(in) :: fname
      real(dp),     intent(in) :: x(:), y(:)
      integer,      intent(in) :: n
      integer :: u,i
      open(newunit=u, file=fname, status='replace', action='write')
         write(u,'(A)') '# Time    Value'
         do i=1,n
               write(u,*) x(i), y(i)
         end do
      close(u)
      print *, 'saved to:  ', fname
   end subroutine write_two_col

end program CSD_ML_v1



