program CMD_ind_v4
   use iso_fortran_env
   use omp_lib
   implicit none

   integer,  parameter :: dp    = real64
   integer,  parameter :: n     = 128                ! Number of atoms
   real(dp), parameter :: pi    = 3.14159265358979323846264338327950288419_dp
   real(dp), parameter :: temp  = 0.325_dp
   real(dp), parameter :: g     = 4.3_dp            ! Anharmonic (quartic) strength
   real(dp), parameter :: h     = 0.845_dp           ! On-site quadratic coefficient
   real(dp), parameter :: tau   = 5.0_dp             ! Thermostat time constant
   real(dp), parameter :: dt    = 0.1_dp     
   integer,  parameter :: steps = 2400         
   integer,  parameter :: ntau  = steps/2            
   integer,  parameter :: equil_steps = 4800    
   integer,  parameter :: freq_steps = 5500
   integer,  parameter :: ntraj = 2000
   real(dp), parameter :: eta = 0.02_dp ! Artificial damping

   real(dp)    :: q(n), p(n), f(n)
   real(dp)    :: vacf_k0(ntau), vacf_kpi(ntau), vacf_k0_sum(ntau), vacf_kpi_sum(ntau)
   real(dp)    :: time_lag(ntau), vdos_k0(freq_steps), vdos_kpi(freq_steps), omega(freq_steps)
   complex(dp), allocatable :: u_k0(:), u_kpi(:) 
   real(dp),    allocatable :: p_series(:, :)
   integer     :: i, j, k, traj, k_pi_idx
   real(dp)    :: theta, sumr, sumi


   call random_seed()

   vacf_k0_sum = 0.0_dp
   vacf_kpi_sum = 0.0_dp
   k_pi_idx = n/2  ! k=pi mode index

   do traj = 1, ntraj
      print *, 'Trajectory', traj, '/', ntraj

      allocate(u_k0(steps))
      allocate(u_kpi(steps))
      allocate(p_series(steps, n))

      call init_config(q, p)
      do i = 1, equil_steps
         call step_vv(q, p, f, dt, g, h)
         call thermostat(p, dt, tau, temp)
      end do

      do i = 1, steps
         call step_vv(q, p, f, dt, g, h)
         p_series(i, :) = p(:)
      end do
     
      ! Compute Fourier transform for k=0 mode
      do i = 1, steps
         sumr = 0.0_dp
         sumi = 0.0_dp
         do j = 1, n
            sumr  = sumr + p_series(i, j)
         end do
         u_k0(i) = cmplx(sumr, sumi, kind=dp) / sqrt(real(n,dp))
      end do

      ! Compute Fourier transform for k=pi mode
      do i = 1, steps
         sumr = 0.0_dp
         sumi = 0.0_dp
         do j = 1, n
            sumr  = sumr + p_series(i, j)*(-1.0_dp)**j
         end do
         u_kpi(i) = cmplx(sumr, sumi, kind=dp) / sqrt(real(n,dp))
      end do
   
      call compute_vacf_single(u_k0, steps, time_lag, vacf_k0, ntau, eta, dt)
      vacf_k0_sum = vacf_k0_sum + vacf_k0
      call compute_vacf_single(u_kpi, steps, time_lag, vacf_kpi, ntau, eta, dt)
      vacf_kpi_sum = vacf_kpi_sum + vacf_kpi

      deallocate(u_k0)
      deallocate(u_kpi)
      deallocate(p_series)
   end do

   ! Average over trajectories
   vacf_k0 = vacf_k0_sum / real(ntraj, dp)
   vacf_kpi = vacf_kpi_sum / real(ntraj, dp)

   ! Write VACF data
   call write_two_col('AnewCMD_ind_v3_k0_VACF.dat', time_lag, vacf_k0, ntau)
   call write_two_col('AnewCMD_ind_v3_kpi_VACF.dat', time_lag, vacf_kpi, ntau)

   ! Compute and write VDOS
   call compute_vdos(vacf_k0, vdos_k0, omega, ntau, dt)
   call write_two_col('AnewCMD_ind_v3_k0_VDOS.dat', omega, vdos_k0, freq_steps)

   call compute_vdos(vacf_kpi, vdos_kpi, omega, ntau, dt)
   call write_two_col('AnewCMD_ind_v3_kpi_VDOS.dat', omega, vdos_kpi, freq_steps)


contains

   subroutine init_config(q, p)
      real(dp), intent(out) :: q(:), p(:)
      integer :: i, nn
      nn = size(q)
      do i = 1, nn
         q(i) = real(i-1, dp)
         p(i) = sqrt(temp) * randn()
      end do
   end subroutine init_config

   pure subroutine force(q, f, g, h)
      real(dp), intent(in)  :: q(:), g, h
      real(dp), intent(out) :: f(:)
      real(dp), parameter :: c2 = 4.9638734265770097_dp
      real(dp), parameter :: c4 = 2.6603906937773196_dp
      real(dp), parameter :: c6 = 0.64207861953736323_dp
      integer :: jp, jm, j, nn
      real(dp) :: delta_q(size(q)), delta_qp, delta_qm
      nn = size(q)
    
      do j = 1, nn
         delta_q(j) = q(j) - real(j-1, dp)
         if (delta_q(j) > real(nn/2, dp)) then
            delta_q(j) = delta_q(j) - real(nn, dp)
         else if (delta_q(j) < -real(nn/2, dp)) then
            delta_q(j) = delta_q(j) + real(nn, dp)
         end if
      end do
    
      do j = 1, nn
         jm = j-1; if (jm < 1)  jm = nn
         jp = j+1; if (jp > nn) jp = 1
       
         delta_qp = delta_q(jp)
         delta_qm = delta_q(jm)
       
         f(j) = -(2.0_dp*c2*delta_q(j) + 4.0_dp*c4*delta_q(j)**3 + 6.0_dp*c6*delta_q(j)**5) + delta_qp + delta_qm
      end do
   end subroutine force

   pure subroutine step_vv(q, p, f, dt, g, h)
      real(dp), intent(inout) :: q(:), p(:), f(:)
      real(dp), intent(in)    :: dt, g, h
      real(dp) :: halfdt
      integer  :: i, nn
      nn = size(q)
      halfdt = 0.5_dp*dt
      call force(q, f, g, h)
      p = p + halfdt * f
      q = q + dt * p
      do i = 1, nn
         q(i) = q(i) - real(nn,dp) * floor( q(i) / real(nn,dp) )
      end do
      call force(q, f, g, h)
      p = p + halfdt * f
   end subroutine step_vv

   subroutine thermostat(p, dt, tau, T)
      real(dp), intent(inout) :: p(:)
      real(dp), intent(in)    :: dt, tau, T

      integer  :: np, dof, kshape
      real(dp) :: K, Kbar, c, s, r1, sum_r2, alpha2, alpha, factor

      np   = size(p)
      dof  = np
      K    = 0.5_dp * sum(p*p)
      Kbar = 0.5_dp * real(dof,dp) * T

      c = exp(-dt/tau)
      s = 1.0_dp - c

      r1     = randn()
      kshape = (dof-2)/2
      sum_r2 = rand_gamma(kshape) + randn()**2

      factor = (Kbar / (real(dof,dp)*max(K, tiny(1.0_dp))))
      alpha2 = c                                            &
            + factor * s * (r1*r1 + sum_r2)               &
            + 2.0_dp * exp(-0.5_dp*dt/tau)                &
               * sqrt( factor * s ) * r1

      alpha2 = max(alpha2, 0.0_dp)
      alpha  = sqrt(alpha2)
      p      = alpha * p
   end subroutine thermostat

   real(dp) function rand_gamma(k)
      integer, intent(in) :: k
      real(dp) :: v(k)
      integer  :: i
      call random_number(v)
      do i = 1, k
         v(i) = max(v(i), 1.0e-12_dp)
      end do
      v = -log(v)
      rand_gamma = 2.0_dp * sum(v)
   end function rand_gamma

   real(dp) function randn()
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1    = max(u1, 1.0e-12_dp)
      randn = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pi*u2)
   end function randn

   subroutine compute_vacf_single(u_k, nsteps, time_lag, vacf_out, ntau_out, eta, dt)
      complex(dp), intent(in) :: u_k(:)
      integer, intent(in) :: nsteps, ntau_out
      real(dp), intent(in) :: eta, dt
      real(dp), intent(out) :: vacf_out(:), time_lag(:)
      integer :: i, j, n_origins

      vacf_out = 0.0_dp
      time_lag = 0.0_dp

      do i = 1, nsteps/2
         do j = 1, nsteps/2
            vacf_out(i) = vacf_out(i) + real(conjg(u_k(j))*u_k(j+(i-1)), dp)
         end do
         vacf_out(i) = vacf_out(i) / real(nsteps/2, dp)  
         vacf_out(i) = vacf_out(i) * exp(-eta*real(i-1,dp)*dt)
         time_lag(i) = (i-1) * dt 
      end do
   end subroutine compute_vacf_single


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
         vdos(i) = vdos(i) + 0.5_dp*vacf(1)
         do j = 2, nvacf
            vdos(i) = vdos(i) + vacf(j) * cos(omega(i) * real(j-1, dp) * dt_in)
         end do
      end do
      vdos = (vdos+temp) * (2.0_dp * dt_in) / (pi*temp)

   end subroutine compute_vdos

   subroutine write_two_col(fname, x, y, n)
      character(*), intent(in) :: fname
      real(dp),     intent(in) :: x(:), y(:)
      integer,      intent(in) :: n
      integer :: u,i
      open(newunit=u, file=fname, status='replace', action='write')
         write(u,'(A)') '# Time    Value'
         do i=1,n
               write(u,'(2ES16.8)') x(i), y(i)
         end do
      close(u)
      print *, 'saved to:  ', fname
   end subroutine write_two_col

end program CMD_ind_v4


