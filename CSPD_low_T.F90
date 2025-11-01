program CSPD_low_T
	implicit real(8) (a-h, o-z)


	allocatable :: S(:), omega(:)
	integer :: i, k, nsteps, n
	real(8) :: omega_0, eta, temp, beta, freq_max, domega, pi

	omega_0 = 3.289d0
	eta = 0.02d0
	n = 128
	temp = 0.325d0
	beta = 1.0d0/temp 
	freq_max = 5.5d0
	nsteps = 1000000
	pi = acos(-1.0d0)
	domega = freq_max/real(nsteps)

	allocate (S(nsteps), omega(nsteps))

	S = 0.0d0
	do i = 1, nsteps
   	omega(i) = (i-1)*domega
   	S(i) = 0.0d0
   	do k = 1, n
      	wk = sqrt(omega_0**2 - 2.0d0*cos(2.0d0*pi*real(k,8)/real(n,8)))
      	bk = (exp(beta*wk)-1)/wk
      	S(i) = S(i) + eta / ( eta**2 + (omega(i) - wk)**2 )
   	end do
	end do


	S = S/(real(n)*beta)
	S = S / (pi*temp)


	call write_two_col('Analytical_g43_T0325.dat', omega, S, nsteps)


contains

subroutine write_two_col(fname, x, y, n)
      character(*), intent(in) :: fname
      real(8),      intent(in) :: x(:), y(:)
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


end program CSPD_low_T
