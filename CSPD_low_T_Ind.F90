program CSPD_low_T_Ind
	implicit real(8) (a-h, o-z)


	allocatable :: S(:), omega(:)
	integer :: i, k, nsteps, n
	real(8) :: omega_0, eta, temp, beta, freq_max, domega, pi

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
      S(i) = eta / ( eta**2 + (omega(i) - 3.58016d0)**2 )
	end do


	S = S/beta
	S = S / (pi*temp)


	call write_two_col('Analytical_kpi_g43_T0325.dat', omega, S, nsteps)


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


end program CSPD_low_T_Ind
