program Effective_potential_v1
	use iso_fortran_env
	implicit none

	integer,  parameter :: dp = real64
	real(dp), parameter :: temp = 0.325_dp
	real(dp), parameter :: g = 4.3_dp
	real(dp), parameter :: beta = 1.0_dp/temp
	integer,  parameter :: steps = 2001
	real(dp), parameter :: xinter = 0.001_dp
	real(dp) :: x0, x(steps), omega(steps), a2(steps), w(steps), truepot(steps)
	integer :: i

	do i = 1, steps
		x(i) = real(i,dp)*xinter - real(steps-1,dp)*xinter/2.0_dp - xinter
		x0 = x(i)
		omega(i) = bisection(f2, 0.0001_dp, 100.0_dp, 1.0e-12_dp)	! avoid trivial root at omega=0
		a2(i) = a2f(omega(i), x0)
		w(i) = wf(omega(i), a2(i), x0)
		truepot(i) = fex(x0)
	end do

	call write_two_col('g43_T0325_effective_pot.dat', x, w, steps)
	call write_two_col('g43_Tinfty_true_pot.dat', x, truepot, steps)

contains

	function bisection(f0, a_in, b_in, err) result(root)
		real(dp), intent(in) :: a_in, b_in, err
		real(dp) :: a, b, c
		real(dp) :: root
		integer :: i
		integer, parameter :: max_iter = 100000

		interface
			function f0(x) result(y)
				use iso_fortran_env
				integer, parameter :: dp = real64
				real(dp), intent(in) :: x
				real(dp) :: y
			end function f0
		end interface

		a = a_in
		b = b_in

		do i = 1, max_iter
			c = (a+b)/2.0_dp	
			if (abs(b - a) < err) exit
			
			if (f0(c) < 0.0_dp) then
				a = c 
			else
				b = c
			end if
		end do

		root = c
	end function bisection

	function f1(x, y) result(z)
		real(dp), intent(in) :: x, y
		real(dp) :: z
		real(dp) :: coth_term
		
		if (abs(x) < 1.0e-12_dp) then	! coth x goes as 1/x as x -> 0, prevent dividing by 0
			coth_term = 1.0_dp
		else
			coth_term = (beta*x/2.0_dp)*(1.0_dp/tanh(beta*x/2.0_dp))
		end if
		
		z = x**4 - (1.69_dp+12.0_dp*g*y**2)*x**2 - (12.0_dp*g/beta)*(coth_term-1.0_dp)
	end function f1

	function f2(x) result(y)
		real(dp), intent(in) :: x 
		real(dp) :: y
		y = f1(x, x0)
	end function f2

	function a2f(omega,x) result(z)
		real(dp), intent(in) :: omega, x
		real(dp) :: z
		z = (1.0_dp/(12.0_dp*g))*(omega**2-1.69_dp)-x**2
	end function a2f

	function wf(omega, a2, x) result(y)
		real(dp), intent(in) :: omega, a2, x
		real(dp) :: y
		y = temp * log(sinh(beta*omega/2.0_dp)/(beta*omega/2.0_dp)) - omega**2*a2/2.0_dp +  &
		    0.5_dp*1.69_dp*a2 + 3.0_dp*g*a2**2 + 0.5_dp*(1.69_dp+12.0_dp*g*a2)*x**2 + g*x**4
	end function wf

	function fex(x)
		real(dp), intent(in) :: x 
		real(dp) :: fex
		fex = 0.845_dp*x**2 + g*x**4
	end function fex

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

end program Effective_potential_v1

