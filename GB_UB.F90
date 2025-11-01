program GB_UB
   use iso_fortran_env
   implicit none

   integer, parameter :: dp = real64
   real(dp), parameter :: Omega_0 = sqrt(3.69_dp)
   real(dp), parameter :: g = 0.25_dp
   real(dp), parameter :: temp = 1.3_dp
   real(dp), parameter :: beta = 1.0_dp/temp
   real(dp) :: z

   z = bisect(f, 0.00001_dp, 10.0_dp, 1e-10_dp)
   print *, z


contains

   function bisect(f0, a_in, b_in, err) result(root)
      real(dp), intent(in) :: a_in, b_in, err
      real(dp) :: a, b, c
      real(dp) :: root
      integer :: i
      integer, parameter :: max_iter = 100000000

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
   end function bisect

   function f(x) result(y)
      real(dp), intent(in) :: x 
      real(dp) :: y
      y =  0.5_dp / tanh(beta*x/2.0_dp)                                                                 &  
    - 0.25_dp*(Omega_0**2/x**2 + 1.0_dp) * ( exp(-beta*x/2.0_dp)/sinh(beta*x/2.0_dp) + 1.0_dp )    &  
    - 0.125_dp*beta*(Omega_0**2/x - x) * ( exp(-beta*x/2.0_dp)/sinh(beta*x/2.0_dp)         &
                                     * ( 1.0_dp + 1.0_dp/tanh(beta*x/2.0_dp) ) )                    & 
    - 0.5_dp*g/x**3 * ( 6.0_dp*exp(-beta*x)/tanh(beta*x/2.0_dp)                                     &
                      + 3.0_dp*exp(-beta*x/2.0_dp)/sinh(beta*x/2.0_dp) + 3.0_dp )                   & 
    - g*beta/(4.0_dp*x**2) * ( 6.0_dp*exp(-beta*x)/tanh(beta*x/2.0_dp)                              &
                             + 3.0_dp*exp(-beta*x)/sinh(beta*x/2.0_dp)**2                            &
                             + 1.5_dp*exp(-beta*x/2.0_dp)/sinh(beta*x/2.0_dp)                       &
                               * ( 1.0_dp + 1.0_dp/tanh(beta*x/2.0_dp) ) )
     
   end function f

end program GB_UB
