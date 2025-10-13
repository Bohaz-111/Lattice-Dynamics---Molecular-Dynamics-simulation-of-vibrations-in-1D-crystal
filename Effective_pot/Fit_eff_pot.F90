program Fit_eff_pot
    use iso_fortran_env
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: max_points = 10000
    integer, parameter :: degree = 6  ! Maximum even degree 
    integer, parameter :: steps = 4001
    real(dp), parameter :: xinter = 0.01
    real(dp), parameter :: g = 4.3_dp  ! Anharmonic strength
    
    real(dp), allocatable :: x(:), w(:)
    real(dp), allocatable :: A(:,:), b(:), coeffs(:)
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: x_fit(:), w_fit(:)
    
    integer :: n_points, n_coeffs
    integer :: i, j, info, lwork
    real(dp) :: force(steps), tforce(steps), kf(steps), tfk(steps)
    real(dp) :: xax(steps)
    character(len=256) :: input_file
    character(len=256) :: output_coeff_file, output_fit_file, force_file, tforce_file, kf_file, tkf_file
    
    input_file = 'g43_T0325_effective_pot.dat'
    output_coeff_file = 'g43_T0325_fitted_coefficients.dat'
    output_fit_file = 'g43_T0325_fitted_potential.dat'
    force_file = 'g43_T0325_force.dat'
    kf_file = 'g43_T_0325kf.dat'
    tforce_file = 'g43_tforce.dat'
    tkf_file = 'g43_kf.dat'
    n_coeffs = degree/2 + 1
    
    call read_data(input_file, x, w, n_points)
    
    allocate(A(n_points, n_coeffs))
    allocate(b(n_points))
    allocate(coeffs(n_coeffs))
    
    do i = 1, n_points
        do j = 1, n_coeffs
            A(i, j) = x(i)**(2*(j-1))
        end do
    end do
    
    b = w
    
    ! Solve least squares problem using LAPACK dgels
    ! Query for optimal workspace size
    lwork = -1
    allocate(work(1))
    call dgels('N', n_points, n_coeffs, 1, A, n_points, b, n_points, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    
    call dgels('N', n_points, n_coeffs, 1, A, n_points, b, n_points, work, lwork, info)
    
    if (info /= 0) then
        print *, 'ERROR: LAPACK dgels failed with info = ', info
        stop
    end if
    
    coeffs = b(1:n_coeffs)
    
    print *, 'Fitted coefficients:'
    do i = 1, n_coeffs
        print '(A,I2,A,ES16.8)', 'a', 2*(i-1), ' = ', coeffs(i)
    end do
    print *, ''
    
    allocate(x_fit(n_points))
    allocate(w_fit(n_points))
    x_fit = x
    
    do i = 1, n_points
        w_fit(i) = 0.0_dp
        do j = 1, n_coeffs
            w_fit(i) = w_fit(i) + coeffs(j) * x(i)**(2*(j-1))
        end do
    end do

    do i = 1, steps
        xax(i) = real(i,dp)*xinter - real(steps-1,dp)*xinter/2.0_dp - xinter
        force(i) = -(2.0_dp*coeffs(2)*xax(i) + 4.0_dp*coeffs(3)*xax(i)**3 + 6.0_dp*coeffs(4)*xax(i)**5 + 8.0_dp*coeffs(5)*xax(i)**7) 
        tforce(i) = -(1.69_dp*xax(i) + 4.0_dp*g*xax(i)**3)
        kf(i) = sqrt(2.0_dp*coeffs(2)+12.0_dp*coeffs(3)*abs(xax(i))**2 + 30.0_dp*coeffs(4)*abs(xax(i))**4 + 56.0_dp*coeffs(5)*abs(xax(i))**6)
        tfk(i) = sqrt(1.69_dp + 12.0_dp*g*(xax(i)**2))
    end do
    
    call calculate_goodness_of_fit(w, w_fit, n_points)
    call write_coefficients(output_coeff_file, coeffs, n_coeffs, degree)
    call write_two_col(output_fit_file, x_fit, w_fit, n_points)
    call write_two_col(force_file, xax, force, steps)
    call write_two_col(tforce_file, xax, tforce, steps)
    call write_two_col(kf_file, xax, kf, steps)
    
    deallocate(x, w, A, b, coeffs, work, x_fit, w_fit)



contains

    subroutine read_data(fname, x_out, y_out, n)
        character(*), intent(in) :: fname
        real(dp), allocatable, intent(out) :: x_out(:), y_out(:)
        integer, intent(out) :: n
        
        real(dp) :: temp_x(max_points), temp_y(max_points)
        integer :: u, ios, count
        character(len=256) :: line
        
        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR: Cannot open file ', trim(fname)
            stop
        end if
        
        count = 0
        do
            read(u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            ! Skip comment lines
            if (line(1:1) == '#') cycle
            
            count = count + 1
            if (count > max_points) then
                print *, 'ERROR: Too many data points'
                stop
            end if
            
            read(line, *, iostat=ios) temp_x(count), temp_y(count)
            if (ios /= 0) then
                print *, 'ERROR: Cannot read line: ', trim(line)
                stop
            end if
        end do
        
        close(u)
        
        n = count
        allocate(x_out(n), y_out(n))
        x_out = temp_x(1:n)
        y_out = temp_y(1:n)
    end subroutine read_data

    subroutine calculate_goodness_of_fit(y_obs, y_fit, n)
        real(dp), intent(in) :: y_obs(:), y_fit(:)
        integer, intent(in) :: n
        
        real(dp) :: ss_res, ss_tot, r_squared
        real(dp) :: y_mean
        integer :: i

        y_mean = sum(y_obs) / real(n, dp)
        ss_res = 0.0_dp
        ss_tot = 0.0_dp

        do i = 1, n
            ss_res = ss_res + (y_obs(i) - y_fit(i))**2
            ss_tot = ss_tot + (y_obs(i) - y_mean)**2
        end do

        r_squared = 1.0_dp - ss_res / ss_tot

        print '(A,F10.6)', 'R² = ', r_squared
    end subroutine calculate_goodness_of_fit

    subroutine write_coefficients(fname, coeffs, n, degree)
        character(*), intent(in) :: fname
        real(dp), intent(in) :: coeffs(:)
        integer, intent(in) :: n, degree
        integer :: u, i
        
        open(newunit=u, file=fname, status='replace', action='write')
        write(u,'(A)') '# Even Polynomial Coefficients'
        write(u,'(A,I2)') '# Maximum degree: ', degree
        write(u,'(A)') '# Power    Coefficient'
        do i = 1, n
            write(u,'(I6,4X,ES20.12)') 2*(i-1), coeffs(i)
        end do
        close(u)
        print *, 'Coefficients saved to: ', trim(fname)
    end subroutine write_coefficients

    subroutine write_two_col(fname, x, y, n)
        character(*), intent(in) :: fname
        real(dp),     intent(in) :: x(:), y(:)
        integer,      intent(in) :: n
        integer :: u,i
        open(newunit=u, file=fname, status='replace', action='write')
            write(u,'(A)') '#x        y'
            do i=1,n
                write(u,'(2ES16.8)') x(i), y(i)
            end do
        close(u)
        print *, 'saved to:  ', fname
    end subroutine write_two_col


end program Fit_eff_pot.F90



