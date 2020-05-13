module fft_solver

  use grid_2d
  include 'fftw3.f'

  private
  public :: init_fft_solver, solve_poisson_fft, destroy_fft_solver
  
  integer(8) :: pf, pb
  real, dimension(:), allocatable :: a, b, c, c1, d1, mwn
  double complex, allocatable :: in(:), out(:)

contains

  subroutine init_fft_solver(left,right,top,bottom)

    implicit none

    ! Type of BC
    character(len=*), intent(in) :: left, right, top, bottom

    ! Local variables
    integer :: k, j

    allocate(in(nx))
    allocate(out(nx))
    if (left /= 'periodic') then
      print *, 'No periodic direction, unable to use FFT solver'
      stop
    endif

    ! Tridiagonal solver coefficient
    allocate(a(ny))
    allocate(b(ny))
    allocate(c(ny))
    allocate(c1(ny))
    allocate(d1(ny))
    do j = 1,ny
      a(j) = 1.0/dy**2 
      b(j) = -2.0/dy**2
      c(j) = 1.0/dy**2
    end do

    c1(1) = c(1)/b(1)
    do j = 2,ny
      c1(j) = c(j)/(b(j) - a(j)*c1(j-1))
    end do

    ! Modified wave number
    allocate(mwn(nx))
    do k = 1,nx
      mwn(k) = 2.0*(cos(2.0*pi*k/nx) - 1.0)/dx**2
    end do 

    ! Create the plan
    call dfftw_plan_dft_1d(pf,nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(pb,nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)

  end subroutine init_fft_solver

  subroutine solve_poisson_fft(rhs)

    ! This subroutine solve the poisson equation
    ! nabal^2 f = rhs
    ! with right-hand side rhs. At the end of the dubroutine 
    ! rhs contains the solution f.

    use boundary_conditions

    implicit none

    real, dimension(:,:), intent(inout) :: rhs

    integer :: j, k
    real, dimension(nx,ny) :: f

    ! For every j compute the FFT of the RHS in the x direction
    ! and store it in f
    do j = 1,ny
      in = rhs(:,j)
      call dfftw_execute_dft(pf, in, out)
      f(:,j) = sqrt(real(out)**2 + imag(out)**2) / ny
    end do

    ! Solve tridiagonal system for every k
    do k = 1,nx 

      d1(1) = f(k,1)/(b(1) + mwn(k))
      do j = 2,ny
        d1(j) = (f(k,j) - a(j)*d1(j-1))/(b(j) - a(j)*c1(j))
      end do

      f(k,ny) = d1(ny)
      do j = ny-1,1,-1
        f(k,j) = d1(j) - c1(j)*f(k,j+1)
      end do
    end do

    ! Compute inverse transform of f
    do j = 1,ny
      in = f(:,j)
      call dfftw_execute_dft(pb, in, out)
      do k = 1,nx
        rhs(k,j) = sqrt(real(out(k))**2 + imag(out(k))**2) / ny
      end do
    end do 

  end subroutine solve_poisson_fft

  subroutine destroy_fft_solver

    implicit none

  end subroutine destroy_fft_solver

end module fft_solver
