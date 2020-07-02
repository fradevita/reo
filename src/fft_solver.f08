module fft_solver

  use grid_2D, only : dp

  include 'fftw3.f'

  private
  public :: init_fft_solver, solve_poisson_fft, destroy_fft_solver
  
  integer(8) :: pf, pb
  real(kind=dp), dimension(:), allocatable :: a, b, c, c1, d1, mwn
  double complex, allocatable :: in(:), out(:)

contains

  !===============================================================================
  subroutine init_fft_solver()

    use grid_2D, only : boxarray, pi, dp

    implicit none

    ! Local variables
    integer :: k, j, nx, ny
    real(kind=dp) :: delta

    nx = boxarray(1)%nx
    ny = boxarray(1)%ny
    delta = boxarray(1)%delta

    allocate(in(nx))
    allocate(out(nx))

    if (boxarray(1)%left /= 'periodic') then
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
      a(j) = 1.0/delta**2 
      b(j) = -2.0/delta**2
      c(j) = 1.0/delta**2
    end do

    c1(1) = c(1)/b(1)
    do j = 2,ny
      c1(j) = c(j)/(b(j) - a(j)*c1(j-1))
    end do

    ! Modified wave number
    allocate(mwn(nx))
    do k = 1,nx
      mwn(k) = 2.0*(cos(2.0*pi*k/nx) - 1.0)/delta**2
    end do

    ! Create the plan
    call dfftw_plan_dft_1d(pf,nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(pb,nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)

  end subroutine init_fft_solver
  !===============================================================================

  !===============================================================================
  subroutine solve_poisson_fft(rhs)

    ! This subroutine solve the poisson equation
    ! nabal^2 f = rhs
    ! with right-hand side rhs. At the end of the dubroutine 
    ! rhs contains the solution f.

    use grid_2D, only : boxarray, myarray, dp
    use boundary_conditions

    implicit none

    ! Input / Output variables
    type(myarray), intent(inout) :: rhs(:)

    ! Local variables
    integer :: j, k, nx, ny
    real(kind=dp), dimension(:,:), allocatable :: f

    nx = boxarray(1)%nx
    ny = boxarray(1)%ny
  
    allocate(f(nx,ny))

    ! For every j compute the FFT of the RHS in the x direction
    ! and store it in f
    do j = 1,ny
      in = rhs(1)%f(:,j)
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
        rhs(1)%f(k,j) = sqrt(real(out(k))**2 + imag(out(k))**2) / ny
      end do
    end do 

    deallocate(f)

  end subroutine solve_poisson_fft
  !===============================================================================

  !===============================================================================
  subroutine destroy_fft_solver

    implicit none

    deallocate(a, b, c, c1, d1, mwn, in, out)

  end subroutine destroy_fft_solver
  !===============================================================================

end module fft_solver
