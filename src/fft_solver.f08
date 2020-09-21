module fft_solver

  use grid_2D, only : dp

  include 'fftw3.f'

  private
  public :: init_fft_solver, solve_poisson_fft, destroy_fft_solver
  
  integer(8) :: pf, pb
  integer :: nxtot, nytot
  real(kind=dp), dimension(:), allocatable :: a, b, c, c1, d1, mwn
  double complex, allocatable :: in(:), out(:)

contains

  !===============================================================================
  subroutine init_fft_solver()

    use grid_2D, only : boxarray, pi, dp, mpi_common_world, ierr
    use mpi

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
    call mpi_allreduce(ny,nytot,1,mpi_integer,mpi_sum,mpi_common_world,ierr)
    allocate(a(nytot))
    allocate(b(nytot))
    allocate(c(nytot))
    allocate(c1(nytot))
    allocate(d1(nytot))
    do j = 1,nytot
      a(j) = 1.0/delta**2 
      b(j) = -2.0/delta**2
      c(j) = 1.0/delta**2
    end do

    c1(1) = c(1)/b(1)
    do j = 2,nytot
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

    use grid_2D, only : boxarray, myarray, dp, num_procs, pid
    use boundary_conditions

    implicit none

    ! Input / Output variables
    type(myarray), intent(inout) :: rhs(:)

    ! Local variables
    integer :: k, j, nx, ny, i
    real(kind=dp), dimension(:,:), allocatable :: f, fx, fy

    nxtot = boxarray(1)%nx
    nx = nxtot/num_procs
    ny = boxarray(1)%ny

    allocate(f(nxtot,ny))

    ! For every j compute the FFT of the RHS in the x direction
    ! and store it in f
    do j = 1,ny
      in = rhs(1)%f(:,j)
      call dfftw_execute_dft(pf, in, out)
      f(:,j) = sqrt(real(out)**2 + imag(out)**2) / nxtot
      do i = 1,nxtot
        write(10+pid,*) boxarray(1)%p0(1) + (i - 0.5)*boxarray(1)%delta, &
                        boxarray(1)%p0(2) + (j - 0.5)*boxarray(1)%delta, f(i,j), rhs(1)%f(i,j)
      end do
      write(10+pid,*) ''
    end do
    write(10+pid,*) ''
    write(10+pid,*) ''

    ! If parallel transpose
    if (num_procs > 1) then
      allocate(fx(nxtot,ny))
      allocate(fy(nx,nytot))
      fx = f
      call transposexy(fx, fy)
      deallocate(f)
      allocate(f(nx,nytot))
      f = fy
    endif

    ! Solve tridiagonal system for every i
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

    ! If parallel transpose
    if (num_procs > 1) then
      fy = f
      call transposeyx(fy, fx)
      deallocate(f)
      allocate(f(nxtot,ny))
      f = fx
    endif

    ! For every j compute inverse transform of f
    do j = 1,ny
      in = f(:,j)
      call dfftw_execute_dft(pb, in, out)
      do k = 1,nxtot
        rhs(1)%f(k,j) = sqrt(real(out(k))**2 + imag(out(k))**2) / nxtot
        write(20+pid,*) boxarray(1)%p0(1) + (k - 0.5)*boxarray(1)%delta, &
                        boxarray(1)%p0(2) + (j - 0.5)*boxarray(1)%delta, rhs(1)%f(k,j)
      end do
      write(20+pid,*) ''
    end do
    write(20+pid,*) ''
    write(20+pid,*) ''
    deallocate(f)

  end subroutine solve_poisson_fft
  !===============================================================================

  !===============================================================================
  subroutine transposexy(fx, fy)

    ! Transpose data from x distribution to y distribution

    use grid_2d, only : boxarray, num_procs, pid, ierr, mpi_common_world
    use mpi

    implicit none

    real(kind=dp), intent(in) :: fx(:,:)
    real(kind=dp), intent(out) :: fy(:,:)

    integer :: nx, ny, i, ii, j, jj, n, p, send_to, recv_from
    integer :: status(MPI_STATUS_SIZE)

    nx = boxarray(1)%nx/num_procs
    ny = boxarray(1)%ny

    do p = 0,num_procs-1
      do n = 0,num_procs-1

        send_to = n
        recv_from = p

        if (pid == p) then
          do j = 1,ny
            do i = 1+nx*n,nx*(n + 1)
              call mpi_send(fx(i,j), 1, mpi_real8, send_to, n, &
                            mpi_common_world, ierr)
            end do
          end do
        elseif (pid == n) then
          do j = 1+nx*n,nx*(n + 1)
            do i = 1,nx
              call mpi_recv(fy(i,j), 1, mpi_real8, recv_from, n, &
                            mpi_common_world, status, ierr)
            end do
          end do
        end if

      end do
    end do

    jj = 0
    do j = pid*nx+1,nx*(pid+1)
      jj = jj + 1
      ii = pid*nx+1
      do i = 1,nx
        fy(i,j) = fx(ii,jj)
       ii = ii + 1
      end do
    end do

  end subroutine transposexy
  !===============================================================================

  !===============================================================================
  subroutine transposeyx(fy, fx)

    ! Transpose data from y distribution to x distribution

    use grid_2d, only : boxarray, num_procs, pid, ierr, mpi_common_world
    use mpi

    implicit none

    real(kind=dp), intent(in) :: fy(:,:)
    real(kind=dp), intent(out) :: fx(:,:)

    integer :: nx, ny, i, ii, j, jj, n, p, send_to, recv_from
    integer :: status(MPI_STATUS_SIZE)

    nx = boxarray(1)%nx/num_procs
    ny = boxarray(1)%ny

    do p = 0,num_procs-1
      do n = 0,num_procs-1

        send_to = n
        recv_from = p

        if (pid == p) then
          do j = 1+nx*n,nx*(n+1)
            do i = 1,nx
              call mpi_send(fy(i,j), 1, mpi_real8, send_to, n, &
                            mpi_common_world, ierr)
            end do
          end do
        elseif (pid == n) then
          do j = 1,ny
            do i = 1+nx*n,nx*(n + 1)
              call mpi_recv(fx(i,j), 1, mpi_real8, recv_from, n, &
                            mpi_common_world, status, ierr)
            end do
          end do
        end if

      end do
    end do

    jj = pid*nx
    do j = 1,ny
      jj = jj + 1
      ii = 1
      do i = pid*nx+1,nx*(pid+1)
        fx(i,j) = fy(ii,jj)
       ii = ii + 1
      end do
    end do

  end subroutine transposeyx
  !===============================================================================

  !===============================================================================
  subroutine destroy_fft_solver

    implicit none

    deallocate(a, b, c, c1, d1, mwn, in, out)

  end subroutine destroy_fft_solver
  !===============================================================================

end module fft_solver
