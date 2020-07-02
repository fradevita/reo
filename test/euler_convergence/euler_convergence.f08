program euler

  use grid_2d
  use navier_stokes

  implicit none

  include 'mpif.h'

  integer :: i, j, n, nx, ny
  real(kind=dp) :: x, y
  character(len=4) :: arg

  ! Initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  ! Set the number of points and the domain size
  call getarg(1,arg)
  read(arg,*) n
  nx = 2**n
  ny = 2**n

  ! Create the grid
  allocate(boxarray(1))
  call create_box(1, nx, ny, 1.d0, 1.d0, -0.5d0, -0.5d0)
  
  ! Boundary conditions
  boxarray(1)%left_boundary = 'periodic'
  boxarray(1)%right_boundary = 'periodic'
  boxarray(1)%top_boundary = 'periodic'
  boxarray(1)%bottom_boundary = 'periodic'

  ! Select the Poisson solver
  poisson_solver_type = 'itr'

  ! Set the viscosity to zero
  viscosity = 0.

  ! Initialize the solver
  call init_ns_solver()

  ! Set the initial velocity field
  do j = u(1)%lo(2),u(1)%up(2)
    y = boxarray(1)%p0(2) + (j - 0.5) * boxarray(1)%delta
    do i = u(1)%lo(1),u(1)%up(1)
      x = boxarray(1)%p0(1) + (i - 1.0) * boxarray(1)%delta
      u(1)%f(i,j) = 1.0 - 2.0*cos(2.0*pi*x)*sin(2.0*pi*y)
    end do
  end do

  do j = v(1)%lo(2),v(1)%up(2)
    y = boxarray(1)%p0(2) + (j - 1.0) * boxarray(1)%delta
    do i = v(1)%lo(1),v(1)%up(1)
      x = boxarray(1)%p0(1) + (i - 0.5) * boxarray(1)%delta
      v(1)%f(i,j) = 1.0 + 2.0*sin(2.0*pi*x)*cos(2.0*pi*y)
    end do
  end do

  ! Compute the solution up to time T = 0.5 and then compare the computed solution with
  ! the analytical solution
  dt = 0.35 * boxarray(1)%delta / 4.0
  Tmax = 0.5

  ! Compute the error at the end of the simulation
  event_end => output_error

  ! We run the simulation
  call solve()

contains
    
  subroutine output_error()
    
    ! Compute the error with respect to the analytical solution

    implicit none    

    real(kind=dp) :: s, e, e2, delta

    delta = boxarray(1)%delta

    e2 = 0.0
    do j = u(1)%lo(2),u(1)%up(2)
      y = boxarray(1)%p0(2) + (j - 0.5)*delta
      do i = u(1)%lo(1),u(1)%up(1)
         x = boxarray(1)%p0(1) + (i - 1.0)*delta
        s = 1.0 - 2.0 * cos(2.0*pi*(x - t))*sin(2.0*pi*(y - t))
        e = abs(s-u(1)%f(i,j))
        e2 = e2 + e**2 * delta**2
      end do
    end do
    
    print *, nx, sqrt(e2)
    
  end subroutine output_error
  
end program euler
