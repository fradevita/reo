program poiseuille

  use grid_2d
  use navier_stokes
  
  implicit none
  include 'mpif.h'

  integer :: j, n, nx, ny
  real(kind=dp) :: y
  character(len=4) :: arg

  ! We need an auxiliary field to check steady state
  type(field), dimension(1) :: uold

  ! First we need to initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD
  
  ! Set the number of points and the domain size
  call getarg(1,arg)
  read(arg,*) n
  nx = 4*(2**n)
  ny = 2**n

  ! Create the grid
  allocate(boxarray(1))
  call create_box(1, nx, ny, 4.d0, 1.d0, -0.5d0, -0.5d0)
  
  ! Boundary conditions
  boxarray(1)%left_boundary = 'inflow'
  boxarray(1)%right_boundary = 'outflow'

  ! Select Poisson solver 
  poisson_solver_type = 'itr'

  ! First we need to initialize the solver
  call init_ns_solver()
  
  ! Boundary condition
  do j = u(1)%lo(2),u(1)%up(2)
    y = boxarray(1)%p0(2) + (j-0.5d0)*boxarray(1)%delta
    u(1)%l(j) = 0.5d0*(0.25d0 - y**2)
  end do

  ! We compute some quantites to check the properties of the scheme
  event_i => e_istep
  event_output => output

  ! We run the simulation
  dt = 0.5*dt
  call solve()

contains
  
  subroutine e_istep()

    implicit none

    uold = u
    
  end subroutine e_istep
    
  subroutine output()

    implicit none
    
    ! Check for steady-state
    real(kind=dp) :: diff, s, e, emax, y
    logical :: steady
    integer :: i, j
    
    steady = .true.
    
    do j = u(1)%lo(2),u(1)%up(2)
      do i = u(1)%lo(1),u(1)%up(1)
        diff = abs(u(1)%f(i,j) - uold(1)%f(i,j))
        if (diff > 1.0e-8) steady = .false.
      end do
    end do
   
    if (steady) then
      emax = 0.0
      ! Compute the error with respect to the analytical profile
      do j = u(1)%lo(2),u(1)%up(2)
        y = -0.5 + (j-0.5)*boxarray(1)%delta
        s = 0.5*(0.25 - y**2)
        e = abs(s-u(1)%f(nx/2,j))
        if (e > emax) emax = e
      end do
      print *, 2**n, emax
      call destroy_ns_solver()
    end if
    
  end subroutine output
  
end program poiseuille
