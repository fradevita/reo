program poiseuille

  use grid_2d
  use navier_stokes
  
  implicit none
  include 'mpif.h'

  integer :: j, n, nx, ny
  real :: y
  character(len=4) :: arg

  ! We need an auxiliary field to check steady state
  type(field), dimension(:), allocatable :: uold

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
  nbox = 2
  allocate(boxarray(nbox))
  call create_box(1, nx/2, ny, 2.0, 1.0, -0.5, -0.5)
  call create_box(2, nx/2, ny, 2.0, 1.0,  1.5, -0.5)
  
  ! Boundary conditions
  boxarray(1)%left_boundary = 'inflow'
  boxarray(1)%right_boundary = 'internal'
  boxarray(1)%rb = 2
  boxarray(2)%left_boundary = 'internal'
  boxarray(2)%right_boundary = 'outflow'
  boxarray(2)%lb = 1
  boxarray(2)%ilower = [nx/2+1, 1]
  boxarray(2)%iupper = [nx, ny]

  allocate(uold(nbox))
 
  ! Select Poisson solver 
  poisson_solver_type = 'itr'

  ! First we need to initialize the solver
  call init_ns_solver()
  
  ! We set the viscosity to 1
  viscosity = 1.0

  ! Boundary condition for the left box
  do j = u(1)%lo(2),u(1)%up(2)
    y = boxarray(1)%p0(2) + (j-0.5)*boxarray(1)%delta
    u(1)%l(j) = 0.5*(0.25 - y**2)
  end do

  ! Set the timestep and maximum number of iterations
  dt = 0.05*(boxarray(1)%delta**2 + boxarray(1)%delta**2)/viscosity
  nstep = 1000000

  ! We compute some quantites to check the properties of the scheme
  event_i => e_istep
  event_output => output

  ! We run the simulation
  call solve()

contains
  
  subroutine e_istep()

    implicit none

    uold = u
    
  end subroutine e_istep
    
  subroutine output()

    implicit none
    
    ! Check for steady-state
    real :: diff, s, e, emax, y
    logical :: steady
    integer :: i, j, nb

    steady = .true.

    do nb = 1,nbox
      do j = u(nb)%lo(2),u(nb)%up(2)
        do i = u(nb)%lo(1),u(nb)%up(1)
          diff = abs(u(nb)%f(i,j) - uold(nb)%f(i,j))
          if (diff > 1.0e-8) steady = .false.
        end do
      end do
    end do

    if (steady) then
      emax = 0.0
      do nb = 1,nbox
        ! Compute the error with respect to the analytical profile
        do j = u(nb)%lo(2),u(nb)%up(2)
          y = -0.5 + (j-0.5)*boxarray(1)%delta
          do i = u(nb)%lo(1),u(nb)%up(1)-1
            s = 0.5*(0.25 - y**2)
            e = abs(s-u(nb)%f(i,j))
            if (e > emax) emax = e
          end do
        end do
      end do
      print *, 2**n, emax
      call destroy_ns_solver()
    end if
    
  end subroutine output
  
end program poiseuille
