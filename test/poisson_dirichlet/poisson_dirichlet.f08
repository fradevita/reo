program poisson_dirichlet

  use grid_2d  
  ! We want to solve one simple poisson equation using the hypre solver
  ! We need to include the module
  use hypre_solver
  
  implicit none
  include 'mpif.h'  
  integer :: i, j, n
  real :: s, e, emax
  real, dimension(:), allocatable :: rhs
  
  open(unit = 1, file = 'error')
  
  ! First we need to initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  do n = 3,8

    ! We set the number of points and the domain size
    nx = 2**n
    ny = 2**n
    Lx = 1.0
    Ly = 1.0

    ! and create the computational grid
    call create_grid()

    ! Then we initialize the hypre solver
    left_boundary = 'outflow'
    right_boundary = 'outflow'
    top_boundary = 'outflow'
    bottom_boundary = 'outflow'
    call init_hypre_solver()

    ! We give the RHS of the poisson equation
    allocate(rhs(nx*ny))
    do j = 1,ny
      do i = 1,nx
        rhs(i + (j-1)*ny) = -pi*pi*18.*sin(3.*pi*x(i,j))*sin(3.*pi*y(i,j))
      end do
    end do

    ! We set the tolerance to 1.0e-30 and we add more information in output
    tolerance = 1.0e-30
    
    ! solve the poisson equation
    call solve_poisson(rhs)

    ! and compute the error with analytical solution
    emax = 0.0
    do j = 1,ny
      do i = 1,nx
        s = sin(pi*3*x(i,j))*sin(pi*3*y(i,j))
        e = abs(s - rhs(i+(j-1)*ny))
        if (e .gt. emax) emax = e
      end do
    end do
    write(1,*) nx, emax

    ! Free the memory
    deallocate(rhs)
    call destroy_hypre_solver
    call destroy_grid
    
  end do
  close(1)
  
  ! Finalize the simulation
  call MPI_finalize(ierr)

end program poisson_dirichlet
