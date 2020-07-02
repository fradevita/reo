program poisson_periodic

  use grid_2d  
  use hypre_solver
  
  implicit none

  include 'mpif.h'  

  integer :: i, j, n, nx, ny
  real(kind=dp) :: s, e, emax, Lx, Ly, x0, y0
  type(myarray), dimension(:), allocatable :: rhs

  open(unit = 1, file = 'error')

  ! First we need to initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  do n = 3,8
    ! We set the number of points and the domain size
    nx = 2**n
    ny = 2**n
    Lx = 1.0
    Ly = 1.0
    x0 = 0.
    y0 = 0.
    
    ! Specify that the solution is periodic
    periodic(1) = nx
    periodic(2) = ny
    
    ! and create the computational grid
    allocate(boxarray(1))
    call create_box(1, nx, ny, Lx, Ly, x0, y0)

    ! Set boundary condition on the box
    boxarray(1)%left = 'periodic'
    boxarray(1)%right = 'periodic'
    boxarray(1)%top = 'periodic'
    boxarray(1)%bottom = 'periodic'

    ! Then we initialize the hypre solver
    call init_hypre_solver()

    ! We give the RHS of the poisson equation
    allocate(rhs(nbox))
    allocate(rhs(1)%f(nx,ny))
    do j = 1,ny
      do i = 1,nx
        rhs(1)%f(i,j) = -8.0*pi*pi*sin(2.*pi*boxarray(1)%x(i,j))*sin(2.*pi*boxarray(1)%y(i,j))
      end do
    end do
 
    ! solve the poisson equation
    call solve_poisson_hypre(rhs)

    ! and compute the error with analytical solution
    emax = 0.0
    do j = 1,ny
      do i = 1,nx
        s = sin(2.*pi*boxarray(1)%x(i,j))*sin(2.*pi*boxarray(1)%y(i,j))
        e = abs(s - rhs(1)%f(i,j))
        if (e .gt. emax) emax = e
      end do
    end do
    write(1,*) nx, emax

    ! Free the memory
    deallocate(rhs)
    call destroy_hypre_solver
    call destroy_boxes
    
  end do
  close(1)
  
  ! Finalize the simulation
  call MPI_finalize(ierr)

end program poisson_periodic
