program poisson_dirichlet

  use grid_2d 
  use hypre_solver
  
  implicit none

  include 'mpif.h'  
  
  integer :: i, j, n, nx, ny
  real(kind=dp) :: Lx, Ly, s, e, emax, x0, y0
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
    x0 = 0.0
    y0 = 0.0

    ! and create the computational grid
    allocate(boxarray(1))
    call create_box(1, nx, ny, Lx, Ly, x0, y0)

    ! Set boundary conditions for the Poisson solver on the box
    boxarray(1)%left = 'dirichlet'
    boxarray(1)%right = 'dirichlet'
    boxarray(1)%top = 'dirichlet'
    boxarray(1)%bottom = 'dirichlet'

    ! Then initialize the hypre solver
    call init_hypre_solver()

    ! Set the RHS of the poisson equation
    allocate(rhs(nbox))
    allocate(rhs(1)%f(nx,ny))
    do j = 1,ny
      do i = 1,nx
        rhs(1)%f(i,j) = -pi*pi*18.*sin(3.*pi*boxarray(1)%x(i,j))*sin(3.*pi*boxarray(1)%y(i,j))
      end do
    end do
    
    ! solve the poisson equation
    call solve_poisson_hypre(rhs)

    ! and compute the error with analytical solution
    emax = 0.0
    do j = 1,ny
      do i = 1,nx
        s = sin(pi*3*boxarray(1)%x(i,j))*sin(pi*3*boxarray(1)%y(i,j))
        e = abs(s - rhs(1)%f(i,j))
        if (e .gt. emax) emax = e
      end do
    end do
    write(1,*) nx, emax

    ! Free the memory
    deallocate(rhs(1)%f)
    deallocate(rhs)
    call destroy_hypre_solver
    call destroy_boxes

  end do
  close(1)
  
  ! Finalize the simulation
  call MPI_finalize(ierr)

end program poisson_dirichlet
