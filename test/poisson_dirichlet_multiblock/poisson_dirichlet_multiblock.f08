program poisson_dirichlet_multiblock

  use grid_2d  
  use hypre_solver
  
  implicit none

  include 'mpif.h'  

  integer :: i, j, n, imax, jmax, nx, ny, nb
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

    ! Divide the domain into two boxes in the x direction
    nbox = 2
    allocate(boxarray(2))
    call create_box(1, nx/2, ny, Lx/2, Ly, x0, y0)
    boxarray(1)%left = 'dirichlet'
    boxarray(1)%right = 'internal'
    boxarray(1)%top = 'dirichlet'
    boxarray(1)%bottom = 'dirichlet'

    call create_box(2, nx/2, ny, Lx/2, Ly, Lx/2.0, y0)
    boxarray(2)%ilower = [nx/2+1, 1]
    boxarray(2)%iupper = [nx, ny]
    boxarray(2)%left = 'internal'
    boxarray(2)%right = 'dirichlet'
    boxarray(2)%top = 'dirichlet'
    boxarray(2)%bottom = 'dirichlet'

    ! Then we initialize the hypre solver
    call init_hypre_solver()

    ! We give the RHS of the poisson equation
    allocate(rhs(nbox))
    do nb = 1,nbox
      allocate(rhs(nb)%f(boxarray(nb)%nx,boxarray(nb)%ny))
      do j = 1,boxarray(nb)%ny
        do i = 1,boxarray(nb)%nx
          rhs(nb)%f(i,j) = -pi*pi*18.*sin(3.*pi*boxarray(nb)%x(i,j))*sin(3.*pi*boxarray(nb)%y(i,j))
        end do
      end do
    end do

    ! solve the poisson equation
    call solve_poisson_hypre(rhs)

    ! and compute the error with analytical solution
    emax = 0.0
    do nb = 1,nbox
      do j = 1,boxarray(nb)%ny
        do i = 1,boxarray(nb)%nx
          s = sin(pi*3*boxarray(nb)%x(i,j))*sin(pi*3*boxarray(nb)%y(i,j))
          e = abs(s - rhs(nb)%f(i,j))
          if (e .gt. emax) then 
            emax = e
            imax = i
            jmax = j
          endif
        end do
      end do
    end do
    write(1,*) nx, emax, imax, jmax

    ! Free the memory
    do nb = 1,nbox
      deallocate(rhs(nb)%f)
    end do
    deallocate(rhs)
    call destroy_hypre_solver
    call destroy_boxes
    
  end do
  close(1)
  
  ! Finalize the simulation
  call MPI_finalize(ierr)

end program poisson_dirichlet_multiblock
