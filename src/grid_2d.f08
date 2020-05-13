module grid_2d

  ! Number of grid points in one direction
  integer :: nx, ny

  ! Size of the domain and grid spacing
  real :: Lx = 1.0, Ly = 1.0, x0 = 0., y0 = 0., dx, dy

  ! Center of each grid cell
  real, dimension(:,:), allocatable :: x, y
  
  ! MPI variables
  integer :: ierr, mpi_common_world, myid, num_procs

  ! Pi
  real, parameter :: pi = acos(-1.0)
    
  ! Boundary conditions type, can be:
  ! no-slip: Dirichlet 0 on normal and tangential velocity, Neumann 0 on pressure
  ! free-slip: Dirichlet 0 on normal velocity, Dirichlet on tangential velocity,
  !             Neumann 0 on pressure
  ! inflow: Dirichlet on normal velocity, Dirichlet 0 on tangential velocity,
  !           Neumann 0 on pressure
  ! outflow: Neumann 0 on normal velocity, Dirichlet 0 on tangential velocity,
  !          Dirichlet 0 on pressure
  ! periodic: periodic in that direction
  ! By defaul all boundary are set to no-slip
  character(len=9) :: left_boundary = 'no-slip', right_boundary = 'no-slip'
  character(len=9) :: top_boundary = 'no-slip', bottom_boundary = 'no-slip'

contains

  subroutine create_grid()

    ! Generate the 2D computational grid
    
    implicit none
    
    integer :: i, j

    dx = Lx / nx
    dy = Ly / ny

    allocate(x(nx,ny))
    allocate(y(nx,ny))

    do j = 1,ny
      do i = 1,nx
        x(i,j) = x0 + (i - 0.5) * dx
        y(i,j) = y0 + (j - 0.5) * dy
      end do
    end do
        
  end subroutine create_grid

  subroutine destroy_grid()

    implicit none

   deallocate(x,y)
  
  end subroutine destroy_grid

end module grid_2d
