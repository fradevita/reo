module grid_2d

  ! Double precision kind
  integer, parameter :: dp = kind(1.d0)

  ! Box definition
  integer :: nbox = 1
  type box
    ! Number of points in each direction
    integer :: nx, ny
    ! Index of the left lower and upper right corners
    integer, dimension(2) :: ilower, iupper
    ! Physical size and grid spacing
    real(kind=dp) :: Lx, Ly, delta
    ! Coordinates of the lower left point
    real(kind=dp), dimension(2) :: p0
    ! Physical position of the cell center of all poinst in the box
    real(kind=dp), dimension(:,:), allocatable :: x, y
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
    ! Boundary condition for the Poisson solver 
    character(len=9) :: left = 'neumann', right = 'neumann', top = 'neumann', bottom = 'neumann'
    ! Neighbour boxes index
    integer :: lb, rb, tb, bb
    ! Neighbour processors index
    integer :: lp, rp, tp, bp
  end type box

  type(box), dimension(:), allocatable :: boxarray

  ! Field variable definition
  type field
    ! Array of the field
    real(kind=dp), dimension(:,:), allocatable :: f
    ! Boundary conditions type
    character(len=9) :: left, right, top, bottom
    ! Boundary conditions values
    real(kind=dp), dimension(:), allocatable :: l, r, t, b
    ! Identifier between center or face location
    ! 0 = cell center
    ! 1 = x face
    ! 2 = y face
    integer :: location = 0
    ! Bound coordinates for domain cylces
    integer, dimension(2) :: lo, up
  end type field

  ! A simple 2D array, needed for creating an array
  type myarray
    real(kind=dp), dimension(:,:), allocatable :: f
  end type myarray

  ! MPI variables
  integer :: ierr, mpi_common_world, pid, num_procs

  ! Pi
  real(kind=dp), parameter :: pi = acos(-1.0)

contains

  !===============================================================================
  subroutine create_box(n, nx, ny, Lx, Ly, x0, y0)

    ! Generate the 2D computational grid

    implicit none

    include 'mpif.h'
    
    ! Input / output variables
    integer, intent(in) :: n, nx, ny
    real(kind=dp), intent(in) :: Lx, Ly, x0, y0

    ! Local variables
    integer :: i, j

    boxarray(n)%nx = nx
    boxarray(n)%ny = ny
    boxarray(n)%ilower = [1, 1]
    boxarray(n)%iupper = [nx, ny]
    boxarray(n)%Lx = Lx
    boxarray(n)%Ly = Ly
    boxarray(n)%delta = boxarray(n)%Lx / boxarray(n)%nx
    if (boxarray(n)%Lx / boxarray(n)%nx /= boxarray(n)%Ly / boxarray(n)%ny) then
      print *, 'ERROR: grid spacing not equal in x and y'
      stop
    endif
    boxarray(n)%p0 = [x0, y0]

    allocate(boxarray(n)%x(boxarray(n)%nx,boxarray(n)%ny))
    allocate(boxarray(n)%y(boxarray(n)%nx,boxarray(n)%ny))

    do j = 1,boxarray(n)%ny
      do i = 1,boxarray(n)%nx
        boxarray(n)%x(i,j) = boxarray(n)%p0(1) + (i - 0.5) * boxarray(n)%delta
        boxarray(n)%y(i,j) = boxarray(n)%p0(2) + (j - 0.5) * boxarray(n)%delta
      end do
    end do

    boxarray(n)%lp = mpi_proc_null
    boxarray(n)%rp = mpi_proc_null
    boxarray(n)%tp = mpi_proc_null
    boxarray(n)%bp = mpi_proc_null

  end subroutine create_box
  !===============================================================================

  !===============================================================================
  subroutine destroy_boxes()
    ! Free the memory allocated for each box

    implicit none

    integer :: n

    do n = 1,nbox
      deallocate(boxarray(n)%x,boxarray(n)%y)
    end do

    deallocate(boxarray)
 
  end subroutine destroy_boxes
  !===============================================================================

end module grid_2d
