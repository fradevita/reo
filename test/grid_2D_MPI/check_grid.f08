program check_grid

  ! Check consistency one generating 2D grids

  use grid_2D

  implicit none
  
  include 'mpif.h'
  
  integer :: n, i, j
  real(kind=dp) :: a, b
  character(len=1) :: sn
  character(len=4) :: filename

  ! Initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  ! Now do the same for three grids 
  nbox = 1
  allocate(boxarray(nbox))
  do n = 1,nbox
    call create_box(n, 20, 20, 1.d0, 1.d0, 1.d0*pid, 1.d0*pid)

    print *, 'Grid ', n, 'of ', nbox, 'on processor ', pid
    print *, 'Sizes: ', boxarray(n)%Lx, boxarray(n)%Ly, boxarray(n)%nx, boxarray(n)%ny, boxarray(n)%delta
    print *, 'Origin: ', boxarray(n)%p0
    print *, ' x and y array size: ', size(boxarray(n)%x), size(boxarray(n)%y)
    print *, 'Physical BC: '
    print *, '             left   ', boxarray(n)%left_boundary
    print *, '             right  ', boxarray(n)%right_boundary
    print *, '             top    ', boxarray(n)%top_boundary
    print *, '             bottom ', boxarray(n)%bottom_boundary
    print *, 'Interanl BC: '
    print *, '             left   ', boxarray(n)%left
    print *, '             right  ', boxarray(n)%right
    print *, '             top    ', boxarray(n)%top
    print *, '             bottom ', boxarray(n)%bottom
    print *, ''

    ! Try to access all position in the grid
    write(sn,'(I1)') pid+1
    filename = 'box'//sn
    open(unit = pid+1, file = filename) 
    do j = 1,boxarray(n)%ny
      do i = 1,boxarray(n)%nx
        a = boxarray(n)%x(i,j)
        b = boxarray(n)%y(i,j)
        write(pid+1,*) boxarray(n)%x(i,j), boxarray(n)%y(i,j)
      end do
    end do
    flush(pid+1)
    close(pid+1)

  end do

  call destroy_boxes()

  call MPI_finalize(ierr)

end program

