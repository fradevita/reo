program reconstruction

  use grid_2d
  use volume_of_fluid
  use mpi

  implicit none

  integer :: i, j, np, n
  real(kind=dp) :: L, delta, r
  real(kind=dp) :: x, y, xc, yc
  character(len=1) :: sn
  character(len=4) :: filename

  ! First we need to initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  np = 16
  L = 0.5
  delta = L / np
  r = 0.3
  xc = 0.5
  yc = 0.5

  nbox = 4 
  allocate(boxarray(nbox))
  call create_box(1, np, np, L, L, 0.d0, 0.d0)
  call create_box(2, np, np, L, L, L, 0.d0)
  call create_box(3, np, np, L, L, 0.d0, L)
  call create_box(4, np, np, L, L, L, L)

  call init_vof()
  call set_vof('circle', xc, yc, r)

  ! Reconstruct the interface
  call reconstruct()

  ! Print the interface
  do n = 1,nbox
    write(sn,'(I1)') n
    filename = 'box'//sn
    open(10+n, file = filename)
    do j = 1,np
      y = boxarray(n)%p0(2) + (j-0.5)*delta
      do i = 1,np
        x = boxarray(n)%p0(1) + (i-0.5)*delta
        write(10+n,*) x, y, vof(n)%f(i,j), h(n)%f(i,j)
      end do
      write(10+n,*) ''
    end do
    close(10+n)
  end do

  call mpi_finalize(ierr)

end program reconstruction
