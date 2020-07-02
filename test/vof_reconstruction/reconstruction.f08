program reconstruction

  use grid_2d
  use volume_of_fluid
  use mpi

  implicit none

  integer :: i, j, np
  real(kind=dp) :: Lx, Ly, delta, r
  real(kind=dp) :: xc, yc, modnorm

  ! First we need to initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  np = 32
  Lx = 1.
  Ly = 1.
  delta = Lx/np
  r = 0.433
  xc = 0.5
  yc = 0.5
 
  allocate(boxarray(1))
  call create_box(1, np, np, Lx, Ly, 0.d0, 0.d0)

  call init_vof()
  call set_vof('circle', xc, yc, r)
  call compute_norm()

  ! Override the norms and its derivatives
  !do j = 1,np
  !  y = (j - 0.5)*delta
  !  do i = 1,np
  !    x = (i - 0.5)*delta
      !nx(1)%f(i,j) = (x - xc)/sqrt((x - xc)**2 + (y - yc)**2) 
      !ny(1)%f(i,j) = (y - yc)/sqrt((x - xc)**2 + (y - yc)**2)
      !lxx(1)%f(i,j) = (y - yc)**2/((x - xc)**2 + (y - yc)**2)**(3./2.)
      !lyy(1)%f(i,j) = (x - xc)**2/((x - xc)**2 + (y - yc)**2)**(3./2.)
  !  end do
  !end do

  ! Reconstruct the interface
  call reconstruct()

  ! Print the interface
  open(unit = 1, file = 'h')
  do j = 1,np
    do i = 1,np
      write(1,*) (i-0.5d0)*delta, (j-0.5d0)*delta, h(1)%f(i,j)
    end do
    write(1,*) ''
  end do
  close(1)

  ! Finalize the simulation
  call MPI_finalize(ierr)

end program reconstruction
