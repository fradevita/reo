program reconstruction

  use grid_2d
  use volume_of_fluid
  use mpi

  implicit none

  integer :: i, j, np, s, smax
  real(kind=dp) :: Lx, Ly, delta, r
  real(kind=dp) :: x, y, xc, yc, dt
  type(field) :: u(1), v(1)

  ! Initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  np = 100
  Lx = pi
  Ly = pi
  delta = Lx / np
  r = 0.2*pi
  xc = 0.5*pi
  yc = 0.2*(pi + 1)
 
  allocate(boxarray(1))
  call create_box(1, np, np, Lx, Ly, 0.d0, 0.d0)

  call init_vof()
  call set_vof('circle', xc, yc, r)

  ! Define velocity field
  allocate(u(1)%f(np+1,np))
  allocate(v(1)%f(np,np+1))
  do j = 1,np
    do i = 1,np+1
      u(1)%f(i,j) = sin((i-1)*delta)*cos((j-0.5)*delta)
    end do
  end do
  do j = 1,np+1
    do i = 1,np
      v(1)%f(i,j) = -cos((i-0.5)*delta)*sin((j-1)*delta)
    end do
  end do

  ! Timestep
  dt = 0.0025*pi

  ! Number of steps
  smax = int(10*pi/dt)

  ! Print the interface at t = 0
  open(unit = 1, file = 'vof')
  do j = 1,np
    y = (j-0.5)*delta
    do i = 1,np
      x = (i-0.5)*delta
      write(1,*) x, y, vof(1)%f(i,j)
    end do
    write(1,*) ''
  end do
  write(1,*) ''
  write(1,*) ''

  ! Time cycle
  do s = 1,smax
    ! When t = T/2 reverse the velocity field and save the profile 
    if (s == smax/2) then
      do j = 1,np
        do i = 1,np+1
          u(1)%f(i,j) = -sin((i-1.)*delta)*cos((j-0.5)*delta)
        end do
      end do
      do j = 1,np+1
        do i = 1,np
          v(1)%f(i,j) = cos((i-0.5)*delta)*sin((j-1)*delta)
        end do
      end do

      ! Print the interface at t = T/2
      do j = 1,np
        y = (j-0.5)*delta
        do i = 1,np
          x = (i-0.5)*delta
          write(1,*) x, y, vof(1)%f(i,j)
        end do
        write(1,*) ''
      end do
      write(1,*) ''
      write(1,*) ''
    endif 

    call advect_vof(u, v, dt)
  end do

  ! Print the interface at t = T
  do j = 1,np
    y = (j-0.5)*delta
    do i = 1,np
      x = (i-0.5)*delta
      write(1,*) x, y, vof(1)%f(i,j)
    end do
    write(1,*) ''
  end do
  write(1,*) ''
  write(1,*) ''
  close(1)

  ! Finalize the simulation
  call MPI_finalize(ierr)

end program reconstruction
