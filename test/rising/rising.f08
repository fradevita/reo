program rising

  use grid_2d
  use navier_stokes
  use volume_of_fluid
  use multiphase
  use mpi

  implicit none
 
  integer :: npx, npy
  real(kind=dp) :: Lx, Ly


  ! First we need to initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  allocate(boxarray(1))
  if (num_procs == 1) then
    npx = 128
    npy = 256
    Lx = 1.
    Ly = 2.
    call create_box(1, npx, npy, Lx, Ly, 0.d0, 0.d0)
    boxarray(1)%left_boundary = 'free-slip'
    boxarray(1)%right_boundary = 'free-slip'
  elseif (num_procs == 4) then
    npx = 64
    npy = 128
    Lx = 0.5
    Ly = 1.
    if (pid == 0) then
      call create_box(1, npx, npy, Lx, Ly, 0.d0, 0.d0)
      boxarray(1)%left_boundary = 'free-slip'
      boxarray(1)%right_boundary = 'halo'
      boxarray(1)%top_boundary = 'halo'
      boxarray(1)%rp = 1
      boxarray(1)%tp = 2
    endif
    if (pid == 1) then
      call create_box(1, npx, npy, Lx, Ly, Lx, 0.d0)
      boxarray(1)%left_boundary = 'halo'
      boxarray(1)%right_boundary = 'free-slip'
      boxarray(1)%top_boundary = 'halo'
      boxarray(1)%lp = 0
      boxarray(1)%tp = 3
      boxarray(1)%ilower = [npx + 1, 1]
      boxarray(1)%iupper = [2*npx, npy]
    endif
    if (pid == 2) then
      call create_box(1, npx, npy, Lx, Ly, 0.d0, Ly)
      boxarray(1)%left_boundary = 'free-slip'
      boxarray(1)%right_boundary = 'halo'
      boxarray(1)%bottom_boundary = 'halo'
      boxarray(1)%rp = 3
      boxarray(1)%bp = 0
      boxarray(1)%ilower = [1, npy + 1]
      boxarray(1)%iupper = [npx, 2*npy]
    endif
    if (pid == 3) then
      call create_box(1, npx, npy, Lx, Ly, Lx, Ly)
      boxarray(1)%left_boundary = 'halo'
      boxarray(1)%right_boundary = 'free-slip'
      boxarray(1)%bottom_boundary = 'halo'
      boxarray(1)%lp = 2
      boxarray(1)%bp = 1
      boxarray(1)%ilower = [npx + 1, npy + 1]
      boxarray(1)%iupper = [2*npx, 2*npy]
    endif
  else
    print *, 'Run either with 1 or 4 processors'
    stop
  endif

  ! Select and initialize the Poisson solver 
  poisson_solver_type = 'itr'
  call init_ns_solver()

  ! Set the vof
  call init_vof
  call set_vof('circle', 0.5d0, 0.5d0, 0.25d0) 

  ! Multiphase properties
  rho0 = 1000.
  rho1 = 100.
  mu0 = 10.
  mu1 = 1.
  sigma = 24.5
  solvemultiphase = .true.

  ! Set gravity
  g(2) = -0.98

  ! Set boundary condition on pressure on top
  p(1)%t(1:npx) = rho0*g(2)

  ! Set the maximum number of iterations
  Tmax = 3.

  ! We compute some quantites to check the properties of the scheme
  event_output => output

  ! We run the simulation
  nstep = 100
  call solve()

contains

    subroutine output()

    implicit none

    integer :: i, j
    
    if (t >= Tmax) then
      !open(30,file = 'out')
      do j = 1,npy
        do i = 1,npx
          write(30+pid,*) boxarray(1)%x(i,j), boxarray(1)%y(i,j), vof(1)%f(i,j)
        end do
        write(30+pid,*) ''
      end do
      write(30+pid,*) ''
      write(30+pid,*) ''
      flush(30+pid)
    endif

  end subroutine output

end program rising
