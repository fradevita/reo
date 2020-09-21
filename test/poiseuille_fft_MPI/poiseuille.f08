program poiseuille

  use grid_2d
  use navier_stokes
  
  implicit none
  include 'mpif.h'
  
  integer :: n, nx, ny
  character(len=4) :: arg
  real(kind=dp) :: Lx, Ly
  
  ! We need an auxiliary field to check steady state
  type(field), dimension(1) :: uold
  
  ! First we need to initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  ! Set the number of points and the domain size
  call getarg(1,arg)
  read(arg,*) n

  nx = 2**n
  ny = (2**n)/num_procs
  Lx = 1.
  Ly = 1./float(num_procs)

  ! Create the grid
  allocate(boxarray(1))
  call create_box(1, nx, ny, Lx, Ly, -0.5d0, -0.5d0 + 0.5d0*pid)

  ! Boundary conditions
  boxarray(1)%left_boundary = 'periodic'
  boxarray(1)%right_boundary = 'periodic'
  if (num_procs > 1) then
    boxarray(1)%top_boundary = 'halo'
    boxarray(1)%tp = pid+1
    boxarray(1)%bottom_boundary = 'halo'
    boxarray(1)%bp = pid-1
    boxarray(1)%ilower = [1, 1 + pid*ny]
    boxarray(1)%iupper = [nx, ny + pid*ny]
    if (pid == 0) then 
      boxarray(1)%bottom_boundary = 'no-slip'
      boxarray(1)%bp = mpi_proc_null
    endif
    if (pid == num_procs - 1) then
      boxarray(1)%top_boundary = 'no-slip'
      boxarray(1)%tp = mpi_proc_null
    endif
  endif

  ! Select the Poisson solver
  poisson_solver_type = 'fft'

  ! Initialize the solver
  call init_ns_solver()

  ! Source term
  Sx = 1.0

  ! We compute some quantites to check the properties of the scheme
  event_i => e_istep
  event_output => output

  ! We run the simulation
  nstep = 2
  call solve()

contains
  
  subroutine e_istep()

    implicit none
    integer :: i, j
    uold = u

  end subroutine e_istep
    
  subroutine output()

    implicit none
    
    ! Check for steady-state
    real(kind=dp) :: diff, maxdiff, s, e, emax, y
    logical :: steady
    integer :: i, j, imax, jmax
    
    steady = .true.
    maxdiff = 0.
    do j = u(1)%lo(2),u(1)%up(2)
      do i = u(1)%lo(1),u(1)%up(1)
        diff = abs(u(1)%f(i,j) - uold(1)%f(i,j))
        if (diff > maxdiff) maxdiff = diff
      end do
    end do
    call mpi_allreduce(mpi_in_place,maxdiff,1,mpi_real8,mpi_max,mpi_common_world,ierr)
    if (maxdiff > 1.0e-8) steady = .false.

    if (steady) then
      emax = 0.0
      ! Compute the error with respect to the analytical profile
      do j = u(1)%lo(2),u(1)%up(2)
        y = boxarray(1)%p0(2) + (j-0.5)*boxarray(1)%delta
        do i = u(1)%lo(1),u(1)%up(1)
          s = 0.5*(0.25 - y**2)
          e = abs(s-u(1)%f(i,j))
          if (e > emax) then
            emax = e
            imax = i
            jmax = j
          end if
        end do
      end do
      call mpi_allreduce(mpi_in_place,emax,1,mpi_real8,mpi_max,mpi_common_world,ierr)
      print *, 2**n, emax
      call destroy_ns_solver()
    end if
    
  end subroutine output
  
end program poiseuille
