program poiseuille

  use grid_2d
  use navier_stokes

  implicit none
  include 'mpif.h'

  integer :: j, n, nx, ny
  real(kind=dp) :: y
  character(len=4) :: arg

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
  nx = 4*(2**n) / num_procs
  ny = 2**n

  ! Create the grid
  allocate(boxarray(1))
  call create_box(1, nx, ny, 4.d0 / float(num_procs), 1.d0, -0.5d0 + 2.d0*pid, -0.5d0)
 
  ! Boundary conditions
  if (num_procs == 1) then
    boxarray(1)%left_boundary = 'inflow'
    boxarray(1)%right_boundary = 'outflow'
  elseif (num_procs == 2) then
    if (pid == 0) then
      boxarray(1)%left_boundary = 'inflow'
      boxarray(1)%right_boundary = 'halo'
      boxarray(1)%rp = 1
    elseif (pid == 1) then
      boxarray(1)%left_boundary = 'halo'
      boxarray(1)%right_boundary = 'outflow'
      boxarray(1)%lp = 0
      boxarray(1)%ilower = [nx+1, 1]
      boxarray(1)%iupper = [2*nx, ny]
    endif
  endif

  ! Select Poisson solver 
  poisson_solver_type = 'itr'

  ! Initialize the solver
  call init_ns_solver()
  
  ! Boundary condition
  if (pid == 0) then
    do j = u(1)%lo(2),u(1)%up(2)
      y = boxarray(1)%p0(2) + (j-0.5)*boxarray(1)%delta
      u(1)%l(j) = 0.5*(0.25 - y**2)
    end do
  end if

  ! Check steady state and compute errors
  event_i => e_istep
  event_output => output

  ! We run the simulation
  dt = 0.5*dt
  call solve()

contains
  
  subroutine e_istep()

    implicit none

    uold = u
    
  end subroutine e_istep
    
  subroutine output()

    implicit none
    
    ! Check for steady-state
    real(kind=dp) :: diff, s, e, emax, y
    integer :: i, j, steady
    
    steady = 0
    do j = u(1)%lo(2),u(1)%up(2)
      do i = u(1)%lo(1),u(1)%up(1)
        diff = abs(u(1)%f(i,j) - uold(1)%f(i,j))
        if (diff > 1.0e-8) steady = 1
      end do
    end do
    call mpi_allreduce(mpi_in_place,steady,1,mpi_integer,mpi_sum,mpi_common_world,ierr)
   
    if (steady == 0) then
      emax = 0.0
      ! Compute the error with respect to the analytical profile
      do j = u(1)%lo(2),u(1)%up(2)
        y = -0.5 + (j-0.5)*boxarray(1)%delta
        do i = u(1)%lo(1),u(1)%up(1)-1
          s = 0.5*(0.25 - y**2)
          e = abs(s-u(1)%f(i,j))
          if (e > emax) emax = e
        end do
      end do
      call mpi_allreduce(mpi_in_place,emax,1,mpi_real8,mpi_max,mpi_common_world,ierr)
      if (pid ==0) print *, 2**n, emax
      call destroy_ns_solver()
    end if
    
  end subroutine output
  
end program poiseuille
