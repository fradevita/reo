program poiseuille

  use grid_2d
  use navier_stokes
  
  implicit none
  include 'mpif.h'
  
  integer :: n, nx, ny
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

  nx = 2**n
  ny = 2**n

  ! Create the grid
  allocate(boxarray(1))
  call create_box(1, nx, ny, 1.d0, 1.d0, -0.5d0, -0.5d0)

  ! Boundary conditions
  boxarray(1)%left_boundary = 'periodic'
  boxarray(1)%right_boundary = 'periodic'

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
    logical :: steady
    integer :: i, j, imax, jmax
    
    steady = .true.
    
    do j = u(1)%lo(2),u(1)%up(2)
      do i = u(1)%lo(1),u(1)%up(1)
        diff = abs(u(1)%f(i,j) - uold(1)%f(i,j))
        if (diff > 1.0e-8) steady = .false.
      end do
    end do
   
    if (steady) then
      emax = 0.0
      ! Compute the error with respect to the analytical profile
      do j = u(1)%lo(2),u(1)%up(2)
        y = -0.5 + (j-0.5)*boxarray(1)%delta
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
      print *, 2**n, emax
      call destroy_ns_solver()
    end if
    
  end subroutine output
  
end program poiseuille
