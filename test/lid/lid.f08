program lid

  use grid_2d
  use navier_stokes

  implicit none
  include 'mpif.h'

  integer :: nx, ny

  ! We need an auxiliary field to check steady state
  type(field), dimension(1) :: uold

  ! First we need to initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  ! Set the number of points and the domain size
  nx = 64
  ny = 64

  ! Create the grid
  allocate(boxarray(1))
  call create_box(1, nx, ny, 1.d0, 1.d0, -0.5d0, -0.5d0)

  ! Select Poisson solver
  poisson_solver_type = 'itr'
  
  ! We set the viscosity to 1
  viscosity = 1.0e-3

  ! First we need to initialize the solver
  call init_ns_solver()

  ! Boundary conditions
  u(1)%t = 1.0

  ! We check for steady state
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
    real(kind=dp) :: diff, uc, vc, diffmax, vort, x, y
    logical :: steady
    integer :: i, j, out_id, xprof, yprof
    
    steady = .true.
    diffmax = 0.0
 
    do j = u(1)%lo(2),u(1)%up(2)
      do i = u(1)%lo(1),u(1)%up(1)
        diff = abs(u(1)%f(i,j) - uold(1)%f(i,j))
        if (diff > 1.0e-5) steady = .false.
        if (diff > diffmax) diffmax = diff
      end do
    end do

    if (steady) then
      ! Output the horizontal porfile of y-component of the velocity and
      ! vertical profile of x-component of the velocity on the centerline of the box
      open(newunit = xprof, file = 'xprof')
      i = nx / 2
      do j = 1,ny
        y = boxarray(1)%p0(2) + (j - 0.5)*boxarray(1)%delta
        write(xprof,*) y, 0.5*(u(1)%f(i,j) + u(1)%f(i+1,j))
      end do
      close(xprof)

      open(newunit = yprof, file = 'yprof')
      j = ny / 2
      do i = 1,nx
        x = boxarray(1)%p0(1) + (i - 0.5)*boxarray(1)%delta
        write(yprof,*) x, 0.5*(v(1)%f(i,j) + v(1)%f(i,j+1))
      end do
      close(yprof)

      ! Output also the norm of the velocity field at the steady state
      ! and vorticity
      open(newunit = out_id, file = 'out')
      do j = 1,ny
        do i = 1,nx
          uc = 0.5*(u(1)%f(i,j) + u(1)%f(i+1,j))
          vc = 0.5*(v(1)%f(i,j) + v(1)%f(i,j+1))
          vort = 0.25*(u(1)%f(i,j+1) + u(1)%f(i+1,j+1) - &
                       u(1)%f(i,j-1) - u(1)%f(i+1,j-1)) / boxarray(1)%delta - &
                 0.25*(v(1)%f(i+1,j) + v(1)%f(i+1,j+1) - &
                       v(1)%f(i-1,j+1) - v(1)%f(i-1,j)) / boxarray(1)%delta
          
          write(out_id,*) i, j, sqrt(uc**2 + vc**2), vort
        end do
        write(out_id,*) ''
      end do
      close(out_id)
      print *, 'Steady State'
      call destroy_ns_solver()
    end if
    
  end subroutine output
  
end program lid
