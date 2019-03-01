program euler

  use grid_2d
  use navier_stokes

  implicit none
  integer :: i, j, n
  real :: xl, yl
  character(len=4) :: arg
  
  ! Boundary conditions
  left_boundary = 'periodic'
  right_boundary = 'periodic'
  top_boundary = 'periodic'
  bottom_boundary = 'periodic'

  ! Set the number of points and the domain size
  call getarg(1,arg)
  read(arg,*) n
  nx = 2**n
  ny = 2**n
  Lx = 1.0
  Ly = 1.0
  x0 = -0.5
  y0 = -0.5

  ! Create the grid
  call create_grid()

  ! First we need to initialize the solver
  call init_ns_solver()

  ! We set the initial velocity field
  do j = 1,ny
    yl = y0 + (j - 0.5) * dy
    do i = 1,nx+1
      xl = x0 + (i - 1.0) * dx 
      u(i,j) = 1.0 - 2.0*cos(2.0*pi*xl)*sin(2.0*pi*yl)
    end do
  end do

  do j = 1,ny+1
    yl = y0 + (j - 1.0) * dy
    do i = 1,nx
      xl = x0 + (i - 0.5) * dx
      v(i,j) = 1.0 + 2.0*sin(2.0*pi*xl)*cos(2.0*pi*yl)
    end do
  end do

  ! We want to solve the euler equation so we set the viscosity to zero
  mu = 0.0

  ! We compute the solution up to time T = 0.5 and then compare the computed solution with
  ! the analytical solution
  dt = 0.35 * dx / 4.0
  Tmax = 0.5
  nstep = int(Tmax / dt)

  ! We compute the error at the end of the simulation
  event_end => output_error

  ! We run the simulation
  call solve()

contains

  subroutine output_error()
    
    ! Compute the error with respect to the analytical solution
    
    real :: s, e, e2
    
    e2 = 0.0
    do j = 1,ny
      do i = 1,nx
        s = 1.0 - 2.0 * cos(2.0*pi*(x(i,j) - dx/2.0 - t))*sin(2.0*pi*(y(i,j)-t))
        e = abs(s-u(i,j))
        e2 = e2 + e**2 * dx**2
      end do
    end do
    
    print *, nx, sqrt(e2)
    
  end subroutine output_error
  
end program euler
