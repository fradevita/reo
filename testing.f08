program test

  use grid_2d
  use navier_stokes
  
  implicit none
  integer :: i, j
  real :: xl, yl
  
  ! Set the number of points and the domain size
  nx = 128
  ny = 128
  Lx = 1.0
  Ly = 1.0
  x0 = -0.5
  y0 = -0.5
  
  ! Create the grid
  call create_grid()

  ! Boundary conditions
  p_left = 'Periodic'
  p_right = 'Periodic'
  p_top = 'Periodic'
  p_bottom = 'Periodic'
  
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

  ! We compute the solution up to time T = 3.0 and then compare the computed solution with
  ! the analytical solution
  dt = 0.25 * dx / 4.0
  Tmax = 3.0
  nstep = int(Tmax / dt)

  ! We compute the error at the end of the simulation
  event_end => output_error
  
  ! We run the simulation
  call solve()

contains

  subroutine output_error()

    ! Compute the error with respect to the analytical solution

    real :: s, e, emax

    emax = 0.0
    do j = 1,ny
      do i = 1,nx
        s = 1.0 - 2.0 * cos(2.0*pi*(x(i,j) - dx/2.0 - t))*sin(2.0*pi*(y(i,j)-t))
        e = abs(s-u(i,j))
        if (e > emax) emax = e
      end do
    end do

    print *, nx, emax 

  end subroutine output_error
  
end program test
