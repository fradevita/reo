program poiseuille

  use grid_2d
  include 'navier_stokes.h'
  
  implicit none
  
  integer :: i, j, n
  character(len=4) :: arg

  ! We need an auxiliary field to check steady state
  real, dimension(:,:), allocatable :: uold
  
  ! Boundary conditions
  left_boundary = 'inflow'
  right_boundary = 'outflow'
  
  ! Set the number of points and the domain size
  call getarg(1,arg)
  read(arg,*) n
  print *, n
  nx = 4*(2**n)
  ny = 2**n
  Lx = 4.0
  Ly = 1.0
  x0 = -0.5
  y0 = -0.5

  allocate(uold(nx+1,ny))

  ! Create the grid
  call create_grid()

  ! First we need to initialize the solver
  call init_ns_solver()
  
  ! We set the viscosity to 1
  mu = 1.0

  ! Boundary condition
  do i = u%lo(2),u%up(2)
    u%l(i) = 0.5*(0.25 - (y(1,i))**2)
  end do

  ! Set the timestep and maximum number of iterations
  dt = 0.05*(dx**2 + dy**2) / mu
  nstep = 1000000

  ! We compute some quantites to check the properties of the scheme
  event_i => e_istep
  event_output => output

  ! We run the simulation
  call solve()

contains
  
  subroutine e_istep()

    do j = u%lo(2),u%up(2)
      do i = u%lo(1),u%up(1)
        uold(i,j) = u%f(i,j)
      end do
    end do
    
  end subroutine e_istep
    
  subroutine output()

    implicit none
    
    ! Check for steady-state
    real :: diff, s, e, emax
    logical :: steady
    integer :: i, j
    
    steady = .true.
    
    do j = u%lo(2),u%up(2)
      do i = u%lo(1),u%up(1)
        diff = abs(u%f(i,j) - uold(i,j))
        if (diff > 1.0e-8) steady = .false.
      end do
    end do
   
    if (steady) then
      emax = 0.0
      ! Compute the error with respect to the analytical profile
      do j = u%lo(2),u%up(2)
        do i = u%lo(1),u%up(1)-1
          s = 0.5*(0.25 - (y(i,j))**2)
          e = abs(s-u%f(i,j))
          if (e > emax) emax = e
        end do
      end do
      print *, 2**n, emax
      call destroy_ns_solver()
    end if
    
  end subroutine output
  
end program poiseuille
