module navier_stokes_pub

  ! Public variables for the navier stokes solver.

  interface
    subroutine event()
    end subroutine event
  end interface

  ! Maximum number of iteration
  integer :: nstep = 0, log

  ! Maximum phisycal time, time step and time
  real :: Tmax, dt, t

  ! Phisycal parameters of the solver
  real :: mu = 1.0, rho = 1.0, g = 0.0

  ! Variables of the solver: pressure, velocity field, 
  ! predicted velocity field and projection operator
  type field
    ! Array of the field
    real, dimension(:,:), allocatable :: f
    ! Boundary conditions type: periodic, dirichlet or neumann
    character(len=9) :: left, right, top, bottom
    ! Boundary conditions values
    real, dimension(:), allocatable :: l, r, t, b
    ! Identifier between center or face location
    ! 0 = cell center
    ! 1 = x face
    ! 2 = y face
    integer :: location = 0
    ! Bound coordinates for domain cylces
    integer, dimension(2) :: lo, up
  end type field

  type(field) :: u, v, p, us, vs, phi

  ! Source term in momentum equation
  real :: Sx = 0.0, Sy = 0.0

  ! Divergence tolerance
  real :: divergence_tol = 1.0e-8

  ! Procedure pointer to function event_output
  procedure(event), pointer :: event_output => Null()
  procedure(event), pointer :: event_end => Null()
  procedure(event), pointer :: event_i => Null()

  ! Flag for the Poisson solver
  character(len=3) :: poisson_solver_type

end module navier_stokes_pub
