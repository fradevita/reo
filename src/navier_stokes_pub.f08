module navier_stokes_pub

  ! Public variables for the navier stokes solver.

  interface
    subroutine event()
    end subroutine event
  end interface

  ! Maximum number of iteration
  integer :: nstep = 0, log

  ! Maximum phisycal time
  real :: Tmax, dt, t

  ! Phisycal parameters of the solver
  real :: mu = 1.0, rho = 1.0, g = 0.0

  ! Variables of the solver: pressure, velocity field, 
  ! predicted velocity field and projection operator
  real, dimension(:,:), allocatable :: u, v, us, vs, p

  ! Source term in momentum equation
  real :: Sx = 0.0, Sy = 0.0

  ! Divergence tolerance
  real :: divergence_tol = 1.0e-8

  ! Procedure pointer to function event_output
  procedure(event), pointer :: event_output => Null()
  procedure(event), pointer :: event_end => Null()
  procedure(event), pointer :: event_i => Null()

end module navier_stokes_pub
