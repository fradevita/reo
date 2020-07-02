module navier_stokes

  ! This module perform all the operations for solving the navier stokes equations.
  ! The steps are:
  ! 1-> compute the predicted velocity field
  ! 2-> solve the poisson equation
  ! 3-> correct the velocity field
  ! 4-> update the pressure

  use grid_2d
  use boundary_conditions
  use hypre_solver
  use fft_solver
  use volume_of_fluid, only : advect_vof, destroy_vof, vof, curv
  use multiphase
  use mpi
 
  implicit none

  interface
    subroutine poisson_solver(rhs)
      Import :: myarray
      type(myarray), intent(inout) :: rhs(:)
    end subroutine
  end interface

  interface
    subroutine event()
    end subroutine event
  end interface

  ! Maximum number of iteration
  integer :: nstep = 100000, log

  ! Maximum phisycal time, time step and time
  real(kind=dp) :: Tmax = 1000., dt, t

  ! Phisycal parameters of the solver
  real(kind=dp) :: viscosity = 1., density = 1., g(2) = 0.0
  type(field), dimension(:), allocatable :: mu, rho

  ! Array of fields, the dimension is the number of boxes
  type(field), dimension(:), allocatable :: u, v, p, us, vs, phi

  ! Source term in momentum equation
  real(kind=dp) :: Sx = 0.0, Sy = 0.0

  ! Divergence tolerance
  real(kind=dp) :: divergence_tol = 1.0e-8

  ! Procedure pointer to function event_output
  procedure(event), pointer :: event_output => Null()
  procedure(event), pointer :: event_end => Null()
  procedure(event), pointer :: event_i => Null()

  ! Flag for the Poisson solver
  character(len=3) :: poisson_solver_type
  
  private
  public :: nstep, Tmax, dt, t, mu, rho, g, u, v, p, Sx, Sy, viscosity, density
  public :: event_i, event_output, event_end
  public :: init_ns_solver, solve, destroy_ns_solver, poisson_solver_type
  
  procedure(poisson_solver), pointer :: solve_poisson
  
contains
 
  !===============================================================================
  subroutine init_ns_solver

    implicit none

    integer :: n
    real(kind=dp) :: dtconv, dtvisc, dtsurf

    ! Allocate the field array based on the number of boxes in the domain
    allocate(u(nbox)) 
    allocate(v(nbox))
    allocate(p(nbox))
    allocate(us(nbox))
    allocate(vs(nbox))
    allocate(phi(nbox))
    allocate(mu(nbox))
    allocate(rho(nbox))
 
    ! Cycle over all the boxes in the domain 
    box_cycle : do n = 1,nbox
 
      ! Based on the type of boundary, compute the extension for each field
      u(n)%lo = 1
      us(n)%lo = 1
      if (boxarray(n)%left_boundary == 'no-slip' .or. boxarray(n)%left_boundary == 'free-slip' .or. &
          boxarray(n)%left_boundary == 'inflow'  .or. boxarray(n)%left_boundary == 'periodic' .or. &
          boxarray(n)%left_boundary == 'internal' .or. boxarray(n)%left_boundary == 'halo') then
        u(n)%lo(1) = 2
        us(n)%lo(1) = 2
      endif
      u(n)%up(1) = boxarray(n)%nx + 1
      us(n)%up(1) = boxarray(n)%nx + 1
      u(n)%up(2) = boxarray(n)%ny
      us(n)%up(2) = boxarray(n)%ny
      if (boxarray(n)%right_boundary == 'no-slip' .or. boxarray(n)%right_boundary == 'free-slip' .or. &
          boxarray(n)%right_boundary == 'inflow') then
        u(n)%up(1) = boxarray(n)%nx
        us(n)%up(1) = boxarray(n)%nx
      endif

      v(n)%lo = 1
      vs(n)%lo = 1
      if (boxarray(n)%bottom_boundary == 'no-slip' .or. boxarray(n)%bottom_boundary == 'free-slip' .or. &
          boxarray(n)%bottom_boundary == 'inflow'  .or. boxarray(n)%bottom_boundary == 'periodic' .or. &
          boxarray(n)%bottom_boundary == 'internal' .or. boxarray(n)%bottom_boundary == 'halo') then
        v(n)%lo(2) = 2
        vs(n)%lo(2) = 2
      endif
      v(n)%up(1) = boxarray(n)%nx
      vs(n)%up(1) = boxarray(n)%nx
      v(n)%up(2) = boxarray(n)%ny + 1
      vs(n)%up(2) = boxarray(n)%ny + 1
      if (boxarray(n)%top_boundary == 'no-slip' .or. boxarray(n)%top_boundary == 'free-slip' .or. &
          boxarray(n)%top_boundary == 'inflow') then
        v(n)%up(2) = boxarray(n)%ny
        vs(n)%up(2) = boxarray(n)%ny
      endif

      ! Fields located at the center of the cell have same extension of the box
      p(n)%lo = 1
      p(n)%up(1) = boxarray(n)%nx
      p(n)%up(2) = boxarray(n)%ny
      phi(n)%lo = 1
      phi(n)%up(1) = boxarray(n)%nx
      phi(n)%up(2) = boxarray(n)%ny
      mu(n)%lo = 1
      mu(n)%up(1) = boxarray(n)%nx
      mu(n)%up(2) = boxarray(n)%ny
      rho(n)%lo = 1
      rho(n)%up(1) = boxarray(n)%nx
      rho(n)%up(2) = boxarray(n)%ny

      ! Base on the type of boundary condition on the domain
      ! select the boundary conditions for the fields
      if (boxarray(n)%left_boundary == 'no-slip' .or. boxarray(n)%left_boundary == 'inflow') then
        u(n)%left = 'dirichlet'
        v(n)%left = 'dirichlet'
        us(n)%left = 'dirichlet'
        vs(n)%left = 'dirichlet'
        p(n)%left = 'neumann'
        phi(n)%left = 'neumann'
        mu(n)%left = 'neumann'
        rho(n)%left = 'neumann'
      elseif (boxarray(n)%left_boundary == 'free-slip') then
        u(n)%left = 'dirichlet'
        v(n)%left = 'neumann'
        us(n)%left = 'dirichlet'
        vs(n)%left = 'neumann'
        p(n)%left = 'neumann'
        phi(n)%left = 'neumann'
        mu(n)%left = 'neumann'
        rho(n)%left = 'neumann'
      elseif (boxarray(n)%left_boundary == 'outflow') then
        u(n)%left = 'neumann'
        v(n)%left = 'dirichlet'
        us(n)%left = 'neumann'
        vs(n)%left = 'dirichlet'
        p(n)%left = 'dirichlet'
        phi(n)%left = 'dirichlet'
        mu(n)%left = 'neumann'
        rho(n)%left = 'neumann'
      elseif (boxarray(n)%left_boundary == 'periodic') then
        u(n)%left = 'periodic'
        v(n)%left = 'periodic'
        us(n)%left = 'periodic'
        vs(n)%left = 'periodic'
        p(n)%left = 'periodic'
        phi(n)%left = 'periodic'
        mu(n)%left = 'periodic'
        rho(n)%left = 'periodic'
      elseif (boxarray(n)%left_boundary == 'internal') then
        u(n)%left = 'internal'
        v(n)%left = 'internal'
        us(n)%left = 'internal'
        vs(n)%left = 'internal'
        p(n)%left = 'internal'
        phi(n)%left = 'internal'
        mu(n)%left = 'internal'
        rho(n)%left = 'internal'
      elseif (boxarray(n)%left_boundary == 'halo') then
        u(n)%left = 'halo'
        v(n)%left = 'halo'
        us(n)%left = 'halo'
        vs(n)%left = 'halo'
        p(n)%left = 'halo'
        phi(n)%left = 'halo'
        mu(n)%left = 'halo'
        rho(n)%left = 'halo'
      else
        print *, 'ERROR on left_boundary'
        stop
      endif
      if (boxarray(n)%left_boundary == 'inflow') then
        mu(n)%left = 'dirichlet'
        rho(n)%left = 'dirichlet'
      endif

      if (boxarray(n)%right_boundary == 'no-slip' .or. boxarray(n)%right_boundary == 'inflow') then
        u(n)%right = 'dirichlet'
        v(n)%right = 'dirichlet'
        us(n)%right = 'dirichlet'
        vs(n)%right = 'dirichlet'
        p(n)%right = 'neumann'
        phi(n)%right = 'neumann'
        mu(n)%right = 'neumann'
        rho(n)%right = 'neumann'
      elseif (boxarray(n)%right_boundary == 'free-slip') then
        u(n)%right = 'dirichlet'
        v(n)%right = 'neumann'
        us(n)%right = 'dirichlet'
        vs(n)%right = 'neumann'
        p(n)%right = 'neumann'
        phi(n)%right = 'neumann'
        mu(n)%right = 'neumann'
        rho(n)%right = 'neumann'
      elseif (boxarray(n)%right_boundary == 'outflow') then
        u(n)%right = 'neumann'
        v(n)%right = 'dirichlet'
        us(n)%right = 'neumann'
        vs(n)%right = 'dirichlet'
        p(n)%right = 'dirichlet'
        phi(n)%right = 'dirichlet'
        mu(n)%right = 'neumann'
        rho(n)%right = 'neumann'
      elseif (boxarray(n)%right_boundary == 'periodic') then
        u(n)%right = 'periodic'
        v(n)%right = 'periodic'
        us(n)%right = 'periodic'
        vs(n)%right = 'periodic'
        p(n)%right = 'periodic'
        phi(n)%right = 'periodic'
        mu(n)%right = 'periodic'
        rho(n)%right = 'periodic'
      elseif (boxarray(n)%right_boundary == 'internal') then
        u(n)%right = 'internal'
        v(n)%right = 'internal'
        us(n)%right = 'internal'
        vs(n)%right = 'internal'
        p(n)%right = 'internal'
        phi(n)%right = 'internal'
        mu(n)%right = 'internal'
        rho(n)%right = 'internal'
      elseif (boxarray(n)%right_boundary == 'halo') then
        u(n)%right = 'halo'
        v(n)%right = 'halo'
        us(n)%right = 'halo'
        vs(n)%right = 'halo'
        p(n)%right = 'halo'
        phi(n)%right = 'halo'
        mu(n)%right = 'halo'
        rho(n)%right = 'halo'
      else
        print *, 'ERROR on right_boundary'
        stop
      endif
      if (boxarray(n)%right == 'inflow') then
        mu(n)%right = 'dirichlet'
        rho(n)%right = 'dirichlet'
      endif
  
      if (boxarray(n)%top_boundary == 'no-slip' .or. boxarray(n)%top_boundary == 'inflow') then
        u(n)%top = 'dirichlet'
        v(n)%top = 'dirichlet'
        us(n)%top = 'dirichlet'
        vs(n)%top = 'dirichlet'
        p(n)%top = 'neumann'
        phi(n)%top = 'neumann'
        mu(n)%top = 'neumann'
        rho(n)%top = 'neumann'
      elseif (boxarray(n)%top_boundary == 'free-slip' ) then
        u(n)%top = 'neumann'
        v(n)%top = 'dirichlet'
        us(n)%top = 'neumann'
        vs(n)%top = 'dirichlet'
        p(n)%top = 'neumann'
        phi(n)%top = 'neumann'
        mu(n)%top = 'neumann'
        rho(n)%top = 'neumann'
      elseif (boxarray(n)%top_boundary == 'outflow') then
        u(n)%top = 'dirichlet'
        v(n)%top = 'neumann'
        us(n)%top = 'dirichlet'
        vs(n)%top = 'neumann'
        p(n)%top = 'dirichlet'
        phi(n)%top = 'dirichlet'
        mu(n)%top = 'neumann'
        rho(n)%top = 'neumann'
      elseif (boxarray(n)%top_boundary == 'periodic') then
        u(n)%top = 'periodic'
        v(n)%top = 'periodic'
        us(n)%top = 'periodic'
        vs(n)%top = 'periodic'
        p(n)%top = 'periodic'
        phi(n)%top = 'periodic'
        mu(n)%top = 'periodic'
        rho(n)%top = 'periodic'
      elseif (boxarray(n)%top_boundary == 'internal') then
        u(n)%top = 'internal'
        v(n)%top = 'internal'
        us(n)%top = 'internal'
        vs(n)%top = 'internal'
        p(n)%top = 'internal'
        phi(n)%top = 'internal'
        mu(n)%top = 'internal'
        rho(n)%top = 'internal'
      elseif (boxarray(n)%top_boundary == 'halo') then
        u(n)%top = 'halo'
        v(n)%top = 'halo'
        us(n)%top = 'halo'
        vs(n)%top = 'halo'
        p(n)%top = 'halo'
        phi(n)%top = 'halo'
        mu(n)%top = 'halo'
        rho(n)%top = 'halo'
      else
        print *, 'ERROR on top boundary'
        stop
      endif
      if (boxarray(n)%top == 'inflow') then
        mu(n)%top = 'dirichlet'
        rho(n)%top = 'dirichlet'
      endif
  
      if (boxarray(n)%bottom_boundary == 'no-slip' .or. boxarray(n)%bottom_boundary == 'inflow') then
        u(n)%bottom = 'dirichlet'
        v(n)%bottom = 'dirichlet'
        us(n)%bottom = 'dirichlet'
        vs(n)%bottom = 'dirichlet'
        p(n)%bottom = 'neumann'
        phi(n)%bottom = 'neumann'
        mu(n)%bottom = 'neumann'
        rho(n)%bottom = 'neumann'
      elseif (boxarray(n)%bottom_boundary == 'free-slip') then
        u(n)%bottom = 'neumann'
        v(n)%bottom = 'dirichlet'
        us(n)%bottom = 'neumann'
        vs(n)%bottom = 'dirichlet'
        p(n)%bottom = 'neumann'
        phi(n)%bottom = 'neumann'
        mu(n)%bottom = 'neumann'
        rho(n)%bottom = 'neumann'
      elseif (boxarray(n)%bottom_boundary == 'outflow') then
        u(n)%bottom = 'dirichlet'
        v(n)%bottom = 'neumann'
        us(n)%bottom = 'dirichlet'
        vs(n)%bottom = 'neumann'
        p(n)%bottom = 'dirichlet'
        phi(n)%bottom = 'dirichlet'      
        mu(n)%bottom = 'neumann'
        rho(n)%bottom = 'neumann'
      elseif (boxarray(n)%bottom_boundary == 'periodic') then
        u(n)%bottom = 'periodic'
        v(n)%bottom = 'periodic'
        us(n)%bottom = 'periodic'
        vs(n)%bottom = 'periodic'
        p(n)%bottom = 'periodic'
        phi(n)%bottom = 'periodic'
        mu(n)%bottom = 'periodic'
        rho(n)%bottom = 'periodic'
      elseif (boxarray(n)%bottom_boundary == 'internal') then
        u(n)%bottom = 'internal'
        v(n)%bottom = 'internal'
        us(n)%bottom = 'internal'
        vs(n)%bottom = 'internal'
        p(n)%bottom = 'internal'
        phi(n)%bottom = 'internal'
        mu(n)%bottom = 'internal'
        rho(n)%bottom = 'internal'
      elseif (boxarray(n)%bottom_boundary == 'halo') then
        u(n)%bottom = 'halo'
        v(n)%bottom = 'halo'
        us(n)%bottom = 'halo'
        vs(n)%bottom = 'halo'
        p(n)%bottom = 'halo'
        phi(n)%bottom = 'halo'
        mu(n)%bottom = 'halo'
        rho(n)%bottom = 'halo'
      else
        print *, 'ERROR on bottom boundary'
        stop
      endif
      if (boxarray(n)%bottom == 'inflow') then
        mu(n)%bottom = 'dirichlet'
        rho(n)%bottom = 'dirichlet'
      endif
  
      ! Allocate all the fields
      allocate(u(n)%f(u(n)%lo(1)-1:u(n)%up(1)+1,u(n)%lo(2)-1:u(n)%up(2)+1))
      allocate(v(n)%f(v(n)%lo(1)-1:v(n)%up(1)+1,v(n)%lo(2)-1:v(n)%up(2)+1))
      allocate(p(n)%f(p(n)%lo(1)-1:p(n)%up(1)+1,p(n)%lo(2)-1:p(n)%up(2)+1))
      allocate(us(n)%f(us(n)%lo(1)-1:us(n)%up(1)+1,us(n)%lo(2)-1:us(n)%up(2)+1))
      allocate(vs(n)%f(vs(n)%lo(1)-1:vs(n)%up(1)+1,vs(n)%lo(2)-1:vs(n)%up(2)+1))
      allocate(phi(n)%f(phi(n)%lo(1)-1:phi(n)%up(1)+1,phi(n)%lo(2)-1:phi(n)%up(2)+1))
      allocate(mu(n)%f(mu(n)%lo(1)-1:mu(n)%up(1)+1,mu(n)%lo(2)-1:mu(n)%up(2)+1))
      allocate(rho(n)%f(rho(n)%lo(1)-1:rho(n)%up(1)+1,rho(n)%lo(2)-1:rho(n)%up(2)+1))
  
      ! Allocate boundary values and set to zero
      allocate(u(n)%l(u(n)%lo(2)-1:u(n)%up(2)+1))
      allocate(u(n)%r(u(n)%lo(2)-1:u(n)%up(2)+1))
      allocate(u(n)%t(u(n)%lo(1)-1:u(n)%up(1)+1))
      allocate(u(n)%b(u(n)%lo(1)-1:u(n)%up(1)+1))
      allocate(v(n)%l(v(n)%lo(2)-1:v(n)%up(2)+1))
      allocate(v(n)%r(v(n)%lo(2)-1:v(n)%up(2)+1))
      allocate(v(n)%t(v(n)%lo(1)-1:v(n)%up(1)+1))
      allocate(v(n)%b(v(n)%lo(1)-1:v(n)%up(1)+1))
      allocate(p(n)%l(p(n)%lo(2)-1:p(n)%up(2)+1))
      allocate(p(n)%r(p(n)%lo(2)-1:p(n)%up(2)+1))
      allocate(p(n)%b(p(n)%lo(1)-1:p(n)%up(1)+1))
      allocate(p(n)%t(p(n)%lo(1)-1:p(n)%up(1)+1))
      allocate(mu(n)%l(mu(n)%lo(2)-1:mu(n)%up(2)+1))
      allocate(mu(n)%r(mu(n)%lo(2)-1:mu(n)%up(2)+1))
      allocate(mu(n)%b(mu(n)%lo(1)-1:mu(n)%up(1)+1))
      allocate(mu(n)%t(mu(n)%lo(1)-1:mu(n)%up(1)+1))
      allocate(rho(n)%l(rho(n)%lo(2)-1:rho(n)%up(2)+1))
      allocate(rho(n)%r(rho(n)%lo(2)-1:rho(n)%up(2)+1))
      allocate(rho(n)%b(rho(n)%lo(1)-1:rho(n)%up(1)+1))
      allocate(rho(n)%t(rho(n)%lo(1)-1:rho(n)%up(1)+1))
  
      allocate(us(n)%l(us(n)%lo(2)-1:us(n)%up(2)+1))
      allocate(us(n)%r(us(n)%lo(2)-1:us(n)%up(2)+1))
      allocate(us(n)%t(us(n)%lo(1)-1:us(n)%up(1)+1))
      allocate(us(n)%b(us(n)%lo(1)-1:us(n)%up(1)+1))
      allocate(vs(n)%l(vs(n)%lo(2)-1:vs(n)%up(2)+1))
      allocate(vs(n)%r(vs(n)%lo(2)-1:vs(n)%up(2)+1))
      allocate(vs(n)%t(vs(n)%lo(1)-1:vs(n)%up(1)+1))
      allocate(vs(n)%b(vs(n)%lo(1)-1:vs(n)%up(1)+1))
      allocate(phi(n)%l(phi(n)%lo(2)-1:phi(n)%up(2)+1))
      allocate(phi(n)%r(phi(n)%lo(2)-1:phi(n)%up(2)+1))
      allocate(phi(n)%b(phi(n)%lo(2)-1:phi(n)%up(1)+1))
      allocate(phi(n)%t(phi(n)%lo(2)-1:phi(n)%up(1)+1))

      u(n)%l = 0.0
      u(n)%r = 0.0
      u(n)%t = 0.0
      u(n)%b = 0.0
      v(n)%l = 0.0
      v(n)%r = 0.0
      v(n)%t = 0.0
      v(n)%b = 0.0
      p(n)%l = 0.0
      p(n)%r = 0.0
      p(n)%t = 0.0
      p(n)%b = 0.0
      mu(n)%l = 0.0
      mu(n)%r = 0.0
      mu(n)%t = 0.0
      mu(n)%b = 0.0
      rho(n)%l = 0.0
      rho(n)%r = 0.0
      rho(n)%t = 0.0
      rho(n)%b = 0.0

      us(n)%l = 0.0
      us(n)%r = 0.0
      us(n)%t = 0.0
      us(n)%b = 0.0
      vs(n)%l = 0.0
      vs(n)%r = 0.0
      vs(n)%t = 0.0
      vs(n)%b = 0.0
      phi(n)%l = 0.0
      phi(n)%r = 0.0
      phi(n)%t = 0.0
      phi(n)%b = 0.0

      ! Set all the field to zero
      u(n)%f = 0.0
      v(n)%f = 0.0
      p(n)%f = 0.0
      us(n)%f = 0.0
      vs(n)%f = 0.0
      phi(n)%f = 0.0

      ! Set density and viscosity 
      mu(n)%f(:,:) = viscosity
      rho(n)%f(:,:) = density

      ! Set field type
      u(n)%location = 1
      v(n)%location = 2
      p(n)%location = 0
      us(n)%location = 1
      vs(n)%location = 2
      phi(n)%location = 0
      mu(n)%location = 0
      rho(n)%location = 0
 
      ! Check for periodicity
      if (boxarray(n)%left_boundary == 'periodic' .and. boxarray(n)%right_boundary /= 'periodic') then
        print *, 'ERROR: left and right boundaries must be both periodic or not periodic'
        call destroy_ns_solver()
      elseif (boxarray(n)%left_boundary /= 'periodic' .and. boxarray(n)%right_boundary == 'periodic') then
        print *, 'ERROR: left and right boundaries must be both periodic or not periodic'
        call destroy_ns_solver()
      elseif (boxarray(n)%left_boundary == 'periodic' .and. boxarray(n)%right_boundary == 'periodic') then
        periodic(1) = boxarray(n)%nx
      endif

      if (boxarray(n)%top_boundary == 'periodic' .and. boxarray(n)%bottom_boundary /= 'periodic') then
        print *, 'ERROR: top and bottom boundaries must be both periodic or not periodic'
        call destroy_ns_solver()
      elseif (boxarray(n)%top_boundary /= 'periodic' .and. boxarray(n)%bottom_boundary == 'periodic') then
        print *, 'ERROR: top and bottom bouddaries must be both periodic or not periodic'
        call destroy_ns_solver()
      elseif (boxarray(n)%top_boundary == 'periodic' .and. boxarray(n)%bottom_boundary == 'periodic') then
        periodic(2) = boxarray(n)%ny
      endif

      ! Set the box BC for the poisson solver
      boxarray(n)%left = p(n)%left 
      boxarray(n)%right = p(n)%right 
      boxarray(n)%top = p(n)%top 
      boxarray(n)%bottom = p(n)%bottom 

    end do box_cycle 
 
    ! Select the functions for the Poisson solver
    if (poisson_solver_type == 'itr') then
      call init_hypre_solver()
      solve_poisson => solve_poisson_hypre
    elseif (poisson_solver_type == 'fft') then
      call init_fft_solver()
      solve_poisson => solve_poisson_fft
    else
      print *, 'ERROR: unkown Poisson solver'
      call destroy_ns_solver()
    endif

    ! Set the timestep
    dtconv = 0.3*boxarray(1)%delta/1.0
    dtvisc = 0.125*(boxarray(1)%delta**2 + boxarray(1)%delta**2)/viscosity
    dt = min(dtconv,dtvisc)
    if (solvemultiphase) then
      dtvisc = 0.125*(boxarray(1)%delta**2 + boxarray(1)%delta**2)/(max(mu0/rho0,mu1/rho1))
      dtsurf = sqrt(0.5*(rho0+rho1)*(boxarray(1)%delta**3)/(pi*sigma))
      dt = min(dtconv,dtvisc,dtsurf)
    endif

  end subroutine init_ns_solver
  !===============================================================================

  !===============================================================================
  subroutine solve()

    ! This is the master subroutine called from the main program

    ! Because we are using Adam-Bashfort we need to store the old RHS of the
    ! momentum equation: du_o and dv_o. rhs is the right-hand side of the
    ! poisson equation.

    implicit none

    integer :: istep, i, j, n
    type(myarray), dimension(nbox) :: du_o, dv_o, rhs    

    ! Allocate the old RHS of the momentum equation and the RHS of Poisson equation
    ! for each box
    do n = 1,nbox
      allocate(du_o(n)%f(u(n)%lo(1):u(n)%up(1),u(n)%lo(2):u(n)%up(2)))
      allocate(dv_o(n)%f(v(n)%lo(1):v(n)%up(1),v(n)%lo(2):v(n)%up(2)))
      allocate(rhs(n)%f(p(n)%lo(1):p(n)%up(1),p(n)%lo(2):p(n)%up(2)))
    end do

    ! Overwrtie the boundary on the projected field and projector 
    ! operator with that on velocity and pressure
    do n = 1,nbox
      us(n)%l = u(n)%l
      us(n)%r = u(n)%r
      us(n)%t = u(n)%t
      us(n)%b = u(n)%b
      vs(n)%l = v(n)%l
      vs(n)%r = v(n)%r
      vs(n)%t = v(n)%t
      vs(n)%b = v(n)%b
      phi(n)%l = p(n)%l
      phi(n)%r = p(n)%r
      phi(n)%t = p(n)%t
      phi(n)%b = p(n)%b
    end do

    ! First apply the boundary conditions on the initial fields
    ! and material properties
    call boundary(u)
    call boundary(v)
    call boundary(p)
    call boundary(mu)
    call boundary(rho)

    ! Compute the RHS of momentum equation of the initial velocity field
    do n = 1,nbox
      du_o(n)%f = 0.0
      dv_o(n)%f = 0.0
      call updates(u(n), v(n), du_o(n)%f, dv_o(n)%f, mu(n), rho(n), boxarray(n)%delta, n)
    end do

    open(newunit = log, file = 'log')

    ! We advance in time the solution
    istep = 0

    do while (istep <= nstep .and. t < Tmax)
      
      istep = istep + 1

      if (associated(event_i)) call event_i()

      t = t + dt
      if (pid == 0) write(log,*) 'istep: ', istep, 'dt: ', dt, 't: ', t

      if (solvemultiphase) then
        call advect_vof(u, v, dt)
        call update_material_properties(rho, mu)
      endif

      ! Compute the predicted velocity field
      call compute_predicted_velocity(du_o, dv_o)

      ! Then compute the RHS of the Poisson equation
      call poisson_rhs()

      do n = 1,nbox
        rhs(n)%f(:,:) = phi(n)%f(1:boxarray(n)%nx,1:boxarray(n)%ny)
      end do

      ! Then we solve the poisson equation for the projector operator
      if (solvemultiphase) then
        call destroy_hypre_solver
        call init_hypre_solver
        call update_coefficient_matrix(rho)
      end if
      call solve_poisson(rhs)

      ! We copy the array rhs in the two-dimensional array phi
      do n = 1,nbox
        do j = 1,boxarray(n)%ny
          do i = 1,boxarray(n)%nx
            phi(n)%f(i,j) = rhs(n)%f(i,j)
          end do
        enddo
      end do
      call boundary(phi)

      ! Then project the predicted velocity field to the divergence-free
      ! velocity field
      call project_velocity()

      ! and update the pressure
      do n = 1,nbox
        p(n)%f = p(n)%f + phi(n)%f
      end do
      call boundary(p)

      ! Check divergence of the velocity field and CFL
      call check(istep)

      ! Output
      if (associated(event_output)) call event_output()
    end do

    if (associated(event_end)) call event_end()

    ! Finilize the simulation
    call destroy_ns_solver()

  end subroutine solve
  !===============================================================================

  !===============================================================================
  subroutine updates(un, vn, du, dv, mun, rhon, delta, n)

    implicit none

    ! Input / Output variables
    integer, intent(in) :: n
    type(field), intent(in) :: un, vn, mun, rhon
    real(kind=dp), intent(inout) :: du(un%lo(1):un%up(1),un%lo(2):un%up(2))
    real(kind=dp), intent(inout) :: dv(vn%lo(1):vn%up(1),vn%lo(2):vn%up(2))
    real(kind=dp), intent(in) :: delta

    ! Local variables
    integer :: i, j, ip, im, jp, jm
    real(kind=dp) :: uuip, uuim, uvjp, uvjm, tauxxij, tauxxim, dtauxxdx
    real(kind=dp) :: tauxyjp, tauxyij, dtauxydy, uvip, uvim, vvjp, vvjm
    real(kind=dp) :: tauyyij, tauyyjm, dtauyydy, tauxyip, dtauxydx, muc
    real(kind=dp) :: delta2
    real(kind=dp) :: fsigmax(un%lo(1):un%up(1),un%lo(2):un%up(2))
    real(kind=dp) :: fsigmay(vn%lo(1):vn%up(1),vn%lo(2):vn%up(2))

    delta2 = delta**2

    fsigmax = 0.
    fsigmay = 0.
    if (solvemultiphase) call compute_surface_tension_force(fsigmax,    &
                              [un%lo(1),un%lo(2)], [un%up(1),un%up(2)], &
                              fsigmay, [vn%lo(1),vn%lo(2)],             &
                              [vn%up(1),vn%up(2)], n)

    do j = un%lo(2),un%up(2)
      jp = j + 1
      jm = j - 1
      do i = un%lo(1),un%up(1)
        ip = i + 1
        im = i - 1

        ! Advection terms
        uuip  = 0.25 * ( un%f(ip,j) + un%f(i,j) ) * ( un%f(ip,j) + un%f(i,j)  )
        uuim  = 0.25 * ( un%f(im,j) + un%f(i,j) ) * ( un%f(im,j) + un%f(i,j)  )
        uvjp  = 0.25 * ( un%f(i,jp) + un%f(i,j) ) * ( vn%f(i,jp) + vn%f(im,jp))
        uvjm  = 0.25 * ( un%f(i,jm) + un%f(i,j) ) * ( vn%f(i,j) + vn%f(im,j)  )

        ! Viscous term: xx and yy terms are at cell center, xy ar at the corner
        tauxxij = 2.*mun%f(i,j)*(un%f(ip,j) - un%f(i,j))/delta
        tauxxim = 2.*mun%f(im,j)*(un%f(i,j) - un%f(im,j))/delta
        dtauxxdx = (tauxxij - tauxxim)/delta

        muc = 0.25*(mun%f(im,j) + mun%f(im,jp) + mun%f(i,j) + mun%f(i,jp))
        tauxyjp = muc*((un%f(i,jp) - un%f(i,j))/delta + (vn%f(i,jp) - vn%f(im,jp))/delta)
        muc = 0.25*(mun%f(im,jm) + mun%f(im,j) + mun%f(i,jm) + mun%f(i,j))
        tauxyij = muc*((un%f(i,j) - un%f(i,jm))/delta + (vn%f(i,j) - vn%f(im,j))/delta)
        dtauxydy = (tauxyjp - tauxyij) / delta

        du(i,j) = (-uuip + uuim) / delta + (-uvjp + uvjm) / delta + &
          (dtauxxdx + dtauxydy + fsigmax(i,j))/(0.5*(rhon%f(i,j) + rhon%f(im,j)))
            !mun%f(i,j)*((un%f(ip,j) - 2.0*un%f(i,j) + un%f(im,j))/delta2 + &
            !(un%f(i,jp) - 2.0*un%f(i,j) + un%f(i,jm))/delta2)/ & 
            !(0.5*(rhon%f(i,j) + rhon%f(i-1,j)))
      end do
    end do

    do j = vn%lo(2),vn%up(2)
      jp = j + 1
      jm = j - 1
      do i = vn%lo(1),vn%up(1)
        ip = i + 1
        im = i - 1

        ! Advection terms
        uvip  = 0.25 * ( un%f(ip,jm) + un%f(ip,j) ) * ( vn%f(i,j) + vn%f(ip,j) )
        uvim  = 0.25 * ( un%f(i,jm) + un%f(i,j) ) * ( vn%f(i,j) + vn%f(im,j) )
        vvjp  = 0.25 * ( vn%f(i,j) + vn%f(i,jp) ) * ( vn%f(i,j) + vn%f(i,jp) )
        vvjm  = 0.25 * ( vn%f(i,j) + vn%f(i,jm) ) * ( vn%f(i,j) + vn%f(i,jm) )

        ! Viscous terms
        tauyyij = 2.*mun%f(i,j)*(vn%f(i,jp) - vn%f(i,j))/delta
        tauyyjm = 2.*mun%f(i,jm)*(vn%f(i,j) - vn%f(i,jm))/delta
        dtauyydy = (tauyyij - tauyyjm)/delta

        muc = 0.25*(mun%f(i,jm) + mun%f(i,j) + mun%f(ip,jm) + mun%f(ip,j))
        tauxyip = muc*((un%f(ip,j) - un%f(ip,jm))/delta + (vn%f(ip,j) - vn%f(i,j))/delta)
        muc = 0.25*(mun%f(im,jm) + mun%f(im,j) + mun%f(i,jm) + mun%f(i,j))
        tauxyij = muc*((un%f(i,j) - un%f(i,jm))/delta + (vn%f(i,j) - vn%f(im,j))/delta)
        dtauxydx = (tauxyip - tauxyij) / delta

        dv(i,j) = (-uvip + uvim) / delta + (-vvjp + vvjm) / delta + &
          (dtauyydy + dtauxydx + fsigmay(i,j))/(0.5*(rhon%f(i,j) + rhon%f(i,jm)))
          !mun%f(i,j)*((vn%f(ip,j) - 2.0*vn%f(i,j) + vn%f(im,j))/delta2 + &
          !(vn%f(i,jp) - 2.0*vn%f(i,j) + vn%f(i,jm))/delta2)/ &
          !(0.5*(rhon%f(i,j) + rhon%f(i,j-1)))
      end do
    end do

  end subroutine updates
  !===============================================================================

  !===============================================================================
  subroutine compute_predicted_velocity(du_o, dv_o)

    implicit none

    ! Input / Output variables
    type(myarray), dimension(nbox) :: du_o, dv_o

    ! Local variables
    integer :: i, j, n
    real(kind=dp) :: rhs, delta
    real(kind=dp), dimension(:,:), allocatable :: du, dv

    do n = 1,nbox

      delta = boxarray(n)%delta

      allocate(du(u(n)%lo(1):u(n)%up(1),u(n)%lo(2):u(n)%up(2)))
      allocate(dv(v(n)%lo(1):v(n)%up(1),v(n)%lo(2):v(n)%up(2)))

      ! Compute the new RHS of momentum equation
      call updates(u(n), v(n), du, dv, mu(n), rho(n), boxarray(n)%delta, n)

      ! Compute the predicted velocity field
      ! x direction
      do j = us(n)%lo(2),us(n)%up(2)
        do i = us(n)%lo(1),us(n)%up(1)
          rhs = 1.5*du(i,j) - 0.5*du_o(n)%f(i,j) - (p(n)%f(i,j) - p(n)%f(i-1,j)) / &
                (0.5*(rho(n)%f(i-1,j) + rho(n)%f(i,j))*delta) + Sx + g(1)
          us(n)%f(i,j) = u(n)%f(i,j) + dt*rhs
        end do
      end do

      ! y direction
      do j = vs(n)%lo(2),vs(n)%up(2)
        do i = vs(n)%lo(1),vs(n)%up(1)
          rhs = 1.5*dv(i,j) - 0.5*dv_o(n)%f(i,j) - (p(n)%f(i,j) - p(n)%f(i,j-1)) / &
                (0.5*(rho(n)%f(i,j-1) + rho(n)%f(i,j))*delta) + Sy + g(2)
          vs(n)%f(i,j) = v(n)%f(i,j) + dt*rhs
        end do
      end do

      du_o(n)%f = du
      dv_o(n)%f = dv

      deallocate(du,dv)

    end do

    call boundary(us)
    call boundary(vs)

  end  subroutine compute_predicted_velocity
  !===============================================================================

  !===============================================================================
  subroutine poisson_rhs()

    implicit none

    ! Local variables
    integer :: i, j, n
    real(kind=dp) :: delta

    do n = 1,nbox
      delta = boxarray(n)%delta
      do j = phi(n)%lo(2),phi(n)%up(2)
        do i = phi(n)%lo(1),phi(n)%up(1)
          phi(n)%f(i,j) = (1./dt)*((us(n)%f(i+1,j) - us(n)%f(i,j))/delta + &
                                   (vs(n)%f(i,j+1) - vs(n)%f(i,j))/delta)
        end do
      end do
    end do

  end subroutine poisson_rhs
  !===============================================================================

  !===============================================================================
  subroutine project_velocity()

    implicit none

    ! local variables
    integer :: i, j, n
    real(kind=dp) :: delta

    do n = 1,nbox
      delta = boxarray(n)%delta
      ! x direction
      do j = u(n)%lo(2),u(n)%up(2)
        do i = u(n)%lo(1),u(n)%up(1)
          u(n)%f(i,j) = us(n)%f(i,j) - dt*(phi(n)%f(i,j) - phi(n)%f(i-1,j))/ &
                        (0.5*(rho(n)%f(i-1,j) + rho(n)%f(i,j))*delta)
        end do
      end do

      ! y direction
      do j = v(n)%lo(2),v(n)%up(2)
        do i = v(n)%lo(1),v(n)%up(1)
          v(n)%f(i,j) = vs(n)%f(i,j) - dt*(phi(n)%f(i,j) - phi(n)%f(i,j-1))/ &
                        (0.5*(rho(n)%f(i,j-1) + rho(n)%f(i,j))*delta)
        end do
      end do
    end do

    call boundary(u)
    call boundary(v)

  end subroutine project_velocity
  !===============================================================================

  !===============================================================================
  subroutine check(istep)

    ! This subroutine check the divergence of the velocity field
    ! and the instant CFL of the simulation 

    implicit none
    integer, intent(in) :: istep

    integer :: i, j, n, imax, jmax
    real(kind=dp) :: divergence, maximum_divergence, CFL, maximum_CFL

    maximum_divergence = 0.0
    maximum_CFL = 0.0

    divergence_tol = tolerance

    do n = 1,nbox
      do j = 1,boxarray(n)%ny
        do i = 1,boxarray(n)%nx
          ! Check divergence
          divergence = (u(n)%f(i+1,j) - u(n)%f(i,j)) / boxarray(n)%delta + &
                       (v(n)%f(i,j+1) - v(n)%f(i,j)) / boxarray(n)%delta
          if (divergence > maximum_divergence) then
            imax = i
            jmax = j
            maximum_divergence = divergence
          end if

          ! Check CFL
          CFL = 0.5*dt*( (u(n)%f(i+1,j) + u(n)%f(i,j)) / boxarray(n)%delta + &
                         (v(n)%f(i,j+1) + v(n)%f(i,j)) / boxarray(n)%delta )
          if (CFL > maximum_CFL) maximum_CFL = CFL
        end do
      end do
    end do

    call mpi_allreduce(mpi_in_place,maximum_divergence,1,mpi_real8,mpi_max, &
                       mpi_common_world,ierr)
    call mpi_allreduce(mpi_in_place,maximum_CFL,1,mpi_real8,mpi_max, &
                       mpi_common_world,ierr)

    if (maximum_divergence > divergence_tol) then
      if (pid == 0) write(log,*) 'WARNING: divergence is ', maximum_divergence, 'at istep:', istep
      if (pid == 0) write(log,*) 'in ', imax, jmax
    endif

    if (maximum_CFL > 0.8) then
      if (pid == 0) write(log,*) 'WARNING: CFL is:', maximum_CFL 
    endif

  end subroutine check
  !===============================================================================

  !===============================================================================
  subroutine destroy_ns_solver()

    implicit none
 
    integer :: n

    ! Free the memory and finalize the simulation
    do n = 1,nbox
      deallocate(p(n)%f, phi(n)%f, u(n)%f, us(n)%f, v(n)%f, vs(n)%f)
      deallocate(u(n)%l, u(n)%b, u(n)%t, u(n)%r)
      deallocate(v(n)%l, v(n)%b, v(n)%t, v(n)%r)
      deallocate(p(n)%l, p(n)%b, p(n)%t, p(n)%r)
      deallocate(us(n)%l, us(n)%b, us(n)%t, us(n)%r)
      deallocate(vs(n)%l, vs(n)%b, vs(n)%t, vs(n)%r)
      deallocate(phi(n)%l, phi(n)%b, phi(n)%t, phi(n)%r)
    end do
    if (poisson_solver_type == 'itr') then
      call destroy_hypre_solver
    elseif (poisson_solver_type == 'fft') then
      call destroy_fft_solver
    endif
    if (solvemultiphase) call destroy_vof
    call destroy_boxes
    call MPI_finalize(ierr)
    stop

  end subroutine destroy_ns_solver
  !===============================================================================

end module navier_stokes
