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
  use navier_stokes_pub
  include 'mpif.h'

  type(field) :: phi

  private
  public :: init_ns_solver, solve, destroy_ns_solver

contains

  !===============================================================================
  subroutine init_ns_solver

#ifdef IBM
    use ibm, only : ibm_init
#endif
    
    implicit none

    ! Initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    mpi_common_world = MPI_COMM_WORLD

    ! Based on the type of boundary, compute the extension for each field
    u%lo = 1
    us%lo = 1
    if (left_boundary == 'no-slip' .or. left_boundary == 'free-slip' .or. &
        left_boundary == 'periodic') then
      u%lo(1) = 2
      us%lo(1) = 2
    endif
    u%up(1) = nx + 1
    us%up(1) = nx + 1
    u%up(2) = ny
    us%up(2) = ny
    if (right_boundary == 'no-slip' .or. right_boundary == 'free-slip') then
      u%up(1) = nx
      us%up(1) = nx
    endif

    v%lo = 1
    vs%lo = 1
    if (bottom_boundary == 'no-slip' .or. bottom_boundary == 'free-slip' .or. &
        bottom_boundary == 'periodic') then
      v%lo(2) = 2
      vs%lo(2) = 2
    endif
    v%up(1) = nx
    vs%up(1) = nx
    v%up(2) = ny + 1
    vs%up(2) = ny + 1
    if (top_boundary == 'no-slip' .or. top_boundary == 'free-slip') then
      v%up(2) = ny
      vs%up(2) = ny
    endif
    
    p%lo = 1
    p%up(1) = nx
    p%up(2) = ny
    
    ! Allocate all the fields
    allocate(p%f(p%lo(1)-1:p%up(1)+1,p%lo(2)-1:p%up(2)+1))
    allocate(phi%f(phi%lo(1)-1:phi%up(1)+1,phi%lo(2)-1:phi%up(2)+1))
    allocate(u%f(u%lo(1)-1:u%up(1)+1,u%lo(2)-1:u%up(2)+1))
    allocate(us%f(us%lo(1)-1:us%up(1)+1,us%lo(2)-1:us%up(2)+1))
    allocate(v%f(v%lo(1)-1:v%up(1)+1,v%lo(2)-1:v%up(2)+1))
    allocate(vs%f(vs%lo(1)-1:vs%up(1)+1,vs%lo(2)-1:vs%up(2)+1))

    ! Allocate boundary values
    allocate(u%l(u%up(2)))
    allocate(u%r(u%up(2)))
    allocate(u%t(u%up(1)))
    allocate(u%b(u%up(1)))
    allocate(v%l(v%up(2)))
    allocate(v%r(v%up(2)))
    allocate(v%t(v%up(1)))
    allocate(v%b(v%up(1)))
    allocate(p%l(p%up(2)))
    allocate(p%r(p%up(2)))
    allocate(p%b(p%up(1)))
    allocate(p%t(p%up(1)))

    ! We set all the field to zero
    p%f = 0.0
    phi%f = 0.0
    u%f = 0.0
    us%f = 0.0
    v%f = 0.0
    vs%f = 0.0

    ! Set field type
    u%location = 1
    us%location = 1
    v%location = 2
    vs%location = 2
    p%location = 0
    phi%location = 0

    ! Check for periodicity
    if (left_boundary == 'periodic' .and. right_boundary /= 'periodic') then
      print *, 'ERROR: left and right boundaries must be both periodic or not periodic'
      call destroy_ns_solver()
    elseif (left_boundary /= 'periodic' .and. right_boundary == 'periodic') then
      print *, 'ERROR: left and right boundaries must be both periodic or not periodic'
      call destroy_ns_solver()
    elseif (left_boundary == 'periodic' .and. right_boundary == 'periodic') then
      periodic(1) = nx
    endif

    if (top_boundary == 'periodic' .and. bottom_boundary /= 'periodic') then
      print *, 'ERROR: top and bottom boundaries must be both periodic or not periodic'
      call destroy_ns_solver()
    elseif (top_boundary /= 'periodic' .and. bottom_boundary == 'periodic') then
      print *, 'ERROR: top and bottom bouddaries must be both periodic or not periodic'
      call destroy_ns_solver()
    elseif (top_boundary == 'periodic' .and. bottom_boundary == 'periodic') then
      periodic(2) = ny
    endif

    ! Create the hypre solver for the poisson equation
    call init_hypre_solver

    ! If using IBM initialize solid bodies
#ifdef IBM
    call ibm_init()
#endif

  end subroutine init_ns_solver
  !===============================================================================

  !===============================================================================
  subroutine solve()

    ! This is the master subroutine called from the main program

    ! Because we are using Adam-Bashfort we need to store the old RHS of the
    ! momentum equation: du_o and dv_o. rhs is the right-hand side of the
    ! poisson equation. Is a 1D-array of dimension nx*ny.

#ifdef IBM
    use ibm, only : ibm_tag
#endif
    
    implicit none

    integer :: istep, i, j
    real :: du_o(nx+1,ny), dv_o(nx,ny+1), rhs(nx*ny)

    ! First apply the boundary conditions on the initial fields
    call boundary_u(u)
    call boundary_v(v)
    call boundary_p(p)

    ! Compute the RHS of momentum equation of the initial velocity field
    du_o = 0.0
    dv_o = 0.0
    call updates(du_o, dv_o)

    open(newunit = log, file = 'log')

    ! If using IBM initialize solid bodies
#ifdef IBM
    call ibm_tag()
#endif
    
    ! We advance in time the solution
    do istep = 1, nstep

      if (associated(event_i)) call event_i()

      t = t + dt
      write(log,*) 'istep: ', istep, 'dt: ', dt, 't: ', t

      
      ! First we compute the predicted velocity field
      call compute_predicted_velocity(du_o, dv_o)

      ! We compute the RHS of the poisson equation and store it in the
      ! array rhs
      call poisson_rhs(rhs)

      ! Then we solve the poisson equation for the projector operator
      call solve_poisson(rhs)

      ! We copy the array rhs in the two-dimensional array phi
      do j = 1,ny
        do i = 1,nx
          phi%f(i,j) = rhs(i + (j-1)*ny)
        end do
      end do
      call boundary_p(phi)

      ! Then we project the predicted velocity field to the divergence-free
      ! velocity
      call project_velocity()

      ! We update the pressure
      p%f = p%f + phi%f
      call boundary_p(p)

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
  subroutine updates(du, dv)

    real, intent(inout) :: du(:,:), dv(:,:)

    integer :: i, j, ip, im, jp, jm
    real :: uuip, uuim, uvjp, uvjm, uvip, uvim, vvjp, vvjm
    real :: dxq, dyq

    dxq = dx*dx
    dyq = dy*dy

    do j = 2,ny
      jp = j + 1
      jm = j - 1
      do i = 2,nx
        ip = i + 1
        im = i - 1
        ! x direction
        uuip  = 0.25 * ( u%f(ip,j) + u%f(i,j) ) * ( u%f(ip,j) + u%f(i,j)  )
        uuim  = 0.25 * ( u%f(im,j) + u%f(i,j) ) * ( u%f(im,j) + u%f(i,j)  )
        uvjp  = 0.25 * ( u%f(i,jp) + u%f(i,j) ) * ( v%f(i,jp) + v%f(im,jp))
        uvjm  = 0.25 * ( u%f(i,jm) + u%f(i,j) ) * ( v%f(i,j) + v%f(im,j)  )
        du(i,j) = (-uuip + uuim) / dx + (-uvjp + uvjm) / dy + &
             mu*((u%f(ip,j) - 2.0*u%f(i,j) + u%f(im,j))/dxq + &
             (u%f(i,jp) - 2.0*u%f(i,j) + u%f(i,jm))/dyq)

        ! y direction
        uvip  = 0.25 * ( u%f(ip,jm) + u%f(ip,j) ) * ( v%f(i,j) + v%f(ip,j) )
        uvim  = 0.25 * ( u%f(i,jm) + u%f(i,j) ) * ( v%f(i,j) + v%f(im,j) )
        vvjp  = 0.25 * ( v%f(i,j) + v%f(i,jp) ) * ( v%f(i,j) + v%f(i,jp) )
        vvjm  = 0.25 * ( v%f(i,j) + v%f(i,jm) ) * ( v%f(i,j) + v%f(i,jm) )
        dv(i,j) = (-uvip + uvim) / dx + (-vvjp + vvjm) / dy + &
             mu*((v%f(ip,j) - 2.0*v%f(i,j) + v%f(im,j))/dxq + &
             (v%f(i,jp) - 2.0*v%f(i,j) + v%f(i,jm))/dyq)      
      end do
    end do

    ! Bottom row of the grid for u
    j = 1
    jp = j + 1
    jm = j - 1
    do i = 2,nx
      ip = i + 1
      im = i - 1
      uuip  = 0.25 * ( u%f(ip,j) + u%f(i,j) ) * ( u%f(ip,j) + u%f(i,j)  )
      uuim  = 0.25 * ( u%f(im,j) + u%f(i,j) ) * ( u%f(im,j) + u%f(i,j)  )
      uvjp  = 0.25 * ( u%f(i,jp) + u%f(i,j) ) * ( v%f(i,jp) + v%f(im,jp))
      uvjm  = 0.25 * ( u%f(i,jm) + u%f(i,j) ) * ( v%f(i,j) + v%f(im,j)  )
      du(i,j) = (-uuip + uuim) / dx + (-uvjp + uvjm) / dy + &
           mu*((u%f(ip,j) - 2.0*u%f(i,j) + u%f(im,j))/dxq + &
           (u%f(i,jp) - 2.0*u%f(i,j) + u%f(i,jm))/dyq)
    end do

    ! Left column of the grid for v
    i = 1
    ip = i + 1
    im = i - 1
    do j = 2,ny
      jp = j + 1
      jm = j - 1
      uvip  = 0.25 * ( u%f(ip,jm) + u%f(ip,j) ) * ( v%f(i,j) + v%f(ip,j) )
      uvim  = 0.25 * ( u%f(i,jm) + u%f(i,j) ) * ( v%f(i,j) + v%f(im,j) )
      vvjp  = 0.25 * ( v%f(i,j) + v%f(i,jp) ) * ( v%f(i,j) + v%f(i,jp) )
      vvjm  = 0.25 * ( v%f(i,j) + v%f(i,jm) ) * ( v%f(i,j) + v%f(i,jm) )
      dv(i,j) = (-uvip + uvim) / dx + (-vvjp + vvjm) / dy + &
           mu*((v%f(ip,j) - 2.0*v%f(i,j) + v%f(im,j))/dxq + &
           (v%f(i,jp) - 2.0*v%f(i,j) + v%f(i,jm))/dyq)
    end do

    ! If the domain is periodic in x direction we need to solve
    ! also u on the left and right boundaries
    if (left_boundary == 'periodic') then
      i = 1
      ip = i + 1
      im = i - 1
      do j = 1,ny
        jp = j + 1
        jm = j - 1
        uuip  = 0.25 * ( u%f(ip,j) + u%f(i,j) ) * ( u%f(ip,j) + u%f(i,j)  )
        uuim  = 0.25 * ( u%f(im,j) + u%f(i,j) ) * ( u%f(im,j) + u%f(i,j)  )
        uvjp  = 0.25 * ( u%f(i,jp) + u%f(i,j) ) * ( v%f(i,jp) + v%f(im,jp))
        uvjm  = 0.25 * ( u%f(i,jm) + u%f(i,j) ) * ( v%f(i,j) + v%f(im,j)  )
        du(i,j) = (-uuip + uuim) / dx + (-uvjp + uvjm) / dy + &
             mu*((u%f(ip,j) - 2.0*u%f(i,j) + u%f(im,j))/dxq + &
             (u%f(i,jp) - 2.0*u%f(i,j) + u%f(i,jm))/dyq)
      end do

      i = nx + 1
      ip = i + 1
      im = i - 1
      do j = 1,ny
        jp = j + 1
        jm = j - 1
        uuip  = 0.25 * ( u%f(ip,j) + u%f(i,j) ) * ( u%f(ip,j) + u%f(i,j)  )
        uuim  = 0.25 * ( u%f(im,j) + u%f(i,j) ) * ( u%f(im,j) + u%f(i,j)  )
        uvjp  = 0.25 * ( u%f(i,jp) + u%f(i,j) ) * ( v%f(i,jp) + v%f(im,jp))
        uvjm  = 0.25 * ( u%f(i,jm) + u%f(i,j) ) * ( v%f(i,j) + v%f(im,j)  )
        du(i,j) = (-uuip + uuim) / dx + (-uvjp + uvjm) / dy + &
             mu*((u%f(ip,j) - 2.0*u%f(i,j) + u%f(im,j))/dxq + &
             (u%f(i,jp) - 2.0*u%f(i,j) + u%f(i,jm))/dyq)
      end do
    end if

    ! If the domain is periodic in y direction we need to solve
    ! also v on the bottom and top boundaries
    if (top_boundary == 'periodic') then
      j = 1
      jp = j + 1
      jm = j - 1
      do i = 1,nx
        ip = i + 1
        im = i - 1
        uvip  = 0.25 * ( u%f(ip,jm) + u%f(ip,j) ) * ( v%f(i,j) + v%f(ip,j) )
        uvim  = 0.25 * ( u%f(i,jm) + u%f(i,j) ) * ( v%f(i,j) + v%f(im,j) )
        vvjp  = 0.25 * ( v%f(i,j) + v%f(i,jp) ) * ( v%f(i,j) + v%f(i,jp) )
        vvjm  = 0.25 * ( v%f(i,j) + v%f(i,jm) ) * ( v%f(i,j) + v%f(i,jm) )
        dv(i,j) = (-uvip + uvim) / dx + (-vvjp + vvjm) / dy + &
             mu*((v%f(ip,j) - 2.0*v%f(i,j) + v%f(im,j))/dxq + &
             (v%f(i,jp) - 2.0*v%f(i,j) + v%f(i,jm))/dyq)      
      end do

      j = ny + 1
      jp = j + 1
      jm = j - 1
      do i = 1,nx
        ip = i + 1
        im = i - 1
        uvip  = 0.25 * ( u%f(ip,jm) + u%f(ip,j) ) * ( v%f(i,j) + v%f(ip,j) )
        uvim  = 0.25 * ( u%f(i,jm) + u%f(i,j) ) * ( v%f(i,j) + v%f(im,j) )
        vvjp  = 0.25 * ( v%f(i,j) + v%f(i,jp) ) * ( v%f(i,j) + v%f(i,jp) )
        vvjm  = 0.25 * ( v%f(i,j) + v%f(i,jm) ) * ( v%f(i,j) + v%f(i,jm) )
        dv(i,j) = (-uvip + uvim) / dx + (-vvjp + vvjm) / dy + &
             mu*((v%f(ip,j) - 2.0*v%f(i,j) + v%f(im,j))/dxq + &
             (v%f(i,jp) - 2.0*v%f(i,j) + v%f(i,jm))/dyq)      
      end do
    endif

  end subroutine updates
  !===============================================================================

  !===============================================================================
  subroutine compute_predicted_velocity(du_o, dv_o)

#ifdef IBM
    use ibm, only : ibm_force
#endif
 
    implicit none
    real, intent(inout) :: du_o(:,:), dv_o(:,:)

    integer :: i, j
    real :: du(nx+1,ny), dv(nx,ny+1), rhs

    ! Compute the new RHS of momentum equation
    call updates(du, dv)

    ! Compute the predicted velocity field
    do j = 2,ny
      do i = 2,nx

        ! x direction
        rhs = 1.5 * du(i,j) - 0.5 * du_o(i,j) - (p%f(i,j) - p%f(i-1,j)) / (rho*dx) + Sx
        us%f(i,j) = u%f(i,j) + dt * rhs

#ifdef IBM
        us%f(i,j) = us%f(i,j) + dt * ibm_force(i,j,1,rhs,u%f(i,j),dt)
#endif

        ! y direction
        rhs = 1.5 * dv(i,j) - 0.5 * dv_o(i,j) - (p%f(i,j) - p%f(i,j-1)) / (rho*dy) + Sy
        vs%f(i,j) = v%f(i,j) + dt * rhs

#ifdef IBM
        vs%f(i,j) = vs%f(i,j) + dt * ibm_force(i,j,2,rhs,v%f(i,j),dt)  
#endif

      end do
    end do

    ! bottom row of the grid for us
    j = 1
    do i = 2,nx
      rhs = 1.5 * du(i,j) - 0.5 * du_o(i,j) - (p%f(i,j) - p%f(i-1,j)) / (rho*dx) + Sx
      us%f(i,j) = u%f(i,j) + dt * rhs

#ifdef IBM
        us%f(i,j) = us%f(i,j) + dt * ibm_force(i,j,1,rhs,u%f(i,j),dt)
#endif

    end do

    ! left column of the grid for vs
    i = 1
    do j = 2,ny
      rhs = 1.5 * dv(i,j) - 0.5 * dv_o(i,j) - (p%f(i,j) - p%f(i,j-1)) / (rho*dy) + Sy
      vs%f(i,j) = v%f(i,j) + dt * rhs

#ifdef IBM
        vs%f(i,j) = vs%f(i,j) + dt * ibm_force(i,j,2,rhs,v%f(i,j),dt)  
#endif

    end do

    ! If the domain is periodic in x direction we need to solve
    ! also us on the left and right boundaries
    if (left_boundary == 'periodic') then
      i = 1
      do j = 1,ny
        rhs = 1.5 * du(i,j) - 0.5 * du_o(i,j) - (p%f(i,j) - p%f(i-1,j)) / (rho*dx) + Sx
        us%f(i,j) = u%f(i,j) + dt * rhs

#ifdef IBM
        us%f(i,j) = us%f(i,j) + dt * ibm_force(i,j,1,rhs,u%f(i,j),dt)
#endif

      end do

      i = nx + 1
      do j = 1,ny
        rhs = 1.5 * du(i,j) - 0.5 * du_o(i,j) - (p%f(i,j) - p%f(i-1,j)) / (rho*dx) + Sx
        us%f(i,j) = u%f(i,j) + dt * rhs

#ifdef IBM
        us%f(i,j) = us%f(i,j) + dt * ibm_force(i,j,1,rhs,u%f(i,j),dt)
#endif

      end do
    end if

    ! If the domain is periodic in y direction we need to solve
    ! also v on the bottom and top boundaries
    if (top_boundary == 'periodic') then
      j = 1
      do i = 1,nx
        rhs = 1.5 * dv(i,j) - 0.5 * dv_o(i,j) - (p%f(i,j) - p%f(i,j-1)) / (rho*dy) + Sy
        vs%f(i,j) = v%f(i,j) + dt * rhs

#ifdef IBM
        vs%f(i,j) = vs%f(i,j) + dt * ibm_force(i,j,2,rhs,v%f(i,j),dt)
#endif

      end do

      j = ny + 1
      do i = 1,nx
        rhs = 1.5 * dv(i,j) - 0.5 * dv_o(i,j) - (p%f(i,j) - p%f(i,j-1)) / (rho*dy) + Sy
        vs%f(i,j) = v%f(i,j) + dt * rhs

#ifdef IBM
        vs%f(i,j) = vs%f(i,j) + dt * ibm_force(i,j,2,rhs,v%f(i,j),dt) 
#endif

      end do
    endif

    du_o = du
    dv_o = dv

    call boundary_u(us)
    call boundary_v(vs)

  end  subroutine compute_predicted_velocity
  !===============================================================================

  !===============================================================================
  subroutine poisson_rhs(rhs)

    implicit none

    real, dimension(:), intent(inout) :: rhs

    integer :: i, j

    do j = 1,ny
      do i = 1,nx
        rhs(i+(j-1)*ny) = (rho/dt) * ((us%f(i+1,j) - us%f(i,j)) / dx + &
             (vs%f(i,j+1) - vs%f(i,j)) / dy )
      end do
    end do

  end subroutine poisson_rhs
  !===============================================================================

  !===============================================================================
  subroutine project_velocity

    ! local variables
    integer :: i, j

    do j = 2,ny
      do i = 2,nx
        ! x direction
        u%f(i,j) = us%f(i,j) - (dt / rho) * ( phi%f(i,j) - phi%f(i-1,j) ) / dx

        ! y direction
        v%f(i,j) = vs%f(i,j) - (dt / rho) * ( phi%f(i,j) - phi%f(i,j-1) ) / dy
      end do
    end do

    ! Bottom row for u component
    j = 1
    do i = 2,nx
      u%f(i,j) = us%f(i,j) - (dt / rho) * ( phi%f(i,j) - phi%f(i-1,j) ) / dx
    end do

    ! Left column for v component
    i = 1
    do j = 2,ny
      v%f(i,j) = vs%f(i,j) - (dt / rho) * ( phi%f(i,j) - phi%f(i,j-1) ) / dy
    enddo

    ! If the domain is periodic in x direction we need to correct
    ! also u on the left and right boundaries
    if (left_boundary == 'periodic') then
      i = 1
      do j = 1,ny
        u%f(i,j) = us%f(i,j) - (dt / rho) * ( phi%f(i,j) - phi%f(i-1,j) ) / dx
      end do

      i = nx + 1
      do j = 1,ny
        u%f(i,j) = us%f(i,j) - (dt / rho) * ( phi%f(i,j) - phi%f(i-1,j) ) / dx
      end do
    end if

    ! If the domain is periodic in y direction we need to correct
    ! also v on the bottom and top boundaries
    if (top_boundary == 'periodic') then
      j = 1
      do i = 1,nx
        v%f(i,j) = vs%f(i,j) - (dt / rho) * ( phi%f(i,j) - phi%f(i,j-1) ) / dy
      end do

      j = ny + 1
      do i = 1,nx
        v%f(i,j) = vs%f(i,j) - (dt / rho) * ( phi%f(i,j) - phi%f(i,j-1) ) / dy
      end do
    endif
    
    call boundary_u(u)
    call boundary_v(v)

  end subroutine project_velocity
  !===============================================================================

  !===============================================================================
  subroutine check(istep)

    ! This subroutine check the divergence of the velocity field
    ! and the instant CFL of the simulation 

    implicit none
    integer, intent(in) :: istep

    integer :: i, j, imax, jmax
    real :: divergence, maximum_divergence, CFL, maximum_CFL

    maximum_divergence = 0.0
    maximum_CFL = 0.0

    do j = 1,ny
      do i = 1,nx
        ! Check divergence
        divergence = (u%f(i+1,j) - u%f(i,j)) / dx + (v%f(i,j+1) - v%f(i,j)) / dy
        if (divergence > maximum_divergence) then
          maximum_divergence = divergence
          imax = i
          jmax = j
        end if

        ! Check CFL
        CFL = 0.5*dt*( (u%f(i+1,j) + u%f(i,j)) / dx + (v%f(i,j+1) + v%f(i,j)) / dy)
        if (CFL > maximum_CFL) maximum_CFL = CFL
      end do
    end do

    if (maximum_divergence > divergence_tol) then
      print *, 'WARNING: divergence is ', maximum_divergence, 'at istep:', istep, &
           'in :', imax, jmax
    endif

    if (maximum_CFL > 0.8) print *, 'WARNING: CFL is:', maximum_CFL 

  end subroutine check
  !===============================================================================

  !===============================================================================
  subroutine destroy_ns_solver()

    ! Free the memory and finalize the simulation
    deallocate(p, phi, u, us, v, vs)
    call destroy_hypre_solver
    call destroy_grid
    call MPI_finalize(ierr)
    stop

  end subroutine destroy_ns_solver

end module navier_stokes
