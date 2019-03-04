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

  real, dimension(:,:), allocatable :: phi

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

    ! We need to allocate all the fields
    allocate(p(0:nx+1,0:ny+1))
    allocate(phi(0:nx+1,0:ny+1))
    allocate(u(0:nx+2,0:ny+1))
    allocate(us(0:nx+2,0:ny+1))
    allocate(v(0:nx+1,0:ny+2))
    allocate(vs(0:nx+1,0:ny+2))

    ! We set all the field to zero
    p = 0.0
    phi = 0.0
    u = 0.0
    us = 0.0
    v = 0.0
    vs = 0.0

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
          phi(i,j) = rhs(i + (j-1)*ny)
        end do
      end do
      call boundary_p(phi)

      ! Then we project the predicted velocity field to the divergence-free
      ! velocity
      call project_velocity()

      ! We update the pressure
      !p = p + phi
      p = p + 2.0*phi
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
        uuip  = 0.25 * ( u(ip,j) + u(i,j) ) * ( u(ip,j) + u(i,j)  )
        uuim  = 0.25 * ( u(im,j) + u(i,j) ) * ( u(im,j) + u(i,j)  )
        uvjp  = 0.25 * ( u(i,jp) + u(i,j) ) * ( v(i,jp) + v(im,jp))
        uvjm  = 0.25 * ( u(i,jm) + u(i,j) ) * ( v(i,j) + v(im,j)  )
        du(i,j) = (-uuip + uuim) / dx + (-uvjp + uvjm) / dy + &
             mu*((u(ip,j) - 2.0*u(i,j) + u(im,j))/dxq + &
             (u(i,jp) - 2.0*u(i,j) + u(i,jm))/dyq)

        ! y direction
        uvip  = 0.25 * ( u(ip,jm) + u(ip,j) ) * ( v(i,j) + v(ip,j) )
        uvim  = 0.25 * ( u(i,jm) + u(i,j) ) * ( v(i,j) + v(im,j) )
        vvjp  = 0.25 * ( v(i,j) + v(i,jp) ) * ( v(i,j) + v(i,jp) )
        vvjm  = 0.25 * ( v(i,j) + v(i,jm) ) * ( v(i,j) + v(i,jm) )
        dv(i,j) = (-uvip + uvim) / dx + (-vvjp + vvjm) / dy + &
             mu*((v(ip,j) - 2.0*v(i,j) + v(im,j))/dxq + &
             (v(i,jp) - 2.0*v(i,j) + v(i,jm))/dyq)      
      end do
    end do

    ! Bottom row of the grid for u
    j = 1
    jp = j + 1
    jm = j - 1
    do i = 2,nx
      ip = i + 1
      im = i - 1
      uuip  = 0.25 * ( u(ip,j) + u(i,j) ) * ( u(ip,j) + u(i,j)  )
      uuim  = 0.25 * ( u(im,j) + u(i,j) ) * ( u(im,j) + u(i,j)  )
      uvjp  = 0.25 * ( u(i,jp) + u(i,j) ) * ( v(i,jp) + v(im,jp))
      uvjm  = 0.25 * ( u(i,jm) + u(i,j) ) * ( v(i,j) + v(im,j)  )
      du(i,j) = (-uuip + uuim) / dx + (-uvjp + uvjm) / dy + &
           mu*((u(ip,j) - 2.0*u(i,j) + u(im,j))/dxq + &
           (u(i,jp) - 2.0*u(i,j) + u(i,jm))/dyq)
    end do

    ! Left column of the grid for v
    i = 1
    ip = i + 1
    im = i - 1
    do j = 2,ny
      jp = j + 1
      jm = j - 1
      uvip  = 0.25 * ( u(ip,jm) + u(ip,j) ) * ( v(i,j) + v(ip,j) )
      uvim  = 0.25 * ( u(i,jm) + u(i,j) ) * ( v(i,j) + v(im,j) )
      vvjp  = 0.25 * ( v(i,j) + v(i,jp) ) * ( v(i,j) + v(i,jp) )
      vvjm  = 0.25 * ( v(i,j) + v(i,jm) ) * ( v(i,j) + v(i,jm) )
      dv(i,j) = (-uvip + uvim) / dx + (-vvjp + vvjm) / dy + &
           mu*((v(ip,j) - 2.0*v(i,j) + v(im,j))/dxq + &
           (v(i,jp) - 2.0*v(i,j) + v(i,jm))/dyq)
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
        uuip  = 0.25 * ( u(ip,j) + u(i,j) ) * ( u(ip,j) + u(i,j)  )
        uuim  = 0.25 * ( u(im,j) + u(i,j) ) * ( u(im,j) + u(i,j)  )
        uvjp  = 0.25 * ( u(i,jp) + u(i,j) ) * ( v(i,jp) + v(im,jp))
        uvjm  = 0.25 * ( u(i,jm) + u(i,j) ) * ( v(i,j) + v(im,j)  )
        du(i,j) = (-uuip + uuim) / dx + (-uvjp + uvjm) / dy + &
             mu*((u(ip,j) - 2.0*u(i,j) + u(im,j))/dxq + &
             (u(i,jp) - 2.0*u(i,j) + u(i,jm))/dyq)
      end do

      i = nx + 1
      ip = i + 1
      im = i - 1
      do j = 1,ny
        jp = j + 1
        jm = j - 1
        uuip  = 0.25 * ( u(ip,j) + u(i,j) ) * ( u(ip,j) + u(i,j)  )
        uuim  = 0.25 * ( u(im,j) + u(i,j) ) * ( u(im,j) + u(i,j)  )
        uvjp  = 0.25 * ( u(i,jp) + u(i,j) ) * ( v(i,jp) + v(im,jp))
        uvjm  = 0.25 * ( u(i,jm) + u(i,j) ) * ( v(i,j) + v(im,j)  )
        du(i,j) = (-uuip + uuim) / dx + (-uvjp + uvjm) / dy + &
             mu*((u(ip,j) - 2.0*u(i,j) + u(im,j))/dxq + &
             (u(i,jp) - 2.0*u(i,j) + u(i,jm))/dyq)
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
        uvip  = 0.25 * ( u(ip,jm) + u(ip,j) ) * ( v(i,j) + v(ip,j) )
        uvim  = 0.25 * ( u(i,jm) + u(i,j) ) * ( v(i,j) + v(im,j) )
        vvjp  = 0.25 * ( v(i,j) + v(i,jp) ) * ( v(i,j) + v(i,jp) )
        vvjm  = 0.25 * ( v(i,j) + v(i,jm) ) * ( v(i,j) + v(i,jm) )
        dv(i,j) = (-uvip + uvim) / dx + (-vvjp + vvjm) / dy + &
             mu*((v(ip,j) - 2.0*v(i,j) + v(im,j))/dxq + &
             (v(i,jp) - 2.0*v(i,j) + v(i,jm))/dyq)      
      end do

      j = ny + 1
      jp = j + 1
      jm = j - 1
      do i = 1,nx
        ip = i + 1
        im = i - 1
        uvip  = 0.25 * ( u(ip,jm) + u(ip,j) ) * ( v(i,j) + v(ip,j) )
        uvim  = 0.25 * ( u(i,jm) + u(i,j) ) * ( v(i,j) + v(im,j) )
        vvjp  = 0.25 * ( v(i,j) + v(i,jp) ) * ( v(i,j) + v(i,jp) )
        vvjm  = 0.25 * ( v(i,j) + v(i,jm) ) * ( v(i,j) + v(i,jm) )
        dv(i,j) = (-uvip + uvim) / dx + (-vvjp + vvjm) / dy + &
             mu*((v(ip,j) - 2.0*v(i,j) + v(im,j))/dxq + &
             (v(i,jp) - 2.0*v(i,j) + v(i,jm))/dyq)      
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
        rhs = 1.5 * du(i,j) - 0.5 * du_o(i,j) - (p(i,j) - p(i-1,j)) / (rho*dx) + Sx
        us(i,j) = u(i,j) + dt * rhs

#ifdef IBM
        us(i,j) = us(i,j) + dt * ibm_force(i,j,1,rhs,u(i,j),dt)
#endif

        ! y direction
        rhs = 1.5 * dv(i,j) - 0.5 * dv_o(i,j) - (p(i,j) - p(i,j-1)) / (rho*dy) + Sy
        vs(i,j) = v(i,j) + dt * rhs

#ifdef IBM
        vs(i,j) = vs(i,j) + dt * ibm_force(i,j,2,rhs,v(i,j),dt)  
#endif

      end do
    end do

    ! bottom row of the grid for us
    j = 1
    do i = 2,nx
      rhs = 1.5 * du(i,j) - 0.5 * du_o(i,j) - (p(i,j) - p(i-1,j)) / (rho*dx) + Sx
      us(i,j) = u(i,j) + dt * rhs

#ifdef IBM
        us(i,j) = us(i,j) + dt * ibm_force(i,j,1,rhs,u(i,j),dt)
#endif

    end do

    ! left column of the grid for vs
    i = 1
    do j = 2,ny
      rhs = 1.5 * dv(i,j) - 0.5 * dv_o(i,j) - (p(i,j) - p(i,j-1)) / (rho*dy) + Sy
      vs(i,j) = v(i,j) + dt * rhs

#ifdef IBM
        vs(i,j) = vs(i,j) + dt * ibm_force(i,j,2,rhs,v(i,j),dt)  
#endif

    end do

    ! If the domain is periodic in x direction we need to solve
    ! also us on the left and right boundaries
    if (left_boundary == 'periodic') then
      i = 1
      do j = 1,ny
        rhs = 1.5 * du(i,j) - 0.5 * du_o(i,j) - (p(i,j) - p(i-1,j)) / (rho*dx) + Sx
        us(i,j) = u(i,j) + dt * rhs

#ifdef IBM
        us(i,j) = us(i,j) + dt * ibm_force(i,j,1,rhs,u(i,j),dt)
#endif

      end do

      i = nx + 1
      do j = 1,ny
        rhs = 1.5 * du(i,j) - 0.5 * du_o(i,j) - (p(i,j) - p(i-1,j)) / (rho*dx) + Sx
        us(i,j) = u(i,j) + dt * rhs

#ifdef IBM
        us(i,j) = us(i,j) + dt * ibm_force(i,j,1,rhs,u(i,j),dt)
#endif

      end do
    end if

    ! If the domain is periodic in y direction we need to solve
    ! also v on the bottom and top boundaries
    if (top_boundary == 'periodic') then
      j = 1
      do i = 1,nx
        rhs = 1.5 * dv(i,j) - 0.5 * dv_o(i,j) - (p(i,j) - p(i,j-1)) / (rho*dy) + Sy
        vs(i,j) = v(i,j) + dt * rhs

#ifdef IBM
        vs(i,j) = vs(i,j) + dt * ibm_force(i,j,2,rhs,v(i,j),dt)
#endif

      end do

      j = ny + 1
      do i = 1,nx
        rhs = 1.5 * dv(i,j) - 0.5 * dv_o(i,j) - (p(i,j) - p(i,j-1)) / (rho*dy) + Sy
        vs(i,j) = v(i,j) + dt * rhs

#ifdef IBM
        vs(i,j) = vs(i,j) + dt * ibm_force(i,j,2,rhs,v(i,j),dt) 
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
        rhs(i+(j-1)*ny) = (rho/dt) * ((us(i+1,j) - us(i,j)) / dx + &
             (vs(i,j+1) - vs(i,j)) / dy )
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
        u(i,j) = us(i,j) - (dt / rho) * ( phi(i,j) - phi(i-1,j) ) / dx

        ! y direction
        v(i,j) = vs(i,j) - (dt / rho) * ( phi(i,j) - phi(i,j-1) ) / dy
      end do
    end do

    ! Bottom row for u component
    j = 1
    do i = 2,nx
      u(i,j) = us(i,j) - (dt / rho) * ( phi(i,j) - phi(i-1,j) ) / dx
    end do

    ! Left column for v component
    i = 1
    do j = 2,ny
      v(i,j) = vs(i,j) - (dt / rho) * ( phi(i,j) - phi(i,j-1) ) / dy
    enddo

    ! If the domain is periodic in x direction we need to correct
    ! also u on the left and right boundaries
    if (left_boundary == 'periodic') then
      i = 1
      do j = 1,ny
        u(i,j) = us(i,j) - (dt / rho) * ( phi(i,j) - phi(i-1,j) ) / dx
      end do

      i = nx + 1
      do j = 1,ny
        u(i,j) = us(i,j) - (dt / rho) * ( phi(i,j) - phi(i-1,j) ) / dx
      end do
    end if

    ! If the domain is periodic in y direction we need to correct
    ! also v on the bottom and top boundaries
    if (top_boundary == 'periodic') then
      j = 1
      do i = 1,nx
        v(i,j) = vs(i,j) - (dt / rho) * ( phi(i,j) - phi(i,j-1) ) / dy
      end do

      j = ny + 1
      do i = 1,nx
        v(i,j) = vs(i,j) - (dt / rho) * ( phi(i,j) - phi(i,j-1) ) / dy
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
        divergence = (u(i+1,j) - u(i,j)) / dx + (v(i,j+1) - v(i,j)) / dy
        if (divergence > maximum_divergence) then
          maximum_divergence = divergence
          imax = i
          jmax = j
        end if

        ! Check CFL
        CFL = 0.5*dt*( (u(i+1,j) + u(i,j)) / dx + (v(i,j+1) + v(i,j)) / dy)
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
