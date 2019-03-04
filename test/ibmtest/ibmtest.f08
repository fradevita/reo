program lid

  use grid_2d
  include "navier_stokes.h"
  use ibm

  implicit none
  integer :: i, j

  ! We need an auxiliary field to check steady state
  real, dimension(:,:), allocatable :: uold
  
  ! Boundary conditions
  u_t = 1.0
  
  ! Set the number of points and the domain size
  nx = 64
  ny = 64
  Lx = 1.0
  Ly = 1.0
  x0 = -0.5
  y0 = -0.5

  allocate(uold(nx,ny))

  ! Create the grid
  call create_grid()

  ! First we need to initialize the solver
  call init_ns_solver()

#ifdef IBM
  ! Location and size of the cylinder
  nc = 1
  xc(1) = 0.0
  yc(1) = 0.0
  radius(1) = 0.2
#endif
  
  ! We set the viscosity to 1
  mu = 1.0e-3

  ! We compute the solution up to time T = 300
  dt = min(0.3 * dx / 1.0, 0.1*(dx**2 + dy**2) / mu)
  !dt = 1.0e-4
  Tmax = 300
  nstep = int(Tmax / dt)

  ! We check for steady state
  event_i => e_istep
  event_output => output

  ! We run the simulation
  call solve()

contains
  
  subroutine e_istep()

    do j = 1,ny
      do i = 2,nx
        uold(i,j) = u(i,j)
      end do
    end do
    
  end subroutine e_istep
    
  subroutine output()

    implicit none
    
    ! Check for steady-state
    real :: diff, uc, vc, diffmax, vort
    logical :: steady
    integer :: i, j, out_id, xprof, yprof
    
    steady = .true.
    diffmax = 0.0
 
    do j = 1,ny
      do i = 2,nx
        diff = abs(u(i,j) - uold(i,j))
        if (diff > 1.0e-8) steady = .false.
        if (diff > diffmax) diffmax = diff
      end do
    end do

    write(log,*) 'maximum difference: ', diffmax
   
    if (steady) then
      ! We output the horizontal porfile of y-component of the velocity and
      ! vertical profile of x-component of the velocity on the centerline of the box
      open(newunit = xprof, file = 'xprof')
      i = nx / 2
      do j = 1,ny
        write(xprof,*) y(i,j), u(i,j)
      end do
      close(xprof)

      open(newunit = yprof, file = 'yprof')
      j = ny / 2
      do i = 1,nx
        write(yprof,*) x(i,j), v(i,j)
      end do
      close(yprof)

      ! We output also the norm of the velocity field at the steady state
      open(newunit = out_id, file = 'out')
      do j = 1,ny
        do i = 1,nx
          uc = 0.5*(u(i,j) + u(i+1,j))
          vc = 0.5*(v(i,j) + v(i,j+1))
          vort = 0.25*(u(i,j+1) - u(i,j-1) + u(i+1,j+1) - u(i+1,j-1)) / dx + &
                 0.25*(v(i+1,j) - v(i-1,j) + v(i+1,j+1) - v(i+1,j+1)) / dy
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
