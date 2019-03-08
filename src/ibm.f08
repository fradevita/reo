module ibm
  
  use grid_2d

  ! Number of cylinders
  integer :: nc = 1

  ! Center of cylinders and radius
  real, dimension(:), allocatable :: xc, yc, radius

  ! Variable to tag cells, have dimensions: [nx (+1 in x), ny (+1 in y), 2]
  ! (:,:,1) = tag the grid points based on fluid (2), interface (1) or solid(0) region
  ! (:,:,2) = store the body associated to the grid point
  integer, dimension(:,:,:), allocatable :: ibm_tag_u, ibm_tag_v

  ! Distance function for interpolation
  real, dimension(:,:), allocatable :: phi_u, phi_v, phi_p

  ! Norm to the solid body
  real, dimension(:,:,:), allocatable :: norm_u, norm_v

contains

  subroutine ibm_init()
    ! This subroutine initialize variables for IBM force computation

    implicit none

    ! Allocate all variables
    allocate(xc(nc))
    allocate(yc(nc))
    allocate(radius(nc))
    allocate(ibm_tag_u(nx+1,ny,2))
    allocate(ibm_tag_v(nx,ny+1,2))
    allocate(phi_u(0:nx+2,0:ny+1))
    allocate(phi_v(0:nx+1,0:ny+2))
    allocate(phi_p(nx,ny))
    allocate(norm_u(nx+1,ny,2))
    allocate(norm_v(nx,ny+1,2))

    ! Initialize ibm_index to 0 (solid domain)
    ibm_tag_u = 0
    ibm_tag_v = 0

    ! Initialize distance function to zero
    phi_u = 0.0
    phi_v = 0.0
    phi_p = 0.0

  end subroutine ibm_init
  !===============================================================================

  !===============================================================================
  subroutine ibm_tag()
    ! This subroutine tag the cells based on the location of the cylinders

    implicit none

    ! Local variables
    integer :: i, j, imin, n
    real    :: distance(nc), eps, modnorm

    eps = sqrt(2.0)*dx/2.0
    
    ! Compute the distance function from the solid body in the points of u
    do j = 1,ny
      do i = 1,nx
        do n = 1,nc
          distance(n) = sqrt( (x(i,j) - dx*0.5 - xc(n))**2 + (y(i,j) - yc(n))**2 ) - radius(n)
        end do
        phi_u(i,j) = minval(distance)
        imin = minloc(distance, DIM=1)
        ibm_tag_u(i,j,2) = imin
      end do
    end do

    ! Rigth boundary
    i = nx + 1
    do j = 1,ny
      do n = 1,nc
        distance(n) = sqrt( (x(i-1,j) + dx*0.5 - xc(n))**2 + (y(i-1,j) - yc(n))**2 ) - radius(n)
      end do
      phi_u(i,j) = minval(distance)
      imin = minloc(distance, DIM=1)
      ibm_tag_u(i,j,2) = imin
    end do

    ! Boundary condition on the distance function. Needed for tagging
    call bound_phiu(phi_u)

    ! Tag the cells
    do j = 1,ny
      do i = 1,nx+1
        if ( phi_u(i,j)*phi_u(i-1,j) < 0.0 .or. phi_u(i,j)*phi_u(i+1,j) < 0.0 .or. &
             phi_u(i,j)*phi_u(i,j-1) < 0.0 .or. phi_u(i,j)*phi_u(i,j+1) < 0.0 ) then
          ibm_tag_u(i,j,1) = 1.0 ! interface
        else
          if (phi_u(i,j) > 0.0) then
            ibm_tag_u(i,j,1) = 2.0 ! fluid
          else
            ibm_tag_u(i,j,1) = 0.0 ! solid
          endif
        endif
      end do
    end do

    ! Right boundary
    i = nx + 1
    do j = 1,ny
      if ( phi_u(i,j)*phi_u(i-1,j) < 0.0 .or. phi_u(i,j)*phi_u(i+1,j) < 0.0 .or. &
           phi_u(i,j)*phi_u(i,j-1) < 0.0 .or. phi_u(i,j)*phi_u(i,j+1) < 0.0 ) then
        ibm_tag_u(i,j,1) = 1.0 ! interface
      else
        if (phi_u(i,j) > 0.0) then
          ibm_tag_u(i,j,1) = 2.0 ! fluid
        else
          ibm_tag_u(i,j,1) = 0.0 ! solid
        endif
      endif
    end do
 
    ! Compute the normal vector in the points of u
    do j = 1,ny
      do i = 1,nx
        ! x-component
        norm_u(i,j,1) = ( x(i,j) - dx*0.5 - xc(ibm_tag_u(i,j,2)) ) / &
            sqrt( ( x(i,j) - dx*0.5 - xc(ibm_tag_u(i,j,2)) )**2 +    &
                  ( y(i,j) - yc(ibm_tag_u(i,j,2)) )**2 )

        ! y-component
        norm_u(i,j,2) = ( y(i,j) - yc(ibm_tag_u(i,j,2)) ) /          &
            sqrt( ( x(i,j) - dx*0.5 - xc(ibm_tag_u(i,j,2)) )**2 +    &
                  ( y(i,j) - yc(ibm_tag_u(i,j,2)) )**2 )

        ! Normalization
        modnorm = sqrt( norm_u(i,j,1)**2 + norm_u(i,j,2)**2 ) + 1.0e-16
        norm_u(i,j,1) = norm_u(i,j,1) / modnorm
        norm_u(i,j,2) = norm_u(i,j,2) / modnorm
      end do
    end do

    ! Right boundary
    i = nx + 1
    do j = 1,ny
      ! x-component
      norm_u(i,j,1) = ( x(i-1,j) + dx*0.5 - xc(ibm_tag_u(i,j,2)) ) / &
          sqrt( ( x(i-1,j) + dx*0.5 - xc(ibm_tag_u(i,j,2)) )**2 +    &
                ( y(i-1,j) - yc(ibm_tag_u(i,j,2)) )**2 )
      
      ! y-component
      norm_u(i,j,2) = ( y(i-1,j) - yc(ibm_tag_u(i,j,2)) ) /          &
          sqrt( ( x(i-1,j) + dx*0.5 - xc(ibm_tag_u(i,j,2)) )**2 +    &
                ( y(i-1,j) - yc(ibm_tag_u(i,j,2)) )**2 )
      
      ! Normalization
      modnorm = sqrt( norm_u(i,j,1)**2 + norm_u(i,j,2)**2 ) + 1.0e-16
      norm_u(i,j,1) = norm_u(i,j,1) / modnorm
      norm_u(i,j,2) = norm_u(i,j,2) / modnorm
    end do
    
    ! Compute the distance function from the solid body in the points of v
    do j = 1,ny
      do i = 1,nx
        do n = 1,nc
          distance(n) = sqrt( (x(i,j) - xc(n))**2 + (y(i,j) - dy*0.5 - yc(n))**2 ) - radius(n)
        end do
        phi_v(i,j) = minval(distance)
        imin = minloc(distance, DIM=1)
        ibm_tag_v(i,j,2) = imin
      end do
    end do

    ! Top boundary
    j = ny + 1
    do i = 1,nx
      do n = 1,nc
        distance(n) = sqrt( (x(i,j-1) - xc(n))**2 + (y(i,j-1) + dy*0.5 - yc(n))**2 ) - radius(n)
      end do
      phi_v(i,j) = minval(distance)
      imin = minloc(distance, DIM=1)
      ibm_tag_v(i,j,2) = imin
    end do

    ! Apply boundary condition on distance function
    call bound_phiv(phi_v)
    
    ! Tag the cells
    do j = 1,ny
      do i = 1,nx
        if ( phi_v(i,j)*phi_v(i-1,j) < 0.0 .or. phi_v(i,j)*phi_v(i+1,j) < 0.0 .or. &
             phi_v(i,j)*phi_v(i,j-1) < 0.0 .or. phi_v(i,j)*phi_v(i,j+1) < 0.0 ) then
          ibm_tag_v(i,j,1) = 1.0 ! interface
        else
          if (phi_v(i,j) > 0.0) then
            ibm_tag_v(i,j,1) = 2.0 ! fluid
          else
            ibm_tag_v(i,j,1) = 0.0 ! solid
          endif
        endif
      end do
    end do

    ! Top boundary
    j = ny + 1
    do i = 1,nx
      if ( phi_v(i,j)*phi_v(i-1,j) < 0.0 .or. phi_v(i,j)*phi_v(i+1,j) < 0.0 .or. &
           phi_v(i,j)*phi_v(i,j-1) < 0.0 .or. phi_v(i,j)*phi_v(i,j+1) < 0.0 ) then
        ibm_tag_v(i,j,1) = 1.0 ! interface
      else
        if (phi_v(i,j) > 0.0) then
          ibm_tag_v(i,j,1) = 2.0 ! fluid
        else
          ibm_tag_v(i,j,1) = 0.0 ! solid
        endif
      endif
    end do
    
    ! Compute the normal vector in the points of v
    do j = 1,ny
      do i = 1,nx

        ! x-component
        norm_v(i,j,1) = ( x(i,j) - xc(ibm_tag_v(i,j,2)) ) / &
            sqrt( ( x(i,j) - xc(ibm_tag_v(i,j,2)) )**2 +    &
                  ( y(i,j) - dy*0.5 - yc(ibm_tag_v(i,j,2)) )**2 )

        ! y-component
        norm_v(i,j,2) = ( y(i,j) - dy*0.5 - yc(ibm_tag_v(i,j,2)) ) /          &
            sqrt( ( x(i,j) - xc(ibm_tag_v(i,j,2)) )**2 +    &
                  ( y(i,j) - dy*0.5 - yc(ibm_tag_v(i,j,2)) )**2 )

        ! Normalization
        modnorm = sqrt( norm_v(i,j,1)**2 + norm_v(i,j,2)**2 ) + 1.0e-16
        norm_v(i,j,1) = norm_v(i,j,1) / modnorm
        norm_v(i,j,2) = norm_v(i,j,2) / modnorm
        
      end do
    end do

    ! Top boundary
    j = ny + 1
    do i = 1,nx

      ! x-component
      norm_v(i,j,1) = ( x(i,j-1) - xc(ibm_tag_v(i,j,2)) ) / &
          sqrt( ( x(i,j-1) - xc(ibm_tag_v(i,j,2)) )**2 +    &
                ( y(i,j-1) + dy*0.5 - yc(ibm_tag_v(i,j,2)) )**2 )

      ! y-component
      norm_v(i,j,2) = ( y(i,j-1) + dy*0.5 - yc(ibm_tag_v(i,j,2)) ) /          &
          sqrt( ( x(i,j-1) - xc(ibm_tag_v(i,j,2)) )**2 +    &
                ( y(i,j-1) + dy*0.5 - yc(ibm_tag_v(i,j,2)) )**2 )

      ! Normalization
      modnorm = sqrt( norm_v(i,j,1)**2 + norm_v(i,j,2)**2 ) + 1.0e-16
      norm_v(i,j,1) = norm_v(i,j,1) / modnorm
      norm_v(i,j,2) = norm_v(i,j,2) / modnorm
      
    end do
    
    ! Compute the distance function from the solid body in the points of p
    do j = 1,ny
      do i = 1,nx
        do n = 1,nc
          distance(n) = sqrt( (x(i,j) - xc(n))**2 + (y(i,j) - yc(n))**2 ) - radius(n)
        enddo
        phi_p(i,j) = minval(distance)
      end do
    end do
    
  end subroutine ibm_tag
  !===============================================================================

  !===============================================================================
  function ibm_force(i, j, dir, rhs, v, dt) result(f)
    ! Compute the ibm_force in direction dir in the point i, j

    implicit none

    ! Input variables
    integer, intent(in) :: i, j, dir
    real   , intent(in) :: rhs, v, dt

    ! Output variable
    real :: f

    select case(dir)
    case(1)
      if (ibm_tag_u(i,j,1) == 2) then ! Inside the fluid the force is zero
        f = 0.0
      elseif (ibm_tag_u(i,j,1) == 1) then ! On the interface, compute interpolated velocity
        f = -rhs + (interpolate(i,j,dir,ibm_tag_u(i,j,2)) - v) / dt
      else ! inside the solid
        f = -rhs - v / dt
      endif
    case(2)
      if (ibm_tag_v(i,j,1) == 2) then ! Inside the fluid the force is zero
        f = 0.0
      elseif (ibm_tag_v(i,j,1) == 1) then ! On the interface, compute interpolated velocity
        f = -rhs + (interpolate(i,j,dir,ibm_tag_v(i,j,2)) - v) / dt
      else ! inside the solid
        f = -rhs - v / dt
      endif
    end select

  end function ibm_force
  !===============================================================================

  !===============================================================================
  function interpolate(i, j, dir, n) result(v_int)
    ! Compute the value of the velocity that we need to impose with the ibm force

    use navier_stokes_pub, only : u, v
    
    implicit none

    ! Input variables
    integer, intent(in) :: i, j, dir, n

    ! Output variable
    real :: v_int

    ! Local variables
    real    :: b, c, xl, yl, vx, vy, discr, xx, yy

    v_int = 0.0
    vx = 0.0
    vy = 0.0

    ! ************************* WEIGHTED INTERPOLATION IN Y AND Z *****************************
    select case(dir)
    case(1) ! u
      ! Local coordinates on the point on u
      if( i < nx+1 ) then
        xx = x(i,j) - dx * 0.5
        yy = y(i,j)
      else
        xx = x(i-1,j) + dx * 0.5
        yy = y(i-1,j)
      endif
      
      ! First compute interpolation in i
      b = -2.0*yc(n)
      c = xc(n)**2 + yy**2 + yc(n)**2 - 2.0*yy*yc(n) - radius(n)**2
      discr = b**2 - 4.0*c
 
      if (discr <= 0.0) then ! No intersection
        vx = 0.0
      else
        if (phi_u(i,j) > 0.0) then ! On liquid side of the interface
          if (norm_u(i,j,1) > 0.0) then
            xl = (-b + sqrt(b**2 - 4.0*c))*0.5
            vx = u%f(i+1,j)*(xx-xl)/(xx+dx-xl)
          else
            xl = (-b - sqrt(b**2 - 4.0*c))*0.5
            vx = u%f(i-1,j)*(xx-xl)/(xx-dx-xl)
          endif
        else
          vx = 0.0
        endif
      endif

      ! Then compute interpolation in j
      b = -2.0*yc(n)
      c = yc(n)**2 + xx**2 + xc(n)**2 - 2.0*xx*xc(n) - radius(n)**2
      discr = b**2 - 4.0*c

      if (discr <= 0.0) then ! No intersection
        vy = 0.0
      else
        if (phi_u(i,j) > 0.0) then ! On liquid side
          if (norm_u(i,j,2) > 0.0) then
            yl = (-b + sqrt(b**2 - 4*c))*0.5
            vy = u%f(i,j+1)*(yy-yl)/(yy+dy-yl)
          else
            yl = (-b - sqrt(b**2 - 4*c))*0.5
            vy = u%f(i,j-1)*(yy-yl)/(yy-dy-yl)
          endif
        else
          vy = 0.0
        endif
      endif

      ! Weighted value of the two interpolation
      v_int = vx*norm_u(i,j,1)**2 + vy*norm_u(i,j,2)**2
 
    case(2) ! v
      ! Local coordinates on the point of v
      if (j < ny + 1) then
        xx = x(i,j)
        yy = y(i,j) - dy * 0.5
      else
        xx = x(i,j-1)
        yy = y(i,j-1) + dy * 0.5
      endif
      
      ! First interpolate in j
      b = -2.0*xc(n)
      c = xc(n)**2 + yy**2 + yc(n)**2 - 2.0*yy*yc(n) - radius(n)**2
      discr = b**2 - 4.0*c

      if (discr <= 0.0) then ! No intersection
        vx = 0.0
      else
        if (phi_v(i,j) >0.0) then ! On liquid side
          if (norm_v(i,j,1) > 0.0) then
            xl = (-b + sqrt(b**2 - 4*c))*0.5
            vx = v%f(i+1,j)*(xx-xl)/(xx+dx-xl)
          else
            xl = (-b - sqrt(b**2 - 4*c))*0.5
            vx = v%f(i-1,j)*(xx-xl)/(xx-dx-xl)
          endif
        else
          vx = 0.0
        endif
      endif

      ! Then interpolate in j
      b = -2.0*yc(n)
      c = yc(n)**2 + xx**2 + xc(n)**2 - 2.0*xx*xc(n) - radius(n)**2
      discr = b**2 - 4.0*c

      if (discr <= 0.0) then ! No intersection
        vy = 0.0
      else
        if(phi_v(i,j) > 0.0) then ! On liquid side
          if (norm_v(i,j,2) > 0.0) then
            yl = (-b + sqrt(b**2 - 4*c))*0.5
            vy = v%f(i,j+1)*(yy-yl)/(yy+dy-yl)
          else
            yl = (-b - sqrt(b**2 - 4*c))*0.5
            vy = v%f(i,j-1)*(yy-yl)/(yy-dy-yl)
          endif
        else
          vy = 0.0
        endif
      endif

      ! Weighted value of the two interpolation
      v_int = vx*norm_v(i,j,1)**2 + vy*norm_v(i,j,2)**2

    end select
    
  end function interpolate
  !===============================================================================
  
  !===============================================================================
  subroutine bound_phiu(phi)

    implicit none
    
    real, intent(inout) :: phi(0:,0:)

    integer :: i, j

    ! Left boundary
    if (left_boundary == 'periodic') then
      do j = 0,ny+1
        phi(0,j) = phi(nx,j)
      end do
    else
      do j = 0,ny+1
        phi(0,j) = phi(1,j)
      end do
    end if

    ! Right boundary
    if (right_boundary == 'periodic') then
      do j = 0,ny+1
        phi(nx+2,j) = phi(2,j)
      end do
    else
      do j = 0,ny+1
        phi(nx+2,j) = phi(nx+1,j)
      end do
    end if

    ! Top boundary
    if (top_boundary == 'periodic') then
      do i = 0,nx+2
        phi(i,ny+1) = phi(i,1)
      end do
    else
      do i = 0,nx+2
        phi(i,ny+1) = phi(i,ny)
      end do
    end if

    ! Bottom boundary
    if (bottom_boundary == 'periodic') then
      do i = 0,nx+2
        phi(i,0) = phi(i,ny)
      end do
    else
      do i = 0,nx+2
        phi(i,0) = phi(i,1)
      end do
    end if

  end subroutine bound_phiu
  !===============================================================================

  !===============================================================================
  subroutine bound_phiv(phi)

    implicit none

    real, intent(inout) :: phi(0:,0:)

    integer :: i, j
    
   ! Left boundary
    if (left_boundary == 'periodic') then
      do j = 0,ny+2
        phi(0,j) = phi(nx,j)
      end do
    else
      do j = 0,ny+2
        phi(0,j) = phi(1,j)
      end do
    end if

    ! Right boundary
    if (right_boundary == 'periodic') then
      do j = 0,ny+2
        phi(nx+1,j) = phi(1,j)
      end do
    else
      do j = 0,ny+2
        phi(nx+1,j) = phi(nx,j)
      end do
    end if

    ! Top boundary
    if (top_boundary == 'periodic') then
      do i = 0,nx+1
        phi(i,ny+2) = phi(i,2)
      end do
    else
      do i = 0,nx+1
        phi(i,ny+2) = phi(i,ny+2)
      end do
    end if

    ! Bottom boundary
    if (bottom_boundary == 'periodic') then
      do i = 0,nx+1
        phi(i,0) = phi(i,ny)
      end do
    else
      do i = 0,nx+1
        phi(i,0) = phi(i,1)
      end do
    end if
    
  end subroutine bound_phiv
  !===============================================================================
  
end module ibm
