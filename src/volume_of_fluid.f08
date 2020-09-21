module volume_of_fluid

  use grid_2d, only : field, myarray, dp

  implicit none

  type(field), dimension(:), allocatable :: vof, nx, ny, h, d
  type(myarray), dimension(:), allocatable :: lxx, lyy, curv

  real(kind=dp) :: epsilon = 1.0e-16, cx, cy, a10, a01, a20, a02
  real(kind=dp) :: beta = 1.0
  real(kind=dp) :: cut = 1.0e-8
  logical :: quadratic = .true.

  private
  public :: vof, init_vof, reconstruct, h, compute_norm, nx, ny, curv, &
            advect_vof, set_vof, destroy_vof

contains

  !===============================================================================
  subroutine init_vof()

    use grid_2D, only : boxarray, nbox

    implicit none
 
    integer :: n

    allocate(vof(nbox))
    allocate(nx(nbox))
    allocate(ny(nbox))
    allocate(lxx(nbox))
    allocate(lyy(nbox))
    allocate(curv(nbox))
    allocate(h(nbox))
    allocate(d(nbox))

    do n = 1,nbox
      vof(n)%lo = [1, 1]
      vof(n)%up = [boxarray(n)%nx, boxarray(n)%ny]
      nx(n)%lo = [1, 1]
      nx(n)%up = [boxarray(n)%nx, boxarray(n)%ny]
      ny(n)%lo = [1, 1]
      ny(n)%up = [boxarray(n)%nx, boxarray(n)%ny]
      h(n)%lo = [1, 1]
      h(n)%up = [boxarray(n)%nx, boxarray(n)%ny]
      d(n)%lo = [1, 1]
      d(n)%up = [boxarray(n)%nx, boxarray(n)%ny]
      if (boxarray(n)%left_boundary == 'no-slip' .or. &
          boxarray(n)%left_boundary  == 'inflow') then
        vof(n)%left = 'dirichlet'
        nx(n)%left = 'dirichlet'
        ny(n)%left = 'dirichlet'
        h(n)%left = 'dirichlet'
        d(n)%left = 'dirichlet'
      elseif (boxarray(n)%left_boundary == 'outflow' .or. &
              boxarray(n)%left_boundary == 'free-slip') then
        vof(n)%left = 'neumann'
        nx(n)%left = 'neumann'
        ny(n)%left = 'neumann'
        h(n)%left = 'neumann'
        d(n)%left = 'neumann'
      elseif (boxarray(n)%left_boundary == 'internal') then
        vof(n)%left = 'internal'
        nx(n)%left = 'internal'
        ny(n)%left = 'internal'
        h(n)%left = 'internal'
        d(n)%left = 'internal'
      elseif (boxarray(n)%left_boundary == 'halo') then
        vof(n)%left = 'halo'
        nx(n)%left = 'halo'
        ny(n)%left = 'halo'
        h(n)%left = 'halo'
        d(n)%left = 'halo'
      elseif (boxarray(n)%left_boundary == 'periodic') then
        vof(n)%left = 'periodic'
        nx(n)%left = 'periodic'
        ny(n)%left = 'periodic'
        h(n)%left = 'periodic'
        d(n)%left = 'periodic'
      endif
      if (boxarray(n)%right_boundary == 'no-slip' .or. &
          boxarray(n)%right_boundary  == 'inflow') then
        vof(n)%right = 'dirichlet'
        nx(n)%right = 'dirichlet'
        ny(n)%right = 'dirichlet'
        h(n)%right = 'dirichlet'
        d(n)%right = 'dirichlet'
      elseif (boxarray(n)%right_boundary == 'outflow' .or. &
              boxarray(n)%right_boundary == 'free-slip') then
        vof(n)%right = 'neumann'
        nx(n)%right = 'neumann'
        ny(n)%right = 'neumann'
        h(n)%right = 'neumann'
        d(n)%right = 'neumann'
      elseif (boxarray(n)%right_boundary == 'internal') then
        vof(n)%right = 'internal'
        nx(n)%right = 'internal'
        ny(n)%right = 'internal'
        h(n)%right = 'internal'
        d(n)%right = 'internal'
      elseif (boxarray(n)%right_boundary == 'halo') then
        vof(n)%right = 'halo'
        nx(n)%right = 'halo'
        ny(n)%right = 'halo'
        h(n)%right = 'halo'
        d(n)%right = 'halo'
      elseif (boxarray(n)%right_boundary == 'periodic') then
        vof(n)%right = 'periodic'
        nx(n)%right = 'periodic'
        ny(n)%right = 'periodic'
        h(n)%right = 'periodic'
        d(n)%right = 'periodic'
      endif
      if (boxarray(n)%top_boundary == 'no-slip' .or. &
          boxarray(n)%top_boundary  == 'inflow') then
        vof(n)%top = 'dirichlet'
        nx(n)%top = 'dirichlet'
        ny(n)%top = 'dirichlet'
        h(n)%top = 'dirichlet'
        d(n)%top = 'dirichlet'
      elseif (boxarray(n)%top_boundary == 'outflow' .or. &
              boxarray(n)%top_boundary == 'free-slip') then
        vof(n)%top = 'neumann'
        nx(n)%top = 'neumann'
        ny(n)%top = 'neumann'
        h(n)%top = 'neumann'
        d(n)%top = 'neumann'
      elseif (boxarray(n)%top_boundary == 'internal') then
        vof(n)%top = 'internal'
        nx(n)%top = 'internal'
        ny(n)%top = 'internal'
        h(n)%top = 'internal'
        d(n)%top = 'internal'
      elseif (boxarray(n)%top_boundary == 'halo') then
        vof(n)%top = 'halo'
        nx(n)%top = 'halo'
        ny(n)%top = 'halo'
        h(n)%top = 'halo'
        d(n)%top = 'halo'
      elseif (boxarray(n)%top_boundary == 'periodic') then
        vof(n)%top = 'periodic'
        nx(n)%top = 'periodic'
        ny(n)%top = 'periodic'
        h(n)%top = 'periodic'
        d(n)%top = 'periodic'
      endif
      if (boxarray(n)%bottom_boundary == 'no-slip' .or. &
          boxarray(n)%bottom_boundary  == 'inflow') then
        vof(n)%bottom = 'dirichlet'
        nx(n)%bottom = 'dirichlet'
        ny(n)%bottom = 'dirichlet'
        h(n)%bottom = 'dirichlet'
        d(n)%bottom = 'dirichlet'
      elseif (boxarray(n)%bottom_boundary == 'outflow' .or. &
              boxarray(n)%bottom_boundary == 'free-slip') then
        vof(n)%bottom = 'neumann'
        nx(n)%bottom = 'neumann'
        ny(n)%bottom = 'neumann'
        h(n)%bottom = 'neumann'
        d(n)%bottom = 'neumann'
      elseif (boxarray(n)%bottom_boundary == 'internal') then
        vof(n)%bottom = 'internal'
        nx(n)%bottom = 'internal'
        ny(n)%bottom = 'internal'
        h(n)%bottom = 'internal'
        d(n)%bottom = 'internal'
      elseif (boxarray(n)%bottom_boundary == 'halo') then
        vof(n)%bottom = 'halo'
        nx(n)%bottom = 'halo'
        ny(n)%bottom = 'halo'
        h(n)%bottom = 'halo'
        d(n)%bottom = 'halo'
      elseif (boxarray(n)%bottom_boundary == 'periodic') then
        vof(n)%bottom = 'periodic'
        nx(n)%bottom = 'periodic'
        ny(n)%bottom = 'periodic'
        h(n)%bottom = 'periodic'
        d(n)%bottom = 'periodic'
      endif
      allocate(vof(n)%l(vof(n)%lo(2)-1:vof(n)%up(2)+1))
      allocate(vof(n)%r(vof(n)%lo(2)-1:vof(n)%up(2)+1))
      allocate(vof(n)%t(vof(n)%lo(1)-1:vof(n)%up(1)+1))
      allocate(vof(n)%b(vof(n)%lo(1)-1:vof(n)%up(1)+1))
      vof(n)%l = 0.
      vof(n)%r = 0.
      vof(n)%t = 0.
      vof(n)%b = 0.
      vof(n)%location = 0
      allocate(nx(n)%l(nx(n)%lo(2)-1:nx(n)%up(2)+1))
      allocate(nx(n)%r(nx(n)%lo(2)-1:nx(n)%up(2)+1))
      allocate(nx(n)%t(nx(n)%lo(1)-1:nx(n)%up(1)+1))
      allocate(nx(n)%b(nx(n)%lo(1)-1:nx(n)%up(1)+1))
      nx(n)%l = 0.
      nx(n)%r = 0.
      nx(n)%t = 0.
      nx(n)%b = 0.
      nx(n)%location = 0
      allocate(ny(n)%l(ny(n)%lo(2)-1:ny(n)%up(2)+1))
      allocate(ny(n)%r(ny(n)%lo(2)-1:ny(n)%up(2)+1))
      allocate(ny(n)%t(ny(n)%lo(1)-1:ny(n)%up(1)+1))
      allocate(ny(n)%b(ny(n)%lo(1)-1:ny(n)%up(1)+1))
      ny(n)%l = 0.
      ny(n)%r = 0.
      ny(n)%t = 0.
      ny(n)%b = 0.
      ny(n)%location = 0
      allocate(h(n)%l(h(n)%lo(2)-1:h(n)%up(2)+1))
      allocate(h(n)%r(h(n)%lo(2)-1:h(n)%up(2)+1))
      allocate(h(n)%t(h(n)%lo(1)-1:h(n)%up(1)+1))
      allocate(h(n)%b(h(n)%lo(1)-1:h(n)%up(1)+1))
      h(n)%l = 0.
      h(n)%r = 0.
      h(n)%t = 0.
      h(n)%b = 0.
      h(n)%location = 0
      allocate(d(n)%l(d(n)%lo(2)-1:d(n)%up(2)+1))
      allocate(d(n)%r(d(n)%lo(2)-1:d(n)%up(2)+1))
      allocate(d(n)%t(d(n)%lo(1)-1:d(n)%up(1)+1))
      allocate(d(n)%b(d(n)%lo(1)-1:d(n)%up(1)+1))
      d(n)%l = 0.
      d(n)%r = 0.
      d(n)%t = 0.
      d(n)%b = 0.
      d(n)%location = 0
      allocate(vof(n)%f(0:vof(n)%up(1)+1,0:vof(n)%up(2)+1))
      vof(n)%f = 0.
      allocate(nx(n)%f(0:vof(n)%up(1)+1,0:vof(n)%up(2)+1))
      nx(n)%f = 0.
      allocate(ny(n)%f(0:vof(n)%up(1)+1,0:vof(n)%up(2)+1))
      ny(n)%f = 0.
      allocate(lxx(n)%f(0:vof(n)%up(1)+1,0:vof(n)%up(2)+1))
      lxx(n)%f = 0.
      allocate(lyy(n)%f(0:vof(n)%up(1)+1,0:vof(n)%up(2)+1))
      lyy(n)%f = 0.
      allocate(curv(n)%f(0:vof(n)%up(1)+1,0:vof(n)%up(2)+1))
      curv(n)%f = 0.
      allocate(h(n)%f(0:vof(n)%up(1)+1,0:vof(n)%up(2)+1))
      h(n)%f = 0.
      allocate(d(n)%f(0:vof(n)%up(1)+1,0:vof(n)%up(2)+1))
      d(n)%f = 0.
    end do

  end subroutine
  !===============================================================================

  !===============================================================================
  subroutine compute_norm()

    use grid_2d, only : nbox, boxarray, dp
    use boundary_conditions

    implicit none

    integer :: n, i, j, c
    real(kind=dp) :: mx(4), mxc, my(4), myc, normx(4), normy(4), delta

    do n = 1,nbox

      nx(n)%f = 0.
      ny(n)%f = 0.
  
      delta = boxarray(n)%delta
 
      do j = 1,boxarray(n)%ny
        do i = 1,boxarray(n)%nx

          ! X component
          ! i-1/2,j-1/2
          mx(1) = 0.5*( vof(n)%f(i,j-1) + vof(n)%f(i,j) - &
                  vof(n)%f(i-1,j-1) - vof(n)%f(i-1,j) ) / delta 
          ! i-1/2,j+1/2
          mx(2) = 0.5*( vof(n)%f(i,j) + vof(n)%f(i,j+1) - &
                  vof(n)%f(i-1,j) - vof(n)%f(i-1,j+1) ) / delta 
          ! i+1/2,j+1/2
          mx(3) = 0.5*( vof(n)%f(i+1,j) + vof(n)%f(i+1,j+1) - &
                  vof(n)%f(i,j) - vof(n)%f(i,j+1) ) / delta 
          ! i+1/2,j-1/2
          mx(4) = 0.5*( vof(n)%f(i+1,j-1) + vof(n)%f(i+1,j) - &
                  vof(n)%f(i,j-1) - vof(n)%f(i,j) ) / delta 
  
          mxc = 0.25*(mx(1) + mx(2) + mx(3) + mx(4))

          ! Y component
          ! i-1/2,j-1/2
          my(1) = 0.5*( vof(n)%f(i-1,j) + vof(n)%f(i,j) - &
                  vof(n)%f(i-1,j-1) - vof(n)%f(i,j-1) ) / delta 
          ! i-1/2,j+1/2
          my(2) = 0.5*( vof(n)%f(i-1,j+1) + vof(n)%f(i,j+1) - &
                  vof(n)%f(i-1,j) - vof(n)%f(i,j) ) / delta 
          ! i+1/2,j+1/2
          my(3) = 0.5*( vof(n)%f(i,j+1) + vof(n)%f(i+1,j+1) - &
                  vof(n)%f(i,j) - vof(n)%f(i+1,j) ) / delta 
          ! i+1/2,j-1/2
          my(4) = 0.5*( vof(n)%f(i,j) + vof(n)%f(i+1,j) - &
                  vof(n)%f(i,j-1) - vof(n)%f(i+1,j-1) ) / delta 

          myc = 0.25*(my(1) + my(2) + my(3) + my(4))
  
          ! Norm vector components 
          do c = 1,4
            normx(c) = mx(c)/sqrt(mx(c)**2 + my(c)**2 + epsilon)   
            normy(c) = my(c)/sqrt(mx(c)**2 + my(c)**2 + epsilon)   
          end do
          nx(n)%f(i,j) = mxc / sqrt( mxc**2 + myc**2 + epsilon)
          ny(n)%f(i,j) = myc / sqrt( mxc**2 + myc**2 + epsilon)

          ! Curvature components
          if (quadratic) then
            lxx(n)%f(i,j) = 0.5*delta*(normx(4) + normx(3) - normx(2) - normx(1)) 
            lyy(n)%f(i,j) = 0.5*delta*(normy(2) + normy(3) - normy(1) - normy(4))
            curv(n)%f(i,j) = -(lxx(n)%f(i,j) + lyy(n)%f(i,j))/delta**2
          else
            lxx(n)%f(i,j) = 0. 
            lyy(n)%f(i,j) = 0.
            curv(n)%f(i,j) = -(lxx(n)%f(i,j) + lyy(n)%f(i,j))/delta**2
          endif
        end do
      end do
    end do

    call boundary(nx)
    call boundary(ny)

  end subroutine 
  !===============================================================================

  !===============================================================================
  subroutine reconstruct()
    ! Compute the normal using the Young method

    use grid_2D, only : nbox, boxarray, dp
    use boundary_conditions

    implicit none

    integer :: n, i, j
    real(kind=dp) :: delta, rp, rm, A, Bp, Bm, Q
    real(kind=dp) :: aa, bb, cc

    ! Two points Gaussian quadrature in the interval [0,1]
    rp = 0.5*(1. + 1./sqrt(3.))
    rm = 0.5*(1. - 1./sqrt(3.))

    do n = 1,nbox
  
      delta = boxarray(n)%delta
 
      do j = 1,boxarray(n)%ny
        do i = 1,boxarray(n)%nx

          ! Reconstruct only for values larger than the cut value
          if (vof(n)%f(i,j) < cut .or. vof(n)%f(i,j) > (1. - cut) ) then

            h(n)%f(i,j) = vof(n)%f(i,j)
            d(n)%f(i,j) = 0.

          else

            if ( abs(nx(n)%f(i,j)) .eq. max(abs(nx(n)%f(i,j)),abs(ny(n)%f(i,j))) ) then
              cx = 0.
              cy = 1.
            else
              cx = 1.
              cy = 0.
            endif

            a10 = nx(n)%f(i,j) - 0.5*cx*lxx(n)%f(i,j)
            a01 = ny(n)%f(i,j) - 0.5*cy*lyy(n)%f(i,j)
            a20 = 0.5*cx*lxx(n)%f(i,j)
            a02 = 0.5*cy*lyy(n)%f(i,j)

            A = (1. - cx)*exp(2.*beta*a10) + &
                (1. - cy)*exp(2.*beta*a01)

            Bp = (1. - cx)*exp(2.*beta*P(0.d0,rp)) + &
                 (1. - cy)*exp(2.*beta*P(rp,0.d0))

            Bm = (1. - cx)*exp(2.*beta*P(0.d0,rm)) + &
                 (1. - cy)*exp(2.*beta*P(rm,0.d0))

            Q = (1. - cx)*exp(2.*beta*a10*(2.d0*vof(n)%f(i,j) - 1.d0)) + &
                (1. - cy)*exp(2.*beta*a01*(2.d0*vof(n)%f(i,j) - 1.d0))

            aa = A*Bm*Bp*(A - Q)
            bb = A*(Bp + Bm)*(1.d0 - Q)
            cc = 1.d0 - A*Q

            d(n)%f(i,j) = log(solve_quadratic(aa,bb,cc))/(2.d0*beta)
 
            h(n)%f(i,j) = 0.5d0*(1.0d0 + tanh(beta*(P(0.5d0,0.5d0) + d(n)%f(i,j))))
          endif

        end do
      end do
    end do
    call boundary(h)
    call boundary(d)

  end subroutine
  !===============================================================================

  !===============================================================================
  function P(x, y) result(r)

    implicit none

    real(kind=dp), intent(in) :: x, y
    real(kind=dp) :: r

    r = cx*a20*x**2 + cy*a02*y**2 + a10*x + a01*y

  end function
  !===============================================================================

  !===============================================================================
  function solve_quadratic(a, b, c) result(x)

    real(kind=dp), intent(in) :: a, b, c
    real(kind=dp) :: x1, x2, x

    x1 = (-b + sqrt(b**2 - 4.d0*a*c))/(2.d0*a)
    x2 = (-b - sqrt(b**2 - 4.d0*a*c))/(2.d0*a)

    x = max(x1,x2)

  end function
  !===============================================================================

  !===============================================================================
  subroutine advect_vof(u, v, dt)

    use grid_2D, only : nbox, boxarray, myarray, dp
    use boundary_conditions

    implicit none

    real(kind=dp), intent(in) :: dt
    type(field), intent(in) :: u(:), v(:)

    integer :: i, j, n
    real(kind=dp) :: delta, fp, fm
    type(myarray), allocatable :: vof1(:), vof2(:)

    allocate(vof1(nbox))
    allocate(vof2(nbox))

    ! x direction
    do n = 1,nbox
      delta = boxarray(n)%delta

      allocate(vof1(n)%f(1:vof(n)%up(1),1:vof(n)%up(2)))
      do j = 1,vof(n)%up(2)
        do i = 1,vof(n)%up(1)
          fp = compute_flux(n, 1, i  , j, u(n)%f(i+1,j), dt, delta)
          fm = compute_flux(n, 1, i-1, j, u(n)%f(i,j)  , dt, delta)
          vof1(n)%f(i,j) = (vof(n)%f(i,j) - (fp - fm)/delta)/(1. - dt* &
                            (u(n)%f(i+1,j) - u(n)%f(i,j))/delta)
        end do
      end do

      vof(n)%f(1:vof(n)%up(1),1:vof(n)%up(2)) = vof1(n)%f(:,:)
    end do

    ! Get h from vof1
    call boundary(vof)
    call compute_norm()
    call reconstruct()

    ! y direction
    do n = 1,nbox
      allocate(vof2(n)%f(1:vof(n)%up(1),1:vof(n)%up(2)))
      do j = 1,vof(n)%up(2)
        do i = 1,vof(n)%up(1)
          fp = compute_flux(n, 2, i, j  , v(n)%f(i,j+1), dt, delta)
          fm = compute_flux(n, 2, i, j-1, v(n)%f(i,j)  , dt, delta)
          vof2(n)%f(i,j) = (vof1(n)%f(i,j) - (fp - fm)/delta)/(1. - dt* &
                            (v(n)%f(i,j+1) - v(n)%f(i,j))/delta)
        end do
      end do

      ! Compute vof at time n+1
      do j = 1,vof(n)%up(2)
        do i = 1,vof(n)%up(1)
          vof(n)%f(i,j) = vof2(n)%f(i,j) - dt*(vof1(n)%f(i,j)* &
                           (u(n)%f(i+1,j) - u(n)%f(i,j))/delta + &
                           vof2(n)%f(i,j)*(v(n)%f(i,j+1) - v(n)%f(i,j))/delta)
          if (vof(n)%f(i,j) > 1.) then
            vof(n)%f(i,j) = 1.
          elseif (vof(n)%f(i,j) < 0) then 
            vof(n)%f(i,j) = 0.
          endif
        end do
      end do
      deallocate(vof1(n)%f,vof2(n)%f)
    end do
    deallocate(vof1,vof2)

    ! Get h from the new vof
    call boundary(vof)
    call compute_norm()
    call reconstruct()

  end subroutine advect_vof
  !===============================================================================

  !===============================================================================
  function compute_flux(n, dir, i, j, u, dt, delta) result(f)

    implicit none

    integer, intent(in) :: n, dir, i, j
    real(kind=dp), intent(in) :: u, dt, delta

    integer :: sgn, ii, jj
    real(kind=dp) :: f, xa, xb, ya, yb

    ! Select integration intervals
    if (dir == 1) then
      if (u >= 0.)then
        xa = 1. - dt*u/delta
        xb = 1.
        sgn = 1
        ii = i
        jj = j
      else
        xa = 0.
        xb = -dt*u/delta
        sgn = -1
        ii = i + 1
        jj = j
      endif
      ya = 0.
      yb = 1.
    elseif (dir == 2) then
     xa = 0.
     xb = 1.
     if (u >= 0.) then
       ya = 1. - dt*u/delta
       yb = 1.
       sgn = 1
       ii = i
       jj = j
     else
       ya = 0.
       yb = -dt*u/delta
       sgn = -1
       ii = i
       jj = j + 1
     endif
    endif

    if (vof(n)%f(ii,jj) < cut .or. vof(n)%f(ii,jj) > (1. - cut)) then
      f = sgn*delta*vof(n)%f(ii,jj)*(xb - xa)*(yb - ya)
    else
      if (abs(nx(n)%f(ii,jj)) == max(abs(nx(n)%f(ii,jj)),abs(ny(n)%f(ii,jj)))) then
        cx = 0.
        cy = 1.
        a10 = nx(n)%f(ii,jj) - 0.5*cx*lxx(n)%f(ii,jj)
        a01 = ny(n)%f(ii,jj) - 0.5*cy*lyy(n)%f(ii,jj)
        a20 = 0.5*cx*lxx(n)%f(ii,jj)
        a02 = 0.5*cy*lyy(n)%f(ii,jj)
        f = sgn*delta*0.5*(yb - ya)*(0.5*(xb - xa + 1./(beta*a10)*log( &
            cosh(beta*(P(xb,(yb - ya)/(2.*sqrt(3.)) + 0.5*(ya + yb)) + d(n)%f(ii,jj)))/ &
            cosh(beta*(P(xa,(yb - ya)/(2.*sqrt(3.)) + 0.5*(ya + yb)) + d(n)%f(ii,jj))))) + &
                                   0.5*(xb - xa + 1./(beta*a10)*log( &
            cosh(beta*(P(xb,-(yb - ya)/(2.*sqrt(3.)) + 0.5*(ya + yb)) + d(n)%f(ii,jj)))/ &
            cosh(beta*(P(xa,-(yb - ya)/(2.*sqrt(3.)) + 0.5*(ya + yb)) + d(n)%f(ii,jj))))))
      else
        cx = 1.
        cy = 0.
        a10 = nx(n)%f(ii,jj) - 0.5*cx*lxx(n)%f(ii,jj)
        a01 = ny(n)%f(ii,jj) - 0.5*cy*lyy(n)%f(ii,jj)
        a20 = 0.5*cx*lxx(n)%f(ii,jj)
        a02 = 0.5*cy*lyy(n)%f(ii,jj)
        f = sgn*delta*0.5*(xb - xa)*(0.5*(yb - ya + 1./(beta*a01)*log( &
            cosh(beta*(P((xb - xa)/(2.*sqrt(3.)) + 0.5*(xa + xb),yb) + d(n)%f(ii,jj)))/ &
            cosh(beta*(P((xb - xa)/(2.*sqrt(3.)) + 0.5*(xa + xb),ya) + d(n)%f(ii,jj))))) + &
                                 0.5*(yb - ya + 1./(beta*a01)*log( &
            cosh(beta*(P(-(xb - xa)/(2.*sqrt(3.)) + 0.5*(xa + xb),yb) + d(n)%f(ii,jj)))/ &
            cosh(beta*(P(-(xb - xa)/(2.*sqrt(3.)) + 0.5*(xa + xb),ya) + d(n)%f(ii,jj))))))
      endif
    endif

  end function compute_flux
  !===============================================================================

  !===============================================================================
  subroutine set_vof(geometry, xc, yc, r)

    use grid_2d, only : nbox, boxarray

    implicit none

    real(kind=dp), intent(in) :: xc, yc, r
    character(len=6), intent(in) :: geometry

    integer :: n, i, j, counter, ii, jj
    real(kind=dp) :: x, y, delta, distance, x0, y0, xx, yy

    if (geometry == 'circle') then

      do n = 1,nbox
        delta = boxarray(n)%delta

        do j = 1,boxarray(n)%ny
          y = boxarray(n)%p0(2) + (j - 0.5)*delta
          do i = 1,boxarray(n)%nx
            x = boxarray(n)%p0(1) + (i - 0.5)*delta
            distance = sqrt((x - xc)**2 + (y - yc)**2) - r
            if (distance > 0.5*sqrt(2.)*delta) then
              vof(n)%f(i,j) = 0.
            elseif (distance < -0.5*sqrt(2.)*delta) then
              vof(n)%f(i,j) = 1.
            else
              x0 = x - 0.5*delta
              y0 = y - 0.5*delta
              counter = 0
              do jj = 1,100
                yy = y0 + (jj - 0.5)*delta/100.
                do ii = 1,100
                  xx = x0 + (ii - 0.5)*delta/100.
                  distance = (xx - xc)**2 + (yy - yc)**2 - r**2
                  if (distance < -0.5*sqrt(2.)*delta/100.) counter = counter + 1
                end do
              end do
              vof(n)%f(i,j) = counter / float(100*100)
            endif
          end do
        end do
      end do
      call compute_norm
      call reconstruct
    end if

  end subroutine
  !===============================================================================

  !===============================================================================
  subroutine destroy_vof()

    use grid_2D, only : nbox

    implicit none

    integer :: n

    do n = 1,nbox
      deallocate(vof(n)%f,vof(n)%l,vof(n)%r,vof(n)%t,vof(n)%b)
      deallocate(nx(n)%f,nx(n)%l,nx(n)%r,nx(n)%t,nx(n)%b)
      deallocate(ny(n)%f,ny(n)%l,ny(n)%r,ny(n)%t,ny(n)%b)
      deallocate(h(n)%f,h(n)%l,h(n)%r,h(n)%t,h(n)%b)
      deallocate(d(n)%f,d(n)%l,d(n)%r,d(n)%t,d(n)%b)
      deallocate(lxx(n)%f,lyy(n)%f,curv(n)%f)
    end do
    deallocate(vof,nx,ny,h,d,lxx,lyy,curv)

  end subroutine destroy_vof
  !===============================================================================

end module volume_of_fluid
