module boundary_conditions

  use grid_2d
  
contains
 
  subroutine boundary_u(u)

    implicit none
    
    real, intent(inout) :: u(0:,0:)

    integer :: i, j

    ! Left boundary
    if (left_boundary == 'periodic') then
      do j = 0,ny+1
        u(0,j) = u(nx,j)
      end do
    elseif (left_boundary == 'no-slip') then
      do j = 0,ny+1
        u(1,j) = u_l
      end do
    elseif (left_boundary == 'free-slip') then
      do j = 0,ny+1
        u(1,j) = 0.0
      end do
    else
      print *, 'WRONG LEFT BOUNDARY CONDITION ON U'
      call MPI_finalize(ierr)
      stop
    end if

    ! Right boundary
    if (right_boundary == 'periodic') then
      do j = 0,ny+1
        u(nx+2,j) = u(2,j)
      end do
    elseif (right_boundary == 'no-slip') then
      do j = 0,ny+1
        u(nx+1,j) = u_r
      end do
    elseif (right_boundary == 'free-slip') then
      do j = 0,ny+1
        u(nx+1,j) = 0.0
      end do
    else
      print *, 'WRONG RIGHT BOUNDARY CONDITION ON U'
      call MPI_finalize(ierr)
      stop
    end if

    ! Top boundary
    if (top_boundary == 'periodic') then
      do i = 0,nx+2
        u(i,ny+1) = u(i,1)
      end do
    elseif (top_boundary == 'no-slip') then
      do i = 0,nx+2
        u(i,ny+1) = 2.0*u_t - u(i,ny)
      end do
    elseif (top_boundary == 'free-slip') then
      do i = 0,nx+2
        u(i,ny+1) = u(i,ny)
      end do
    else
      print *, 'WRONG TOP BOUNDARY CONDITION ON U'
      call MPI_finalize(ierr)
      stop
    end if

    ! Bottom boundary
    if (bottom_boundary == 'periodic') then
      do i = 0,nx+2
        u(i,0) = u(i,ny)
      end do
    elseif (bottom_boundary == 'no-slip') then
      do i = 0,nx+2
        u(i,0) = 2.0*u_b - u(i,1)
      end do
    elseif (bottom_boundary == 'free-slip') then
      do i = 0,nx+2
        u(i,0) = u(i,1)
      end do
    else
      print *, 'WRONG BOTTOM BOUNDARY CONDITION ON U'
      call MPI_finalize(ierr)
      stop
    end if

  end subroutine boundary_u

  subroutine boundary_v(v)

    implicit none

    real, intent(inout) :: v(0:,0:)

    integer :: i, j

    ! Left boundary
    if (left_boundary == 'periodic') then
      do j = 0,ny+2
        v(0,j) = v(nx,j)
      end do
    elseif (left_boundary == 'no-slip') then
      do j = 0,ny+2
        v(0,j) = 2.0*v_l - v(1,j)
      end do
    elseif (left_boundary == 'free-slip') then
      do j = 0,ny+2
        v(0,j) = v(1,j)
      end do
    else
      print *, 'WRONG LEFT BOUNDARY CONDITION ON V'
      call MPI_finalize(ierr)
      stop
    end if

    ! Right boundary
    if (right_boundary == 'periodic') then
      do j = 0,ny+2
        v(nx+1,j) = v(1,j)
      end do
    elseif (right_boundary == 'no-slip') then
      do j = 0,ny+2
        v(nx+1,j) = 2.0*v_r - v(nx,j)
      end do
    elseif (right_boundary == 'free-slip') then
      do j = 0,ny+2
        v(nx+1,j) = v(nx,j)
      end do
    else
      print *, 'WRONG RIGHT BOUNDARY CONDITION ON V'
      call MPI_finalize(ierr)
      stop
    end if

    ! Top boundary
    if (top_boundary == 'periodic') then
      do i = 0,nx+1
        v(i,ny+2) = v(i,2)
      end do
    elseif (top_boundary == 'no-slip') then
      do i = 0,nx+1
        v(i,ny+1) = v_t
      end do
    elseif (top_boundary == 'free-slip') then
      do i = 0,nx+1
        v(i,ny+1) = 0.0
      end do
    else
      print *, 'WRONG TOP BOUNDARY CONDITION ON V'
      call MPI_finalize(ierr)
      stop
    end if

    ! Bottom boundary
    if (bottom_boundary == 'periodic') then
      do i = 0,nx+1
        v(i,0) = v(i,ny)
      end do
    elseif (bottom_boundary == 'no-slip') then
      do i = 0,nx+1
        v(i,1) = v_b
      end do
    elseif (bottom_boundary == 'free-slip') then
      do i = 0,nx+1
        v(i,0) = 0.0
      end do
    else
      print *, 'WRONG BOTTOM BOUNDARY CONDITION ON V'
      call MPI_finalize(ierr)
      stop
    end if
    
  end subroutine boundary_v

  subroutine boundary_p(p)

    implicit none

    real, intent(inout) :: p(0:,0:)

    integer :: i, j

    ! Left boundary condition
    if (left_boundary == 'periodic') then
      do j = 0,ny+1
        p(0,j) = p(nx,j)
      end do
    elseif (left_boundary == 'no-slip') then
      do j = 0,ny+1
        p(0,j) = p(1,j)
      end do
    elseif (left_boundary == 'free-slip') then
      do j = 0,ny+1
        p(0,j) = p(1,j)
      end do
    else
      print *, 'WRONG LEFT BOUNDARY CONDITION ON P'
      call MPI_finalize(ierr)
      stop
    end if

    ! Right boundary condition
    if (right_boundary == 'periodic') then
      do j = 0,ny+1
        p(nx+1,j) = p(1,j)
      end do
    elseif (right_boundary == 'no-slip') then
      do j = 0,ny+1
        p(nx+1,j) = p(nx,j)
      end do
    elseif (right_boundary == 'free-slip') then
      do j = 0,ny+1
        p(nx+1,j) = p(nx,j)
      end do
    else
      print *, 'WRONG RIGHT BOUNDARY CONDITION ON P'
      call MPI_finalize(ierr)
      stop
    end if
    
    ! Top boundary condition
    if (top_boundary == 'periodic') then
      do i = 0,nx+1
        p(i,ny+1) = p(i,1)
      end do
    elseif (top_boundary == 'no-slip') then
      do i = 0,nx+1
        p(i,ny+1) = p_t
      end do
    elseif (top_boundary == 'free-slip') then
      do i = 0,nx+1
        p(i,ny+1) = p(i,ny)
      end do
    else
      print *, 'WRONG TOP BOUNDARY CONDITION ON P'
      call MPI_finalize(ierr)
      stop
    end if
    
    ! Bottom boundary condition
    if (bottom_boundary == 'periodic') then
      do i = 0,nx+1
        p(i,0) = p(i,ny)
      end do
    elseif (bottom_boundary == 'no-slip') then
      do i = 0,nx+1
        p(i,0) = p(i,1)
      end do
    elseif (bottom_boundary == 'free-slip') then
      do i = 0,nx+1
        p(i,0) = p(i,1)
      end do
    else
      print *, 'WRONG BOTTOM BOUNDARY CONDITION ON P'
      call MPI_finalize(ierr)
      stop
    end if

  end subroutine boundary_p

end module boundary_conditions
