module boundary_conditions

contains

  subroutine boundary(f)

    use navier_stokes_pub

    implicit none

    type(field), intent(inout) :: f

    integer :: i, j
    
    select case(f%location)
    case(0) ! Cell center

      ! Left boundary
      if (f%left == 'neumann') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = f%f(f%lo(1),j)
        end do
      elseif (f%left == 'periodic') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = f%f(f%up(1),j)
        end do
      elseif (f%left == 'dirichlet') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = 2.0*f%l(j) - f%f(f%lo(1),j)
        end do
      else
        print *, 'ERROR 5', f%left
        stop
      end if

      ! Right boundary
      if (f%right == 'neumann') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = f%f(f%up(1),j)
        end do
      elseif (f%right == 'periodic') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = f%f(f%lo(1),j)
        end do
      elseif (f%right == 'dirichlet') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = 2.0*f%r(j) - f%f(f%up(1),j)
        end do
      else
        print *, 'ERROR 6'
        stop
      end if

      ! Top boundary
      if (f%top == 'neumann') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = f%f(i,f%up(2))
        end do
      elseif (f%top == 'periodic') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = f%f(i,f%lo(2))
        end do
      elseif (f%top == 'dirichlet') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = 2.0*f%t(i) - f%f(i,f%up(2))
        end do
      else
        print *, 'ERROR 7'
        stop
      end if

      ! Bottom boundary
      if (f%bottom == 'neumann') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = f%f(i,f%lo(2))
        end do
      elseif (f%bottom == 'periodic') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = f%f(i,f%up(2))
        end do
      elseif (f%bottom == 'dirichlet') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = 2.0*f%b(i) - f%f(i,f%lo(2))
        end do
      else
        print *, 'ERROR 8'
        stop
      end if

    case(1) ! x face

      ! Left boundary
      if (f%left == 'neumann') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = f%f(f%lo(1)+1,j)
        end do
      elseif (f%left == 'periodic') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = f%f(f%up(1),j)
        end do
      elseif (f%left == 'dirichlet') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = f%l(j)
        end do
      else
        print *, 'ERROR 9'
        stop
      end if

      ! Right boundary
      if (f%right == 'neumann') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = f%f(f%up(1)-1,j)
        end do
      elseif (f%right == 'periodic') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = f%f(f%lo(1),j)
        end do
      elseif (f%right == 'dirichlet') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = f%r(j)
        end do
      else
        print *, 'ERROR 10'
        stop
      end if

      ! Top boundary
      if (f%top == 'neumann') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = f%f(i,f%up(2))
        end do
      elseif (f%top == 'periodic') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = f%f(i,f%lo(2))
        end do
      elseif (f%top == 'dirichlet') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = 2.0*f%t(i) - f%f(i,f%up(2))
        end do
      else
        print *, 'ERROR 11'
        stop
      end if

      ! Bottom boundary
      if (f%bottom == 'neumann') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = f%f(i,f%lo(2))
        end do
      elseif (f%bottom == 'periodic') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = f%f(i,f%up(2))
        end do
      elseif( f%bottom == 'dirichlet') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = 2.0*f%b(i) - f%f(i,f%lo(2))
        end do
      else
        print *, 'ERROR 12'
        stop
      end if

    case(2) ! y face

      ! Left boundary
      if (f%left == 'neumann') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = f%f(f%lo(1),j)
        end do
      elseif (f%left == 'periodic') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = f%f(f%up(1),j)
        end do
      elseif (f%left == 'dirichlet') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%lo(1)-1,j) = 2.0*f%l(j) - f%f(f%lo(1),j)
        end do
      else
        print *, 'ERROR 13'
      end if

      ! Right boundary
      if (f%right == 'neumann') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = f%f(f%up(1),j)
        end do
      elseif (f%right == 'periodic') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = f%f(f%lo(1),j)
        end do
      elseif (f%right == 'dirichlet') then
        do j = f%lo(2)-1,f%up(2)+1
          f%f(f%up(1)+1,j) = 2.0*f%r(j) - f%f(f%up(1),j)
        end do
      else
        print *, 'ERROR 14'
        stop
      end if

      ! Top boundary
      if (f%top == 'neumann') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = f%f(i,f%up(2)-1)
        end do
      elseif (f%top == 'periodic') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = f%f(i,f%lo(2))
        end do
      elseif (f%top == 'dirichlet') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%up(2)+1) = 2.0*f%t(i)
        end do
      else
        print *, 'ERROR 15'
        stop
      end if

      ! Bottom boundary
      if (f%bottom == 'neumann') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = f%f(i,f%lo(2)+1)
        end do
      elseif (f%bottom == 'periodic') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = f%f(i,f%up(2))
        end do
      elseif (f%bottom == 'dirichlet') then
        do i = f%lo(1)-1,f%up(1)+1
          f%f(i,f%lo(2)-1) = 2.0*f%b(i)
        end do
      else
        print *, 'ERROR 16'
        stop
      end if

    end select

  end subroutine boundary

end module boundary_conditions
