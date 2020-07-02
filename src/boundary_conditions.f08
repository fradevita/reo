module boundary_conditions

contains

  !===============================================================================
  subroutine boundary(f)

    use grid_2D
!    use navier_stokes_pub

    implicit none

    include 'mpif.h'

    ! Input / Output variable
    type(field), intent(inout) :: f(nbox)

    ! Local variable
    integer :: status(MPI_STATUS_SIZE)
    integer :: i, j, n

    box_cycle : do n = 1,nbox

      select case(f(n)%location)
      case(0) ! Cell center

        ! Left boundary
        if (f(n)%left == 'neumann') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(n)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%left == 'periodic') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(n)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%left == 'dirichlet') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = 2.0*f(n)%l(j) - f(n)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%left == 'internal') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(boxarray(n)%lb)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%left == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 5', f(n)%left
          stop
        end if

        ! Right boundary
        if (f(n)%right == 'neumann') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(n)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%right == 'periodic') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(n)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%right == 'dirichlet') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = 2.0*f(n)%r(j) - f(n)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%right == 'internal') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(boxarray(n)%rb)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%right == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 6'
          stop
        end if

        ! Top boundary
        if (f(n)%top == 'neumann') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(n)%t(i)*boxarray(n)%delta + f(n)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%top == 'periodic') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(n)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%top == 'dirichlet') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = 2.0*f(n)%t(i) - f(n)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%top == 'internal') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(boxarray(n)%tb)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%top == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 7'
          stop
        end if
  
        ! Bottom boundary
        if (f(n)%bottom == 'neumann') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(n)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%bottom == 'periodic') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(n)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%bottom == 'dirichlet') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = 2.0*f(n)%b(i) - f(n)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%bottom == 'internal') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(boxarray(n)%bb)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%bottom == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 8'
          stop
        end if
  
      case(1) ! x face
  
        ! Left boundary
        if (f(n)%left == 'neumann') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(n)%f(f(n)%lo(1)+1,j)
          end do
        elseif (f(n)%left == 'periodic') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(n)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%left == 'dirichlet') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(n)%l(j)
          end do
        elseif (f(n)%left == 'internal') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(boxarray(n)%lb)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%left == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 9'
          stop
        end if
  
        ! Right boundary
        if (f(n)%right == 'neumann') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(n)%f(f(n)%up(1)-1,j)
          end do
        elseif (f(n)%right == 'periodic') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(n)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%right == 'dirichlet') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(n)%r(j)
          end do
        elseif (f(n)%right == 'internal') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(boxarray(n)%rb)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%right == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 10'
          stop
        end if
  
        ! Top boundary
        if (f(n)%top == 'neumann') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(n)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%top == 'periodic') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(n)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%top == 'dirichlet') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = 2.0*f(n)%t(i) - f(n)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%top == 'internal') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(boxarray(n)%tb)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%top == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 11'
          stop
        end if
  
        ! Bottom boundary
        if (f(n)%bottom == 'neumann') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(n)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%bottom == 'periodic') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(n)%f(i,f(n)%up(2))
          end do
        elseif( f(n)%bottom == 'dirichlet') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = 2.0*f(n)%b(i) - f(n)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%bottom == 'internal') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(boxarray(n)%bb)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%bottom == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 12'
          stop
        end if
  
      case(2) ! y face
  
        ! Left boundary
        if (f(n)%left == 'neumann') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(n)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%left == 'periodic') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = f(n)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%left == 'dirichlet') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%lo(1)-1,j) = 2.0*f(n)%l(j) - f(n)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%left == 'internal') then
          do j = f(n)%lo(2)-1,f(n)%up(2)!+1
            f(n)%f(f(n)%lo(1)-1,j) = f(boxarray(n)%lb)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%left == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 13'
        end if
  
        ! Right boundary
        if (f(n)%right == 'neumann') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(n)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%right == 'periodic') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(n)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%right == 'dirichlet') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = 2.0*f(n)%r(j) - f(n)%f(f(n)%up(1),j)
          end do
        elseif (f(n)%right == 'internal') then
          do j = f(n)%lo(2)-1,f(n)%up(2)+1
            f(n)%f(f(n)%up(1)+1,j) = f(boxarray(n)%rb)%f(f(n)%lo(1),j)
          end do
        elseif (f(n)%right == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 14'
          stop
        end if
  
        ! Top boundary
        if (f(n)%top == 'neumann') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(n)%f(i,f(n)%up(2)-1)
          end do
        elseif (f(n)%top == 'periodic') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(n)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%top == 'dirichlet') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = 2.0*f(n)%t(i)
          end do
        elseif (f(n)%top == 'internal') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%up(2)+1) = f(boxarray(n)%tb)%f(i,f(n)%lo(2))
          end do
        elseif (f(n)%top == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 15'
          stop
        end if
  
        ! Bottom boundary
        if (f(n)%bottom == 'neumann') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(n)%f(i,f(n)%lo(2)+1)
          end do
        elseif (f(n)%bottom == 'periodic') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(n)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%bottom == 'dirichlet') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = 2.0*f(n)%b(i)
          end do
        elseif (f(n)%bottom == 'internal') then
          do i = f(n)%lo(1)-1,f(n)%up(1)+1
            f(n)%f(i,f(n)%lo(2)-1) = f(boxarray(n)%bb)%f(i,f(n)%up(2))
          end do
        elseif (f(n)%bottom == 'halo') then
          ! do nothing for now
        else
          print *, 'ERROR 16'
          stop
        end if

      end select

      ! Halo exchage for parallel simulations
      call mpi_barrier(mpi_common_world, ierr)
      if (num_procs > 1) then

        if (f(n)%left == 'halo' .or. f(n)%right == 'halo') then
          do j = f(n)%lo(2),f(n)%up(2)
            ! Mando in avanti
            call mpi_sendrecv(f(n)%f(f(n)%up(1),j), 1, mpi_real8, boxarray(n)%rp, 0, &
                              f(n)%f(f(n)%lo(1)-1,j), 1, mpi_real8, boxarray(n)%lp, 0, &
                              mpi_comm_world, status, ierr)
            ! Mando in dietro
            call mpi_sendrecv(f(n)%f(f(n)%lo(1),j), 1, mpi_real8, boxarray(n)%lp, 0, &
                              f(n)%f(f(n)%up(1)+1,j), 1, mpi_real8, boxarray(n)%rp, 0, &
                              mpi_comm_world, status, ierr)
          end do
        endif

        if (f(n)%top == 'halo' .or. f(n)%bottom == 'halo') then
          do i = f(n)%lo(1),f(n)%up(1)
            ! Mando in su
            call mpi_sendrecv(f(n)%f(i,f(n)%up(2)), 1, mpi_real8, boxarray(n)%tp, 0, &
                              f(n)%f(i,f(n)%lo(2)-1), 1, mpi_real8, boxarray(n)%bp, 0, &
                              mpi_comm_world, status, ierr)
            ! Mando in giù
            call mpi_sendrecv(f(n)%f(i,f(n)%lo(2)), 1, mpi_real8, boxarray(n)%bp, 0, &
                              f(n)%f(i,f(n)%up(2)+1), 1, mpi_real8, boxarray(n)%tp, 0, &
                              mpi_comm_world, status, ierr)
          end do
        endif

      endif  

    end do box_cycle

  end subroutine boundary
  !===============================================================================

end module boundary_conditions
