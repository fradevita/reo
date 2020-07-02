module multiphase

  implicit none

  real :: rho0 = 1., rho1 = 1., mu0 = 1., mu1 = 1., sigma = 0.
  logical :: solvemultiphase = .false.

  private
  public :: update_material_properties, update_coefficient_matrix, compute_surface_tension_force
  public :: rho0, rho1, mu0, mu1, sigma, solvemultiphase

contains

  !===============================================================================
  subroutine update_material_properties(rho, mu)

    use grid_2D, only : nbox, field, boxarray
    use volume_of_fluid, only : vof
    use boundary_conditions

    implicit none

    type(field), intent(inout) :: rho(nbox), mu(nbox)

    integer :: n, i, j

    do n = 1,nbox
      do j = 1,boxarray(n)%ny
        do i = 1,boxarray(n)%nx
          rho(n)%f(i,j) = rho0*(1. - vof(n)%f(i,j)) + rho1*vof(n)%f(i,j)
          mu(n)%f(i,j) = mu0*(1. - vof(n)%f(i,j)) + mu1*vof(n)%f(i,j)
        end do
      end do
    end do

    call boundary(rho)
    call boundary(mu)

  end subroutine update_material_properties
  !===============================================================================

  !===============================================================================
  subroutine compute_surface_tension_force(fisgmax, lox, upx, fisgmay, loy, upy, n)

    use grid_2d, only : boxarray, dp
    use volume_of_fluid, only : vof, curv

    implicit none

    integer, intent(in) :: lox(2), upx(2), loy(2), upy(2), n
    real(kind=dp), intent(out) :: fisgmax(lox(1):upx(1),lox(2):upx(2))
    real(kind=dp), intent(out) :: fisgmay(loy(1):upy(1),loy(2):upy(2))

    integer :: i, j, im, jm

    do j = lox(2),upx(2)
      jm = j - 1
      do i = lox(1),upx(1)
        im = i - 1
        fisgmax(i,j) = 0.5*sigma*(curv(n)%f(i,j) + curv(n)%f(im,j))* &
                       (vof(n)%f(i,j) - vof(n)%f(im,j))/boxarray(n)%delta
      end do
    end do

    do j = loy(2),upy(2)
      jm = j - 1
      do i = loy(1),upy(1)
        im = i - 1
        fisgmay(i,j) = 0.5*sigma*(curv(n)%f(i,j) + curv(n)%f(i,jm))* &
                       (vof(n)%f(i,j) - vof(n)%f(i,jm))/boxarray(n)%delta
      end do
    end do

  end subroutine compute_surface_tension_force
  !===============================================================================

  !===============================================================================
  subroutine update_coefficient_matrix(rho)

    use grid_2D, only : field, box, nbox, boxarray, ierr, dp
    use hypre_solver
    use mpi
 
    implicit none

    ! Input / output variables
    type(field) :: rho(:)

    ! Local variables
    integer :: i, j, nentries, nvalues, n, np, ii, jj
    integer, dimension(2) :: ilower, iupper
    integer, dimension(:), allocatable :: stencil_indices
    real(kind=dp) :: rhoim, rhoip, rhojm, rhojp
    real(kind=dp), dimension(:), allocatable :: values
    type(box) :: lb
 
    ! Set the matrix cefficients for each box.
    box_cycle : do n = 1,nbox

      ! Select the box
      lb = boxarray(n)
      ilower = lb%ilower
      iupper = lb%iupper

      ! nx*ny grid points each with 5 stencil entries
      nentries = 5
      nvalues = lb%nx*lb%ny*nentries
      allocate(values(nvalues))
      allocate(stencil_indices(nentries))
      stencil_indices = [0, 1, 2, 3, 4]
      jj = 1
      do i = 1,nvalues,nentries
        ii = i/nentries + 1 - (jj-1)*lb%nx
        rhoip = 2./(rho(n)%f(ii+1,jj) + rho(n)%f(ii,jj))
        rhoim = 2./(rho(n)%f(ii-1,jj) + rho(n)%f(ii,jj))
        rhojp = 2./(rho(n)%f(ii,jj+1) + rho(n)%f(ii,jj))
        rhojm = 2./(rho(n)%f(ii,jj-1) + rho(n)%f(ii,jj))
        values(i)  = -1.0*(rhoip + rhoim + rhojp + rhojm)
        values(i+1) = 1.0*rhoim
        values(i+2) = 1.0*rhoip
        values(i+3) = 1.0*rhojm
        values(i+4) = 1.0*rhojp
        if (ii == lb%nx) jj = jj + 1
      end do
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, stencil_indices, values, ierr)
      deallocate(values, stencil_indices)
    
      ! Set the coefficients outside of the domain and on the boundaries
      if (periodic(1) == 0) then
        ! Values on the left of the box
        ilower = [lb%ilower(1), lb%ilower(2)]
        iupper = [lb%ilower(1), lb%iupper(2)]
        np = lb%iupper(2) - lb%ilower(2) + 1
        allocate(values(np))
        allocate(stencil_indices(1))
        ! If the left border is not internal or halo set to zero the stencil 1
        if (lb%left /= 'internal' .and. lb%left /= 'halo') then          
          values = 0.0
          stencil_indices(1) = 1
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        endif

        ! Then modify the central point of the stencil
        stencil_indices(1) = 0
        if (lb%left == 'dirichlet') then
          do j = 1,np
            rhoip = 2./(rho(n)%f(2,j) + rho(n)%f(1,j))
            rhoim = 2./(rho(n)%f(0,j) + rho(n)%f(1,j))
            rhojp = 2./(rho(n)%f(1,j+1) + rho(n)%f(1,j))
            rhojm = 2./(rho(n)%f(1,j-1) + rho(n)%f(1,j))
            values(j) = -1.0*(rhoip + 2.*rhoim + rhojp + rhojm)
          end do
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        elseif (lb%left == 'neumann') then
          do j = 1,np
            rhoip = 2./(rho(n)%f(2,j) + rho(n)%f(1,j))
            rhoim = 2./(rho(n)%f(0,j) + rho(n)%f(1,j))
            rhojp = 2./(rho(n)%f(1,j+1) + rho(n)%f(1,j))
            rhojm = 2./(rho(n)%f(1,j-1) + rho(n)%f(1,j))
            values(j) = -1.0*(rhoip + rhojp + rhojm)
          end do
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        elseif (lb%left == 'periodic' .or. lb%left == 'internal' .or. lb%left == 'halo') then
          ! do nothing
        else
          print *, 'WRONG LEFT BC HYPRE'
        end if

        ! Values on the right of the box
        ilower = [lb%iupper(1), lb%ilower(2)]
        iupper = [lb%iupper(1), lb%iupper(2)]
        
        ! If the right border is not internal or halo set to zero the stencil 2
        if (lb%right /= 'internal' .and. lb%right /= 'halo') then
          values = 0.0
          stencil_indices(1) = 2
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        endif

        ! Then modify the central point of the stencil
        stencil_indices(1) = 0
        if (lb%right == 'dirichlet') then
          do j = 1,np
            rhoip = 2./(rho(n)%f(lb%nx+1,j) + rho(n)%f(lb%nx,j))
            rhoim = 2./(rho(n)%f(lb%nx-1,j) + rho(n)%f(lb%nx,j))
            rhojp = 2./(rho(n)%f(lb%nx,j+1) + rho(n)%f(lb%nx,j))
            rhojm = 2./(rho(n)%f(lb%nx,j-1) + rho(n)%f(lb%nx,j))
            values(j) = -1.0*(2.*rhoip + rhoim + rhojp + rhojm)
          end do
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        elseif (lb%right == 'neumann') then
          do j = 1,np
            rhoip = 2./(rho(n)%f(lb%nx+1,j) + rho(n)%f(lb%nx,j))
            rhoim = 2./(rho(n)%f(lb%nx-1,j) + rho(n)%f(lb%nx,j))
            rhojp = 2./(rho(n)%f(lb%nx,j+1) + rho(n)%f(lb%nx,j))
            rhojm = 2./(rho(n)%f(lb%nx,j-1) + rho(n)%f(lb%nx,j))
            values(j) = -1.0*(rhoim + rhojp + rhojm)
          end do
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        elseif (lb%right == 'periodic' .or. lb%right == 'internal' .or. lb%right == 'halo') then
          ! do nothing
        else
          print *, 'WRONG RIGTH BC HYPRE'
        end if

        deallocate(values, stencil_indices)
      end if

      if (periodic(2) == 0) then
        ! Values below the box
        ilower = [lb%ilower(1), lb%ilower(2)]
        iupper = [lb%iupper(1), lb%ilower(2)]
        np = lb%iupper(1) - lb%ilower(1) + 1
        allocate(values(np))
        allocate(stencil_indices(1))
        ! If the bottom border is not internal or halo set to zero the stencil 3
        if (lb%bottom /= 'internal' .and. lb%bottom /= 'halo') then
          values = 0.0
          stencil_indices(1) = 3
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        endif

        ! Then modify the central point of the stencil
        stencil_indices(1) = 0
        if (lb%bottom == 'dirichlet') then
          do i = 1,np
            rhoip = 2./(rho(n)%f(i+1,1) + rho(n)%f(i,1))
            rhoim = 2./(rho(n)%f(i-1,1) + rho(n)%f(i,1))
            rhojp = 2./(rho(n)%f(i,2) + rho(n)%f(i,1))
            rhojm = 2./(rho(n)%f(i,0) + rho(n)%f(i,1))
            values(i) = -1.0*(rhoip + rhoim + rhojp + 2.*rhojm)
          end do
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        elseif (lb%bottom == 'neumann') then
          do i = 1,np
            rhoip = 2./(rho(n)%f(i+1,1) + rho(n)%f(i,1))
            rhoim = 2./(rho(n)%f(i-1,1) + rho(n)%f(i,1))
            rhojp = 2./(rho(n)%f(i,2) + rho(n)%f(i,1))
            rhojm = 2./(rho(n)%f(i,0) + rho(n)%f(i,1))
            values(i) = -1.0*(rhoip + rhoim + rhojp)
          end do
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        elseif (lb%bottom == 'periodic' .or. lb%bottom == 'internal' .or. lb%bottom == 'halo') then
          ! do nothing
        else
          print *, 'WRONG BOTTOM BC HYPRE'
        end if

        ! Values above the box
        ilower = [lb%ilower(1), lb%iupper(2)]
        iupper = [lb%iupper(1), lb%iupper(2)]
        ! If the top border is not internal set to zero the stencil 4
        if (lb%top /= 'internal' .and. lb%top /= 'halo') then
          values = 0.0
          stencil_indices(1) = 4
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        endif

        ! Then modify the central point of the stencil
        stencil_indices(1) = 0
        if (lb%top == 'dirichlet') then
          do i = 1,np
            rhoip = 2./(rho(n)%f(i+1,lb%ny) + rho(n)%f(i,lb%ny))
            rhoim = 2./(rho(n)%f(i-1,lb%ny) + rho(n)%f(i,lb%ny))
            rhojp = 2./(rho(n)%f(i,lb%ny+1) + rho(n)%f(i,lb%ny))
            rhojm = 2./(rho(n)%f(i,lb%ny-1) + rho(n)%f(i,lb%ny))
            values(i) = -1.0*(rhoip + rhoim + 2.*rhojp + rhojm)
          end do
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        elseif (lb%top == 'neumann') then
          do i = 1,np
            rhoip = 2./(rho(n)%f(i+1,lb%ny) + rho(n)%f(i,lb%ny))
            rhoim = 2./(rho(n)%f(i-1,lb%ny) + rho(n)%f(i,lb%ny))
            rhojp = 2./(rho(n)%f(i,lb%ny+1) + rho(n)%f(i,lb%ny))
            rhojm = 2./(rho(n)%f(i,lb%ny-1) + rho(n)%f(i,lb%ny))
            values(i) = -1.0*(rhoip + rhoim + rhojm)
          end do
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        elseif (lb%top == 'periodic' .or. lb%top == 'internal' .or. lb%top == 'halo') then
          ! do nothing
        else
          print *, 'WRONG TOP BC HYPRE'
        end if

        
        deallocate(values, stencil_indices)
      end if

      ! Corners:
      allocate(values(1))
      allocate(stencil_indices(1))
      stencil_indices(1) = 0

      ilower = [lb%ilower(1), lb%ilower(2)]
      iupper = [lb%ilower(1), lb%ilower(2)]
      rhoip = 2./(rho(n)%f(2,1) + rho(n)%f(1,1))
      rhoim = 2./(rho(n)%f(0,1) + rho(n)%f(1,1))
      rhojp = 2./(rho(n)%f(1,2) + rho(n)%f(1,1))
      rhojm = 2./(rho(n)%f(1,0) + rho(n)%f(1,1))
      if (lb%left == 'dirichlet' .and. lb%bottom == 'dirichlet') then
        values(1) = -1.0*(rhoip + 2.*rhoim + rhojp + 2.*rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%left == 'neumann' .and. lb%bottom == 'neumann') then
        values(1) = -1.0*(rhoip + rhojp)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%left == 'dirichlet' .and. (lb%bottom == 'internal' .or. lb%bottom == 'halo')) then
        values(1) = -1.0*(rhoip + 2.*rhoim + rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif ((lb%left == 'internal' .or. lb%left == 'halo') .and. lb%bottom == 'dirichlet') then
        values(1) = -1.0*(rhoip + rhoim + rhojp + 2.*rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%left == 'neumann' .and. (lb%bottom == 'internal' .or. lb%bottom == 'halo')) then
        values(1) = -1.0*(rhoip + rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif ((lb%left == 'internal' .or. lb%left == 'halo') .and. lb%bottom == 'neumann') then
        values(1) = -1.0*(rhoip + rhoim + rhojp)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      else
        ! do nothing
      end if

      ilower = [lb%ilower(1), lb%iupper(2)]
      iupper = [lb%ilower(1), lb%iupper(2)]
      rhoip = 2./(rho(n)%f(2,lb%ny) + rho(n)%f(1,lb%ny))
      rhoim = 2./(rho(n)%f(0,lb%ny) + rho(n)%f(1,lb%ny))
      rhojp = 2./(rho(n)%f(1,lb%ny+1) + rho(n)%f(1,lb%ny))
      rhojm = 2./(rho(n)%f(1,lb%ny-1) + rho(n)%f(1,lb%ny))
      if (lb%left == 'dirichlet' .and. lb%top == 'dirichlet') then
        values(1) = -1.0*(rhoip + 2.*rhoim + 2.*rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)    
      elseif (lb%left == 'neumann' .and. lb%top == 'neumann') then
        values(1) = -1.0*(rhoip + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)    
      elseif (lb%left == 'dirichlet' .and. (lb%top == 'internal' .or. lb%top == 'halo')) then
        values(1) = -1.0*(rhoip + 2.*rhoim + rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)    
      elseif ((lb%left == 'internal'.or. lb%left == 'halo') .and. lb%top == 'dirichlet' ) then
        values(1) = -1.0*(rhoip + rhoim + 2.*rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)    
      elseif (lb%left == 'neumann' .and. (lb%top == 'internal' .or. lb%top == 'halo')) then
        values(1) = -1.0*(rhoip + rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)    
      elseif ((lb%left == 'internal'.or. lb%left == 'halo') .and. lb%top == 'neumann' ) then
        values(1) = -1.0*(rhoip + rhoim + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)    
      else
        ! do nothing
      end if
      
      ilower = [lb%iupper(1), lb%ilower(2)]
      iupper = [lb%iupper(1), lb%ilower(2)]
      rhoip = 2./(rho(n)%f(lb%nx+1,1) + rho(n)%f(lb%nx,1))
      rhoim = 2./(rho(n)%f(lb%nx-1,1) + rho(n)%f(lb%nx,1))
      rhojp = 2./(rho(n)%f(lb%nx,2) + rho(n)%f(lb%nx,1))
      rhojm = 2./(rho(n)%f(lb%nx,0) + rho(n)%f(lb%nx,1))
      if (lb%right == 'dirichlet' .and. lb%bottom == 'dirichlet') then
        values(1) = -1.0*(2.*rhoip + rhoim + rhojp + 2.*rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%right == 'neumann' .and. lb%bottom == 'neumann') then
        values(1) = -1.0*(rhoim + rhojp)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%right == 'dirichlet' .and. (lb%bottom == 'internal' .or. lb%bottom == 'halo')) then
        values(1) = -1.0*(2.*rhoip + rhoim + rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif ((lb%right == 'internal' .or. lb%right == 'halo') .and. lb%bottom == 'dirichlet') then
        values(1) = -1.0*(rhoip + rhoim + rhojp + 2.*rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%right == 'neumann' .and. (lb%bottom == 'internal' .or. lb%bottom == 'halo')) then
        values(1) = -1.0*(rhoim + rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif ((lb%right == 'internal' .or. lb%right == 'halo') .and. lb%bottom == 'neumann') then
        values(1) = -1.0*(rhoip + rhoim + rhojp)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      else
        ! do nothing
      end if

      ilower = [lb%iupper(1), lb%iupper(2)]
      iupper = [lb%iupper(1), lb%iupper(2)]
      rhoip = 2./(rho(n)%f(lb%nx+1,lb%ny) + rho(n)%f(lb%nx,lb%ny))
      rhoim = 2./(rho(n)%f(lb%nx-1,lb%ny) + rho(n)%f(lb%nx,lb%ny))
      rhojp = 2./(rho(n)%f(lb%nx,lb%ny+1) + rho(n)%f(lb%nx,lb%ny))
      rhojm = 2./(rho(n)%f(lb%nx,lb%ny-1) + rho(n)%f(lb%nx,lb%ny))
      if (lb%right == 'dirichlet' .and. lb%top == 'dirichlet') then
        values(1) = -1.0*(2.*rhoip + rhoim + 2.*rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%right == 'neumann' .and. lb%top == 'neumann') then
        values(1) = -1.0*(rhoim + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%right == 'dirichlet' .and. (lb%top == 'internal' .or. lb%top == 'halo')) then
        values(1) = -1.0*(2.*rhoip + rhoim + rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif ((lb%right == 'internal' .or. lb%right == 'halo') .and. lb%top == 'dirichlet') then
        values(1) = -1.0*(rhoip + rhoim + 2.*rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif (lb%right == 'neumann' .and. (lb%top == 'internal' .or. lb%top == 'halo')) then
        values(1) = -1.0*(rhoim + rhojp + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      elseif ((lb%right == 'internal' .or. lb%right == 'halo') .and. lb%top == 'neumann') then
        values(1) = -1.0*(rhoip + rhoim + rhojm)
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
      else
        ! do nothing
      end if
 
      deallocate(values,stencil_indices)
    end do box_cycle
    
    ! This is a collective call finalizing the matrix assembly.
    ! The matrix is now ``ready to be used''
    call HYPRE_StructMatrixAssemble(A, ierr)

    ! To print the matrix A, uncomment the following line
    !call HYPRE_StructMatrixPrint(A, 0, ierr)

  end subroutine update_coefficient_matrix
  !===============================================================================
  
end module multiphase
