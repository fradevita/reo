module hypre_solver

  use grid_2d

  ! Define the grid, the matrix A, the solution array b, the constant array d
  ! the solver and the stencil
  integer(8) :: grid, A, d, b, solver
  integer(8) stencil

  ! Define some properties of the solver
  integer :: verbose = 0, periodic(2) = [0, 0]
  real :: tolerance = 1.0e-6

  private
  public :: init_hypre_solver, solve_poisson_hypre, destroy_hypre_solver, &
            tolerance, verbose, periodic
  
contains
  
  subroutine init_hypre_solver(left,right,top,bottom)

    implicit none

    ! Type of BC
    !character(len=9), intent(in) :: left, right, top, bottom
    character(len=*), intent(in) :: left, right, top, bottom

    ! Local variables
    integer :: entry, i, j, nentries, nvalues, n
    integer, dimension(2) :: ilower, iupper
    integer, dimension(5,2) :: offsets
    integer, dimension(:), allocatable :: stencil_indices
    real :: delta2
    real, dimension(:), allocatable :: values

    !character*32 :: matfile
    
    ! **** 1. Create the grid object ****
    call HYPRE_StructGridCreate(mpi_common_world, 2, grid, ierr)
    call HYPRE_StructGridSetPeriodic(grid, periodic, ierr)

    ! Cycle over all the boxes in the grid

    ! Set grid extents for the box
    ilower = [1, 1]
    iupper = [nx, ny]
    call HYPRE_StructGridSetExtents(grid, ilower, iupper, ierr)

    ! Assemble the grid
    call HYPRE_StructGridAssemble(grid, ierr)

    ! **** 2. Define the disretization stencil ****

    ! Create a 2D 5-pt stencil object
    call HYPRE_StructStencilCreate(2, 5, stencil, ierr)

    ! Define the geometry of the stencil
    offsets = transpose(reshape([0, 0, -1, 0, 1, 0, 0, -1, 0, 1], &
                          shape(transpose(offsets))))
    do entry = 1,5
      call HYPRE_StructStencilSetElement(stencil, entry - 1, offsets(entry,:), ierr) 
    end do

    ! **** 3. Set up a Struct Matrix ****

    ! Create an empty matrix object
    call HYPRE_StructMatrixCreate(mpi_common_world, grid, stencil, A, ierr)

    ! Indicate that the matrix coefficients are ready to be set
    call HYPRE_StructMatrixInitialize(A, ierr)

    ! Set the matrix cefficients. We first set the same stencil entries
    ! for each grid point. Then we make modifications to grid points
    ! near the boundary.

    ! nx*ny grid points each with 5 stencil entries
    ilower = [1, 1]
    iupper = [nx, ny]
    nentries = 5
    nvalues = nx*ny*nentries
    allocate(values(nvalues))
    allocate(stencil_indices(nentries))
    stencil_indices = [0, 1, 2, 3, 4]
    delta2 = dx**2
    do i = 1,nvalues,nentries
      values(i) = -4.0 / delta2
      do j = 1,nentries-1
        values(i+j) = 1.0 / delta2
      end do
    end do
    call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, &
                                       stencil_indices, values, ierr)
    deallocate(values, stencil_indices)
    
    ! Set the coefficients outside of the domain and on the boundaries
    if (periodic(1) == 0) then
      ! Values on the left of the box
      ilower = [1, 1]
      iupper = [1, ny]
      allocate(values(ny))
      allocate(stencil_indices(1))
      do i = 1,ny
        values(i) = 0.0
      end do
      stencil_indices(1) = 1
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
           stencil_indices, values, ierr)

      ! Dirichlet
      if (left == 'dirichlet') then
        do i = 1,ny
          values(i) = -5.0 / delta2
        end do
      elseif (left == 'neumann') then 
        do i = 1,ny
          values(i) = -3.0 / delta2
        end do
      elseif (left == 'periodic') then
       ! do nothing
      else
        print *, 'WRONG LEFT BC HYPRE'
      end if

      stencil_indices(1) = 0
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
                                          stencil_indices, values, ierr)

      ! Values on the right of the box
      ilower = [nx, 1]
      iupper = [nx, ny]
      do i = 1,ny
        values(i) = 0.0
      end do
      stencil_indices(1) = 2
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
                                          stencil_indices, values, ierr)
      ! Dirichlet
      if (right == 'dirichlet') then
        do i = 1,ny
          values(i) = -5.0 / delta2
        end do
      elseif (right == 'neumann') then
        do i = 1,ny
          values(i) = -3.0 / delta2
        end do
      elseif (right == 'periodic') then
        ! do nothing
      else
        print *, 'WRONG RIGTH BC HYPRE'
      end if

      stencil_indices(1) = 0
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
                                          stencil_indices, values, ierr)
      deallocate(values, stencil_indices)
    end if

    if (periodic(2) == 0) then
      ! Values below the box
      ilower = [1, 1]
      iupper = [nx, 1]
      allocate(values(nx))
      allocate(stencil_indices(1))
      stencil_indices(1) = 3
      do i = 1,nx
        values(i) = 0.0
      end do
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
           stencil_indices, values, ierr)

      ! Dirichlet
      if (bottom == 'dirichlet') then
        do i = 1,nx
          values(i) = -5.0 / delta2
        end do
      elseif (bottom == 'neumann') then
        do i = 1,nx
          values(i) = -3.0 / delta2
        end do
      elseif (bottom == 'periodic') then
        ! do nothing
      else
        print *, 'WRONG BOTTOM BC HYPRE'
      end if

      stencil_indices(1) = 0
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
           stencil_indices, values, ierr)

      ! Values above the box
      ilower = [1, ny]
      iupper = [nx, ny]
      stencil_indices(1) = 4
      do i = 1,nx
        values(i) = 0.0
      end do
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
           stencil_indices, values, ierr)

      ! Dirichlet
      if (top == 'dirichlet') then
        do i = 1,nx
          values(i) = -5.0 / delta2
        end do
      elseif (top == 'neumann') then
        do i = 1,nx
          values(i) = -3.0 / delta2
        end do
      elseif (top == 'periodic') then
        ! do nothing
      else
        print *, 'WRONG TOP BC HYPRE'
      end if

      stencil_indices(1) = 0
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
           stencil_indices, values, ierr)
      deallocate(values, stencil_indices)
    end if

    ! Corners
    ilower = [1, 1]
    iupper = [1, 1]
    allocate(values(1))
    allocate(stencil_indices(1))
    stencil_indices(1) = 0
    if (left == 'dirichlet' .and. bottom == 'dirichlet') then
      values(1) = -6.0 / delta2
    elseif (left == 'neumann' .and. bottom == 'neumann') then
      values(1) = -2.0 / delta2
    else
      values(1) = -4.0 / delta2
    end if
    call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
        stencil_indices, values, ierr)    
    ilower = [1, ny]
    iupper = [1, ny]
    stencil_indices(1) = 0
    if (left == 'dirichlet' .and. top == 'dirichlet') then
      values(1) = -6.0 / delta2
    elseif (left == 'neumann' .and. top == 'neumann') then
      values(1) = -2.0 / delta2
    else
      values(1) = -4.0 / delta2
    end if
    call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
        stencil_indices, values, ierr)    
    ilower = [nx, 1]
    iupper = [nx, 1]
    stencil_indices(1) = 0
    if (right == 'dirichlet' .and. bottom == 'dirichlet') then
      values(1) = -6.0 / delta2
    elseif (right == 'neumann' .and. bottom == 'neumann') then
      values(1) = -2.0 / delta2
    else
      values(1) = -4.0 / delta2
    end if
    call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
        stencil_indices, values, ierr)    
    ilower = [nx, ny]
    iupper = [nx, ny]
    stencil_indices(1) = 0
    if (right == 'dirichlet' .and. top == 'dirichlet') then
      values(1) = -6.0 / delta2
    elseif (right == 'neumann' .and. top == 'neumann') then
      values(1) = -2.0 / delta2
    else
      values(1) = -4.0 / delta2
    end if
    call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, & 
        stencil_indices, values, ierr)
    deallocate(values,stencil_indices)
    
    ! This is a collective call finalizing the matrix assembly.
    ! The matrix is now ``ready to be used''
    call HYPRE_StructMatrixAssemble(A, ierr)

!    matfile = 'matrix'
!    matfile(7:7) = char(0)
!    call HYPRE_StructMatrixPrint(matfile, A, 0, ierr)

    ! Create an empty vector object
    call HYPRE_StructVectorCreate(mpi_common_world, grid, b, ierr)
    call HYPRE_StructVectorCreate(mpi_common_world, grid, d, ierr)

    ! Indicate that the vector coefficients are ready to be set
    call HYPRE_StructVectorInitialize(b, ierr)
    call HYPRE_StructVectorInitialize(d, ierr)

    ! Create a PCG solver
   call HYPRE_StructPCGCreate(mpi_common_world, solver, ierr)
   call HYPRE_StructPCGSetTol(solver, tolerance, ierr) ! convergence tolerance
   call HYPRE_StructPCGSetPrintLevel(solver, verbose, ierr) ! amount of info. printed

    ! call HYPRE_StructSMGCreate(mpi_common_world, solver, ierr)
    ! call HYPRE_StructSMGSetMemoryUse(solver, 0, ierr)
    ! call HYPRE_StructSMGSetMaxIter(solver,50,ierr)
    ! call HYPRE_StructSMGSetRelChange(solver, 0, ierr)
    ! call HYPRE_StructSMGSetNumPreRelax(solver, 1, ierr)
    ! call HYPRE_StructSMGSetNumPostRelax(solver, 1, ierr)
    ! call HYPRE_StructSMGSetLogging(solver,1, ierr)

  end subroutine init_hypre_solver
  
  !subroutine solve_poisson(rhs)
  subroutine solve_poisson_hypre(temp)
      
    ! This subroutine solve the poisson equation
    ! nabal^2 f = rhs
    ! with right-hand side rhs. 
    ! RHS is an array of dimension nx*ny and at the end of the subroutine
    ! it will contain the solution array f.
    ! We solve the equation in the form Ab = d
    ! so d is equal to rhs and b is the solution array f.

    use boundary_conditions
    
    implicit none

    integer :: i, j, ilower(2), iupper(2),  num_iterations
    real :: mean
    real, dimension(nx*ny) :: rhs
    real, intent(inout) :: temp(:,:)

    ! We copy the array rhs in the two-dimensional array phi
    do j = 1,ny
      do i = 1,nx
        rhs(i + (j-1)*nx) = temp(i,j)
      end do
    end do

    ! Set the array d equal to the rhs
    ilower = [1, 1]
    iupper = [nx, ny]
    call HYPRE_StructVectorSetBoxValues(d, ilower, iupper, rhs, ierr)
    
    ! Set the solution array b to zero
    rhs = 0.0
    call HYPRE_StructVectorSetBoxValues(b, ilower, iupper, rhs, ierr)

    ! Solve the poisson equation
    call HYPRE_StructPCGSetup(solver, A, d, b, ierr)
    call HYPRE_StructPCGSolve(solver, A, d, b, ierr)

    !call HYPRE_StructSMGSetup(solver, A, d, b, ierr)
    !call HYPRE_StructSMGSolve(solver, A, d, b, ierr)
    !call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
   
    ! Save the solution of the Poisson equation inside the rhs array
    call HYPRE_StructVectorGetBoxValues(b, ilower, iupper, rhs, ierr)
    
    if (periodic(1) /= 0 .and. periodic(2) /= 0) then
      ! Normalize the solution
      mean = 0.0
      do i = 1,nx*ny
        mean = mean + rhs(i)
      end do
      mean = mean / (nx*ny)
      
      do i = 1,nx*ny
        rhs(i) = rhs(i) - mean
      end do
    end if

    ! We copy the array rhs in the two-dimensional array phi
    do j = 1,ny
      do i = 1,nx
        temp(i,j) = rhs(i + (j-1)*nx)
      end do
    end do
        
  end subroutine solve_poisson_hypre

  subroutine destroy_hypre_solver

    implicit none
    
    call HYPRE_StructGridDestroy(grid, ierr)
    call HYPRE_StructStencilDestroy(stencil, ierr)
    call HYPRE_StructMatrixDestroy(A, ierr)
    call HYPRE_StructVectorDestroy(b, ierr)
    call HYPRE_StructVectorDestroy(d, ierr)
    call HYPRE_StructPCGDestroy(solver, ierr)
    !call HYPRE_StructSMGDestroy(solver,ierr)
  end subroutine destroy_hypre_solver

end module hypre_solver
