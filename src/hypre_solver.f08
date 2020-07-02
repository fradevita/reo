module hypre_solver

  use grid_2d, only : box, boxarray, nbox, ierr, mpi_common_world, myarray, dp

  ! Define the grid, the matrix A, the solution array b, the constant array d
  ! the solver and the stencil
  integer(8) :: grid, A, d, b, solver, precond
  integer(8) :: stencil

  ! Define some properties of the solver
  integer :: verbose = 0, periodic(2) = [0, 0]
  real(kind=dp) :: tolerance = 1.0e-6

  type(box) :: lb

  type(myarray), dimension(:), allocatable :: prhs

  character(len=10) :: type_hypre_solver = 'BiCGSTAB'

  private
  public :: init_hypre_solver, solve_poisson_hypre, destroy_hypre_solver, &
            tolerance, verbose, periodic, A, type_hypre_solver
  
contains

  !===============================================================================
  subroutine init_hypre_solver()

    implicit none

    ! Local variables
    integer :: entry, i, j, nentries, nvalues, n, np, nx, ny
    integer, dimension(2) :: ilower, iupper
    integer, dimension(5,2) :: offsets
    integer, dimension(:), allocatable :: stencil_indices
    real(kind=dp), dimension(:), allocatable :: values
 
    ! **** 1. Create the grid object ****
    call HYPRE_StructGridCreate(mpi_common_world, 2, grid, ierr)
    call HYPRE_StructGridSetPeriodic(grid, periodic, ierr)

    ! Create the grid for each box
    do n = 1, nbox
      ! Select the n-th box
      lb = boxarray(n)
      ilower = lb%ilower
      iupper = lb%iupper

      ! Set grid extents for the box
      call HYPRE_StructGridSetExtents(grid, ilower, iupper, ierr)
    end do

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

    ! Set the matrix cefficients for each box.
    box_cycle : do n = 1,nbox

      ! Select the box
      lb = boxarray(n)
      ilower = lb%ilower
      iupper = lb%iupper

      ! We first set the same stencil entries
      ! for each grid point. Then we make modifications to grid points
      ! near the boundary.

      ! nx*ny grid points each with 5 stencil entries
      nentries = 5
      nvalues = lb%nx*lb%ny*nentries
      allocate(values(nvalues))
      allocate(stencil_indices(nentries))
      stencil_indices = [0, 1, 2, 3, 4]
      do i = 1,nvalues,nentries
        values(i) = -4.0
        do j = 1,nentries-1
          values(i+j) = 1.0
        end do
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
        ! If the left border is not internal set to zero the stencil 1
        if (lb%left /= 'internal' .and. lb%left /= 'halo') then          
          values = 0.0
          stencil_indices(1) = 1
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        endif

        ! Then modify the central point of the stencil
        if (lb%left == 'dirichlet') then
          values = -5.0
        elseif (lb%left == 'neumann') then 
          values = -3.0
        elseif (lb%left == 'periodic' .or. lb%left == 'internal' .or. lb%left == 'halo') then
          values = -4.0
        else
          print *, 'WRONG LEFT BC HYPRE'
        end if

        stencil_indices(1) = 0
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

        ! Values on the right of the box
        ilower = [lb%iupper(1), lb%ilower(2)]
        iupper = [lb%iupper(1), lb%iupper(2)]
        
        ! If the right border is not internal set to zero the stencil 2
        if (lb%right /= 'internal' .and. lb%right /= 'halo') then
          values = 0.0
          stencil_indices(1) = 2
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        endif

        ! Then modify the central point of the stencil
        if (lb%right == 'dirichlet') then
          values = -5.0
        elseif (lb%right == 'neumann') then
          values = -3.0
        elseif (lb%right == 'periodic' .or. lb%right == 'internal' .or. lb%right == 'halo') then
          values = -4.0
        else
          print *, 'WRONG RIGTH BC HYPRE'
        end if

        stencil_indices(1) = 0
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

        deallocate(values, stencil_indices)
      end if

      if (periodic(2) == 0) then
        ! Values below the box
        ilower = [lb%ilower(1), lb%ilower(2)]
        iupper = [lb%iupper(1), lb%ilower(2)]
        np = lb%iupper(1) - lb%ilower(1) + 1
        allocate(values(np))
        allocate(stencil_indices(1))
        ! If the bottom border is not internal set to zero the stencil 3
        if (lb%bottom /= 'internal' .and. lb%bottom /= 'halo') then
          values = 0.0
          stencil_indices(1) = 3
          call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        endif

        ! Then modify the central point of the stencil
        if (lb%bottom == 'dirichlet') then
          values = -5.0
        elseif (lb%bottom == 'neumann') then
          values = -3.0
        elseif (lb%bottom == 'periodic' .or. lb%bottom == 'internal' .or. lb%bottom == 'halo') then
          values = -4.0
        else
          print *, 'WRONG BOTTOM BC HYPRE'
        end if

        stencil_indices(1) = 0
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

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
        if (lb%top == 'dirichlet') then
          values = -5.0
        elseif (lb%top == 'neumann') then
          values = -3.0
        elseif (lb%top == 'periodic' .or. lb%top == 'internal' .or. lb%top == 'halo') then
          values = -4.0
        else
          print *, 'WRONG TOP BC HYPRE'
        end if

        stencil_indices(1) = 0
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
        
        deallocate(values, stencil_indices)
      end if

      ! Corners:
      allocate(values(1))
      allocate(stencil_indices(1))
      stencil_indices(1) = 0

      ilower = [lb%ilower(1), lb%ilower(2)]
      iupper = [lb%ilower(1), lb%ilower(2)]
      if (lb%left == 'dirichlet' .and. lb%bottom == 'dirichlet') then
        values(1) = -6.0
      elseif (lb%left == 'neumann' .and. lb%bottom == 'neumann') then
        values(1) = -2.0
      elseif (lb%left == 'dirichlet' .and. (lb%bottom == 'internal' .or. lb%bottom == 'halo') .or. &
              (lb%left == 'internal' .or. lb%left == 'halo') .and. lb%bottom == 'dirichlet' ) then
        values(1) = -5.0
      elseif (lb%left == 'neumann' .and. (lb%bottom == 'internal' .or. lb%bottom == 'halo') .or. &
              (lb%left == 'internal' .or. lb%left == 'halo') .and. lb%bottom == 'neumann' ) then
        values(1) = -3.0
      else
        values(1) = -4.0
      end if
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

      ilower = [lb%ilower(1), lb%iupper(2)]
      iupper = [lb%ilower(1), lb%iupper(2)]
      if (lb%left == 'dirichlet' .and. lb%top == 'dirichlet') then
        values(1) = -6.0
      elseif (lb%left == 'neumann' .and. lb%top == 'neumann') then
        values(1) = -2.0
      elseif (lb%left == 'dirichlet' .and. (lb%top == 'internal' .or. lb%top == 'halo') .or. &
              (lb%left == 'internal'.or. lb%left == 'halo') .and. lb%top == 'dirichlet' ) then
        values(1) = -5.0
      elseif (lb%left == 'neumann' .and. (lb%top == 'internal' .or. lb%top == 'halo') .or. &
              (lb%left == 'internal'.or. lb%left == 'halo') .and. lb%top == 'neumann' ) then
        values(1) = -3.0
      else
        values(1) = -4.0
      end if
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)    
      
      ilower = [lb%iupper(1), lb%ilower(2)]
      iupper = [lb%iupper(1), lb%ilower(2)]
      if (lb%right == 'dirichlet' .and. lb%bottom == 'dirichlet') then
        values(1) = -6.0
      elseif (lb%right == 'neumann' .and. lb%bottom == 'neumann') then
        values(1) = -2.0
      elseif (lb%right == 'dirichlet' .and. (lb%bottom == 'internal' .or. lb%bottom == 'halo') .or. &
              (lb%right == 'internal' .or. lb%right == 'halo') .and. lb%bottom == 'dirichlet' ) then
        values(1) = -5.0
      elseif (lb%right == 'neumann' .and. (lb%bottom == 'internal' .or. lb%bottom == 'halo') .or. &
              (lb%right == 'internal' .or. lb%right == 'halo') .and. lb%bottom == 'neumann' ) then
        values(1) = -3.0
      else
        values(1) = -4.0
      end if
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

      ilower = [lb%iupper(1), lb%iupper(2)]
      iupper = [lb%iupper(1), lb%iupper(2)]
      if (lb%right == 'dirichlet' .and. lb%top == 'dirichlet') then
        values(1) = -6.0
      elseif (lb%right == 'neumann' .and. lb%top == 'neumann') then
        values(1) = -2.0
      elseif (lb%right == 'dirichlet' .and. (lb%top == 'internal' .or. lb%top == 'halo') .or. &
              (lb%right == 'internal' .or. lb%right == 'halo') .and. lb%top == 'dirichlet' ) then
        values(1) = -5.0
      elseif (lb%right == 'neumann' .and. (lb%top == 'internal' .or. lb%top == 'halo') .or. &
              (lb%right == 'internal' .or. lb%right == 'halo') .and. lb%top == 'neumann' ) then
        values(1) = -3.0
      else
        values(1) = -4.0
      end if
      call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
 
      deallocate(values,stencil_indices)
    end do box_cycle
    
    ! This is a collective call finalizing the matrix assembly.
    ! The matrix is now ``ready to be used''
    call HYPRE_StructMatrixAssemble(A, ierr)

    ! To print the matrix A, uncomment the following line
    ! call HYPRE_StructMatrixPrint(A, 0, ierr)

    ! Create an empty vector object
    call HYPRE_StructVectorCreate(mpi_common_world, grid, b, ierr)
    call HYPRE_StructVectorCreate(mpi_common_world, grid, d, ierr)

    ! Indicate that the vector coefficients are ready to be set
    call HYPRE_StructVectorInitialize(b, ierr)
    call HYPRE_StructVectorInitialize(d, ierr)

    ! Allocate the Poisson RHS variable
    allocate(prhs(nbox))
    do n = 1,nbox
      nx = boxarray(n)%iupper(1) - boxarray(n)%ilower(1) + 1
      ny = boxarray(n)%iupper(2) - boxarray(n)%ilower(2) + 1
      allocate(prhs(n)%f(nx,ny))
    end do

  end subroutine init_hypre_solver
  !===============================================================================
  
  !===============================================================================
  subroutine solve_poisson_hypre(temp)
      
    ! This subroutine solve the poisson equation
    ! nabal^2 f = rhs
    ! with right-hand side rhs. 
    ! RHS is an array of dimension nx*ny and at the end of the subroutine
    ! it will contain the solution array f.
    ! We solve the equation in the form Ab = d
    ! so d is equal to rhs and b is the solution array f.

    use grid_2d, only : box

    implicit none

    ! In / out variables
    type(myarray), intent(inout) :: temp (:)

    ! Local variabels
    integer :: i, j, ilower(2), iupper(2),  num_iterations, n
    real(kind=dp) :: mean, meanb, delta2
    real(kind=dp), dimension(:), allocatable :: rhs
    type(box) :: lb

    ! Cycle over boxes
    do n = 1,nbox

      delta2 = boxarray(n)%delta**2
 
      lb = boxarray(n)
      allocate(rhs(lb%nx*lb%ny))

      !  Copy the two-dimensional array temp in the one-dimensional array rhs
      do j = 1,lb%ny
        do i = 1,lb%nx
          rhs(i + (j-1)*lb%nx) = temp(n)%f(i,j)*delta2
        end do
      end do

      ! Set the array d equal to the rhs
      ilower = lb%ilower
      iupper = lb%iupper
      call HYPRE_StructVectorSetBoxValues(d, ilower, iupper, rhs, ierr)

      ! Set the solution array b to zero
      rhs = 0.0
      call HYPRE_StructVectorSetBoxValues(b, ilower, iupper, rhs, ierr)
      deallocate(rhs)

    end do

    ! Assemble the vectors
    call HYPRE_StructVectorAssemble(d, ierr)
    call HYPRE_StructVectorAssemble(b, ierr)

    ! Solve the poisson equation
    if (type_hypre_solver == 'BiCGSTAB') then
      call HYPRE_StructBiCGSTABCreate(mpi_common_world, solver, ierr)
      call HYPRE_StructJacobiCreate(mpi_common_world, precond, ierr)
      call HYPRE_StructJacobiSetMaxIter(precond, 2, ierr)
      call HYPRE_StructJacobiSetTol(precond, 0.0, ierr)
      call HYPRE_StructJacobiSetZeroGuess(precond, ierr)
      call HYPRE_StructBiCGSTABSetup(solver, A, d, b, ierr)
      call HYPRE_StructBiCGSTABSolve(solver, A, d, b, ierr)
      call HYPRE_StructBiCGSTABSetPrecond(solver, 1, precond, ierr)
    elseif (type_hypre_solver == 'SGM') then
      call HYPRE_StructSMGCreate(mpi_common_world, solver, ierr)
      call HYPRE_StructSMGSetMemoryUse(solver, 0, ierr)
      call HYPRE_StructSMGSetMaxIter(solver,100,ierr)
      call HYPRE_StructSMGSetRelChange(solver, 0, ierr)
      call HYPRE_StructSMGSetNumPreRelax(solver, 1, ierr)
      call HYPRE_StructSMGSetNumPostRelax(solver, 1, ierr)
      call HYPRE_StructSMGSetLogging(solver,1, ierr)
      call HYPRE_StructSMGSetup(solver, A, d, b, ierr)
      call HYPRE_StructSMGSolve(solver, A, d, b, ierr)
      call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
    else  
      print *, 'Wrong Hypre Solver'
    endif 

    ! Save the solution of the Poisson equation inside the temp array
    do n = 1,nbox

      lb = boxarray(n)

      allocate(rhs(lb%nx*lb%ny))
      ilower = lb%ilower
      iupper = lb%iupper

      call HYPRE_StructVectorGetBoxValues(b, ilower, iupper, rhs, ierr)
      
      ! We copy the array rhs in the two-dimensional array phi
      do j = 1,lb%ny
        do i = 1,lb%nx
          temp(n)%f(i,j) = rhs(i + (j - 1)*lb%nx)
        end do
      end do
      deallocate(rhs)

    end do

    ! In the case of a periodic domain in both direction substract the mean      
    if (periodic(1) /= 0 .and. periodic(2) /= 0) then
      ! Normalize the solution
      mean = 0.0
      do n = 1,nbox
        meanb = 0.0
        lb = boxarray(n)
        do j = 1,lb%ny
          do i = 1,lb%nx
            meanb = meanb + temp(n)%f(i,j)
          end do
        end do
        mean = mean + meanb / (lb%nx*lb%ny)
      end do

      do n = 1,nbox
        lb = boxarray(n)
        do j = 1,lb%ny
          do i = 1,lb%nx
            temp(n)%f(i,j) = temp(n)%f(i,j) - mean
          end do
        end do
      end do
    end if

    if (type_hypre_solver == 'BiCGSTAB') then
      call HYPRE_StructBiCGSTABDestroy(solver, ierr)
      call HYPRE_StructJacobiDestroy(precond, ierr)
    elseif (type_hypre_solver == 'SGM') then
      call HYPRE_StructSMGDestroy(solver,ierr)
    else  
      print *, 'Wrong Hypre Solver'
    endif 

  end subroutine solve_poisson_hypre
  !===============================================================================

  !===============================================================================
  subroutine destroy_hypre_solver

    implicit none

    integer :: n
    
    call HYPRE_StructGridDestroy(grid, ierr)
    call HYPRE_StructStencilDestroy(stencil, ierr)
    call HYPRE_StructMatrixDestroy(A, ierr)
    call HYPRE_StructVectorDestroy(b, ierr)
    call HYPRE_StructVectorDestroy(d, ierr)
  
    do n = 1,nbox
      deallocate(prhs(n)%f)
    end do
    deallocate(prhs)
  
  end subroutine destroy_hypre_solver
  !===============================================================================

end module hypre_solver
