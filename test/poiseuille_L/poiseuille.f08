program poiseuille

  use grid_2d
  use navier_stokes
  
  implicit none
  include 'mpif.h'

  integer :: j, nx, ny
  real(kind=dp) :: y

  ! Auxiliary field to check steady state
  type(field), dimension(:), allocatable :: uold

  ! First initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
  mpi_common_world = MPI_COMM_WORLD

  ! Set the number of points and the domain size
  nx = 32
  ny = 32

  ! Create the grid
  nbox = 3
  allocate(boxarray(nbox))
  call create_box(1, nx, ny, 1.d0, 1.d0, -0.5d0, -0.5d0)
  call create_box(2, nx, ny, 1.d0, 1.d0,  0.5d0, -0.5d0)
  call create_box(3, nx, ny, 1.d0, 1.d0,  0.5d0,  0.5d0)
  
  ! Boundary conditions
  boxarray(1)%left_boundary = 'inflow'
  boxarray(1)%right_boundary = 'internal'
  boxarray(1)%rb = 2
  boxarray(2)%left_boundary = 'internal'
  boxarray(2)%top_boundary = 'internal'
  boxarray(2)%lb = 1
  boxarray(2)%tb = 3
  boxarray(2)%ilower = [nx+1, 1]
  boxarray(2)%iupper = [2*nx, ny]
  boxarray(3)%bottom_boundary = 'internal'
  boxarray(3)%top_boundary = 'outflow'
  boxarray(3)%bb = 2
  boxarray(3)%ilower = [nx+1, ny+1]
  boxarray(3)%iupper = [2*nx, 2*ny]

  allocate(uold(nbox))
 
  ! Select Poisson solver 
  poisson_solver_type = 'itr'

  ! First we need to initialize the solver
  call init_ns_solver() 

  ! Boundary condition for the left box
  do j = u(1)%lo(2),u(1)%up(2)
    y = boxarray(1)%p0(2) + (j-0.5)*boxarray(1)%delta
    u(1)%l(j) = 0.5*(0.25 - y**2)
  end do

  ! We compute some quantites to check the properties of the scheme
  event_i => e_istep
  event_output => output

  ! We run the simulation
  dt = 0.5*dt
  call solve()

contains
  
  subroutine e_istep()

    implicit none

    uold = u
    
  end subroutine e_istep
    
  subroutine output()

    implicit none
    
    ! Check for steady-state
    real(kind=dp) :: diff, x, y, uc, vc, vort, Qi, Qo
    logical :: steady
    integer :: i, j, nb
    character(len=1) :: sn
    character(len=4) :: filename

    steady = .true.

    do nb = 1,nbox
      do j = u(nb)%lo(2),u(nb)%up(2)
        do i = u(nb)%lo(1),u(nb)%up(1)
          diff = abs(u(nb)%f(i,j) - uold(nb)%f(i,j))
          if (diff > 1.0e-8) steady = .false.
        end do
      end do
    end do

    if (steady) then

      ! Compute the difference between the inflow mass flow and the outflow 
      ! mass flow
      Qi = 0.0
      do j = u(1)%lo(2),u(1)%up(2)
        Qi = Qi + u(1)%f(u(1)%lo(1)-1,j)
      end do
      Qi = Qi / float(u(1)%up(2) - u(1)%lo(2))

      Qo = 0.0
      do i = v(3)%lo(1),v(3)%up(1)
        Qo = Qo + v(3)%f(i,v(3)%up(2))
      end do
      Qo = Qo / float(v(3)%up(1) - v(3)%lo(1))

      print *, Qi, Qo, (Qo - Qi)/Qi

      ! Save the velocity field, pressure and vorticity
      do nb = 1,nbox
        write(sn,'(I1)') nb
        filename = 'box'//sn
        open(10+nb, file = filename)
        do j = p(nb)%lo(2),p(nb)%up(2)
          y = boxarray(nb)%p0(2) + (j - 0.5)*boxarray(nb)%delta
          do i = p(nb)%lo(1),p(nb)%up(1)
            x = boxarray(nb)%p0(1) + (i - 0.5)*boxarray(nb)%delta
            uc = 0.5*(u(nb)%f(i,j) + u(nb)%f(i+1,j))
            vc = 0.5*(v(nb)%f(i,j) + v(nb)%f(i,j+1))
            vort = 0.25*(u(nb)%f(i,j+1) + u(nb)%f(i+1,j+1) - &
                         u(nb)%f(i,j-1) - u(nb)%f(i+1,j-1)) / boxarray(nb)%delta - &
                   0.25*(v(nb)%f(i+1,j) + v(nb)%f(i+1,j) - &
                         v(nb)%f(i-1,j+1) - v(nb)%f(i-1,j+1)) / boxarray(nb)%delta
            write(10+nb,*) x, y, vort, sqrt(uc**2+vc**2)
          end do
          write(10+nb,*) ''
        end do

      end do

      call destroy_ns_solver()
    end if
 
  end subroutine output

end program poiseuille
