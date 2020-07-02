program check_grid

  ! Check consistency one generating 2D grids

  use grid_2D

  implicit none
  
  integer :: n, i, j
  real(kind=dp) :: a, b
  character(len=1) :: sn
  character(len=4) :: filename

  ! First generate one single grid
  n = nbox
  allocate(boxarray(nbox))
  call create_box(1, 20, 20, 1.d0, 1.d0, 0.d0, 0.d0)

  print *, 'Grid ', n, 'of ', nbox
  print *, 'Size: ', boxarray(n)%Lx, boxarray(n)%Ly, boxarray(n)%nx, boxarray(n)%ny, boxarray(n)%delta
  print *, 'Origin :', boxarray(n)%p0
  print *, 'x an y array szie: ', size(boxarray(n)%x), size(boxarray(n)%y)
  print *, 'Physical BC: '
  print *, '             left   ', boxarray(n)%left_boundary 
  print *, '             right  ', boxarray(n)%right_boundary 
  print *, '             top    ', boxarray(n)%top_boundary 
  print *, '             bottom ', boxarray(n)%bottom_boundary 
  print *, 'Interanl BC: '
  print *, '             left   ', boxarray(n)%left 
  print *, '             right  ', boxarray(n)%right 
  print *, '             top    ', boxarray(n)%top 
  print *, '             bottom ', boxarray(n)%bottom 
  print *, ''

  ! Try to access all position in the grid
  do j = 1,boxarray(n)%ny
    do i = 1,boxarray(n)%nx
      a = boxarray(n)%x(i,j)
      b = boxarray(n)%y(i,j)
    end do
  end do

  call destroy_boxes()

  ! Now do the same for three grids 
  nbox = 3
  allocate(boxarray(nbox))
  do n = 1,nbox
    call create_box(n, 20, 20, 1.d0, 1.d0, 1.d0*(n-1), 1.d0*(n-1))

    print *, 'Grid ', n, 'of ', nbox
    print *, 'Sizes: ', boxarray(n)%Lx, boxarray(n)%Ly, boxarray(n)%nx, boxarray(n)%ny, boxarray(n)%delta
    print *, 'Origin: ', boxarray(n)%p0
    print *, ' x and y array size: ', size(boxarray(n)%x), size(boxarray(n)%y)
    print *, 'Physical BC: '
    print *, '             left   ', boxarray(n)%left_boundary 
    print *, '             right  ', boxarray(n)%right_boundary 
    print *, '             top    ', boxarray(n)%top_boundary 
    print *, '             bottom ', boxarray(n)%bottom_boundary 
    print *, 'Interanl BC: '
    print *, '             left   ', boxarray(n)%left 
    print *, '             right  ', boxarray(n)%right 
    print *, '             top    ', boxarray(n)%top 
    print *, '             bottom ', boxarray(n)%bottom 
    print *, ''

    ! Try to access all position in the grid
    write(sn,'(I1)') n
    filename = 'box'//sn
    open(unit = n, file = filename) 
    do j = 1,boxarray(n)%ny
      do i = 1,boxarray(n)%nx
        a = boxarray(n)%x(i,j)
        b = boxarray(n)%y(i,j)
        write(n,*) boxarray(n)%x(i,j), boxarray(n)%y(i,j)
      end do
    end do
    flush(n)
    close(n)

  end do

  call destroy_boxes()

end program

