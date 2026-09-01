program main
  implicit none

  type :: my_type
     integer, pointer :: my_size(:)      ! F95
     !integer, allocatable :: my_size(:) ! F95 + TR 15581 or F2003
  end type my_type

  type(my_type), allocatable :: x(:)

  allocate(x(3))

  allocate(x(1)%my_size(3))
  allocate(x(2)%my_size(2))
  allocate(x(3)%my_size(1))

  print*, x(1)%my_size
  print*, x(2)%my_size
  print*, x(3)%my_size

  deallocate(x(3)%my_size, x(2)%my_size, x(1)%my_size)
  deallocate(x)


subroutine hello()
end subroutine hello
end program main
