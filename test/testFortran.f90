program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: nverts
  real(c_double), dimension(:), allocatable :: array
  integer :: ierr, self
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)
  call polympo_initialize()
  allocate(array(nverts))
  call polympo_modifyArray(nverts,array)
  deallocate(array)
  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
