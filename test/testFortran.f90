program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: nverts = 3
  real(c_double), dimension(:), allocatable :: array
  integer :: ierr, self
  type(c_ptr) :: mpMesh
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)
  call polympo_initialize()
  call polympo_createMpMesh(mpMesh)
  allocate(array(nverts))
  array = 1
  write (*,*) array
  call polympo_modifyArray(nverts,array)
  write (*,*) array
  deallocate(array)
  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
