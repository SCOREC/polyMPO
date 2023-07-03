program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: ierr, self
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)
  call polympo_setCommunicator(MPI_COMM_WORLD)
  call polympo_initialize()

  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
