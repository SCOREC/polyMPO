program main
  use polympo
  use iso_c_binding
  use mpi_f08
  implicit none

  integer :: ierr, self
  TYPE(MPI_Comm) :: mpi_comm_handle = MPI_COMM_WORLD
  
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)
  call polympo_setCommunicator(mpi_comm_handle)
  !call polympo_setCommunicator(mpi_comm_handle%MPI_VAL)
  call polympo_initialize()

  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
