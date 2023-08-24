!---------------------------------------------------------------------------
!> This is a test/example on how the polyMPO work
!---------------------------------------------------------------------------
program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: ierr, self
  integer :: mpi_comm_handle = MPI_COMM_WORLD

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  !polympo start here
  call polympo_setMPICommunicator(mpi_comm_handle)
  call polympo_initialize()

  call polympo_finalize()
  !polympo end
    
  call mpi_finalize(ierr)

  stop
end
