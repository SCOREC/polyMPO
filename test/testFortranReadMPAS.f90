program main
  use polympo
  use iso_c_binding
  use mpi_f08
  implicit none

  integer :: ierr, self
  TYPE(MPI_Comm) :: mpi_comm_handle = MPI_COMM_WORLD
  character (len=128) :: filename = "./grid_full.nc"
  integer :: nCells, nVertices
  integer(c_int), dimension(:), pointer :: nEdgesOnCell
  real(c_double), dimension(:), pointer :: xVertex, yVertex, zVertex
  integer(c_int), dimension(:,:), allocatable, target :: verticesOnCell, cellsOnVertex, cellsOnCell

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)
  call polympo_setCommunicator(mpi_comm_handle%MPI_VAL)
  call polympo_initialize()

  call polympo_readMPASMesh(filename, &
                            nCells, nVertices, nEdgesOnCell, &
                            xVertex, yVertex, zVertex, &
                            verticesOnCell, cellsOnVertex, cellsOnCell)

  deallocate(nEdgesOnCell)
  deallocate(xVertex)
  deallocate(yVertex)
  deallocate(zVertex)
  deallocate(verticesOnCell)
  deallocate(cellsOnVertex)
  deallocate(cellsOnCell)
  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
