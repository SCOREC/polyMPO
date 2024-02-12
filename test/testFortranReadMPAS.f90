!---------------------------------------------------------------------------
!> This is a test on how to use loadMPASMesh
!> For specific usage on setting mesh properties, see test/readMPAS.f90
!---------------------------------------------------------------------------
program main
  use :: polympo
  use :: readMPAS
  use :: iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: ierr, self
  integer :: argc, i, arglen
  integer :: setMeshOption, setMPOption
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  integer :: maxEdges, vertexDegree, nCells, nVertices
  character (len=2048) :: filename
  type(c_ptr) :: mpMesh
  character (len=64) :: onSphere, stringYes = "YES"
  integer :: numMPs
  real(kind=MPAS_RKIND) :: sphereRadius
  integer, dimension(:), pointer :: nEdgesOnCell
  real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  integer, dimension(:,:), pointer :: cellsOnVertex, verticesOnCell, cellsOnCell

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle)
  call polympo_initialize()

  argc = command_argument_count()
  if(argc == 1) then
    call get_command_argument(1, filename)
  else
    write(0, *) "Usage: ./testFortranReadMPAS <path to the nc file>"
  end if

  setMeshOption = 0 !create an empty mesh
  setMPOption = 0   !create an empty set of MPs
  mpMesh = polympo_createMPMesh(setMeshOption, setMPOption)

  call readMPASMeshFromNCFile(filename, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
                        verticesOnCell, cellsOnCell, cellsOnVertex)
  call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, &
                        verticesOnCell, cellsOnCell, cellsOnVertex)

  !todo check the value using get functions. 
  
  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  deallocate(nEdgesOnCell)
  deallocate(xVertex)
  deallocate(yVertex)
  deallocate(zVertex)
  deallocate(verticesOnCell)
  deallocate(cellsOnCell)
  
  stop
end program
