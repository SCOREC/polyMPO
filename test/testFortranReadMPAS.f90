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

  call loadMPASMesh(mpMesh, filename, &
                    maxEdges, vertexDegree, nCells, nVertices)

  !todo check the value using get functions. 
  
  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  stop
end program
