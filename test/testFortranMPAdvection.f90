module mpMesh_ptr
use polympo
use iso_c_binding
implicit none

  type(c_ptr):: mpMesh
contains
  subroutine createMPMesh(setMeshOption, setMPOption)
    integer:: setMeshOption, setMPOption
    mpMesh = polympo_createMPMesh(setMeshOption,setMPOption) !creates test mesh
  end subroutine
end module
!---------------------------------------------------------------------------
!> todo add a discription
!---------------------------------------------------------------------------
program main
  use :: polympo
  use :: readMPAS
  use :: iso_c_binding
  use mpMesh_ptr
  implicit none
  include 'mpif.h'

  integer, parameter :: APP_RKIND = selected_real_kind(15)
  integer :: ierr, self
  integer :: argc, i, arglen
  integer :: setMeshOption, setMPOption
  integer :: maxEdges, vertexDegree, nCells, nVertices
  integer :: nCompsDisp
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  character (len=2048) :: filename
  real(kind=APP_RKIND), dimension(:,:), pointer :: dispIncr
  character (len=64) :: onSphere, stringYes = "YES"
  real(kind=MPAS_RKIND) :: sphereRadius
  integer, dimension(:), pointer :: nEdgesOnCell
  real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
  integer :: numMPs 
  integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle)
  call polympo_initialize()

  call polympo_checkPrecisionForRealKind(APP_RKIND)
  argc = command_argument_count()
  if(argc == 1) then
    call get_command_argument(1, filename)
  else
    write(0, *) "Usage: ./testFortranInterpolatePush <path to the nc file>"
  end if

  setMeshOption = 0 !create an empty mesh
  setMPOption = 0   !create an empty set of MPs
  nCompsDisp = 2
  call createMPMesh(setMeshOption, setMPOption)
  call readMPASMeshFromNCFile(filename, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
                        verticesOnCell, cellsOnCell)
  call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
                        verticesOnCell, cellsOnCell)

  allocate(dispIncr(nCompsDisp,nVertices))
  !createMPs
  numMPs = nCells
  allocate(mpsPerElm(nCells))
  allocate(mp2Elm(numMPs))
  allocate(isMPActive(numMPs))
  isMPActive = 1 !no inactive MPs and some changed below
  mpsPerElm = 1 !all elements have 1 MP and some changed below
  do i = 1,numMPs
    mp2Elm(i) = i
  end do
  call polympo_createMPs(mpMesh,nCells,numMPs,c_loc(mpsPerElm),c_loc(mp2Elm),c_loc(isMPActive))

  ! set disp incr
  dispIncr = 1
  !do i = 1,numMPs
  !  dispIncr(0,i) = 0
  !  dispIncr(1,i) = 0 !average delta lambda over mesh edges
  !end do
  call polympo_setMeshOnSurfDispIncr(mpMesh, nCompsDisp, nVertices, c_loc(dispIncr))
  call polympo_push(mpMesh)
 
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
