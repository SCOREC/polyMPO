subroutine assert(condition,message)
  implicit none
  logical :: condition
  character(*) :: message
  if (condition .neqv. .true.) then
    write (*,*) message
    call exit(1)
  endif
end subroutine

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
!> This is a test of use a module to create mpMesh pointer for polyMPO
!> The test is a demonstration of how to use a module to create mpMesh object
!> For specific usage of binding functions, see src/pmpo_fortran 
!---------------------------------------------------------------------------
program main
  use polympo
  use iso_c_binding
  use mpMesh_ptr
  implicit none
  include 'mpif.h'
    
  integer, parameter :: APP_RKIND = selected_real_kind(15)
  integer :: nverts, numCompsVel, numCompsCoords, numMPs, numElms
  integer :: i, j
  integer :: setMeshOption, setMPOption
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  real(kind=APP_RKIND) :: test_epsilon = 1e-6
  real(kind=APP_RKIND) :: value1, value2
  integer, dimension(:), pointer :: MPElmID
  real(kind=APP_RKIND), dimension(:,:), pointer :: MPPositions
  real(kind=APP_RKIND), dimension(:,:), pointer :: Mesharray
  integer :: ierr, self

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle) !this is not supported yet! only for showing
  call polympo_initialize()

  call polympo_checkPrecisionForRealKind(APP_RKIND)
  setMeshOption = 1 !create a hard coded planar test mesh
  setMPOption = 1 !create some random test MPs that based on the mesh option you give
  !initialize the mpMesh object using a module
  call createMPMesh(setMeshOption, setMPOption) 

  !These are hard coded test mesh values 
  nverts = 19 !todo use getNumVtx from the Mesh object
  numCompsVel = 2 !todo use getNumComps from velocity fields
  numCompsCoords = 3
  numMPs = 49 !todo use getNumMPs from the MaterialPoints object
  numElms = 10
  value1 = 1337
  value2 = 42

  allocate(Mesharray(numCompsVel,nverts))
  allocate(MPElmID(numMPs))
  allocate(MPPositions(numCompsCoords,numMPs))

  MPPositions = 0
  call polympo_getMPPositions(mpMesh, numCompsCoords, numMPs, c_loc(MPPositions))
  call assert(all(MPPositions .ne. 0), "Assert MPPositions failed!")

  Mesharray = value1
  call polympo_setMeshOnSurfVeloIncr(mpMesh, numCompsVel, nverts, c_loc(Mesharray))
  Mesharray = value2
  call polympo_setMeshOnSurfDispIncr(mpMesh, numCompsVel, nverts, c_loc(Mesharray))

  Mesharray = 1
  call polympo_getMeshOnSurfVeloIncr(mpMesh, numCompsVel, nverts, c_loc(Mesharray))
  call assert(all(Mesharray .eq.value1), "Assert OnSurfVeloIncr failed!")
  Mesharray = 1
  call polympo_getMeshOnSurfDispIncr(mpMesh, numCompsVel, nverts, c_loc(Mesharray))
  call assert(all(Mesharray .eq.value2), "Assert OnSurfDispIncr failed!")


  deallocate(Mesharray)

  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  stop
end
