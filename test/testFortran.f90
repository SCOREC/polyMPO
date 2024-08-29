subroutine assert(condition,message)
  implicit none
  logical :: condition
  character(*) :: message
  if (condition .neqv. .true.) then
    write (*,*) message
    call exit(1)
  endif
end subroutine

!---------------------------------------------------------------------------
!> This is a test on the set/get APIs provided in polyMPO
!> These fields can be set and get at anytime
!> For specific usage, see src/pmpo_fortran 
!---------------------------------------------------------------------------
program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'
    
  integer, parameter :: APP_RKIND = selected_real_kind(15)
  integer :: nverts, nvertsGet, nElmsGet, numCompsVel, numCompsCoords, numMPs, numElms
  integer :: i, j
  integer :: setMeshOption, setMPOption
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  real(kind=APP_RKIND) :: test_epsilon = 1e-6
  real(kind=APP_RKIND) :: value1, value2
  integer, dimension(:), pointer :: MPElmID
  real(kind=APP_RKIND), dimension(:,:), pointer :: MParray
  real(kind=APP_RKIND), dimension(:,:), pointer :: MPPositions
  real(kind=APP_RKIND), dimension(:,:), pointer :: Mesharray
  real(kind=APP_RKIND), dimension(:), pointer :: xArray, yArray, zArray
  integer :: ierr, self
  type(c_ptr) :: mpMesh

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle) !this is not supported yet! only for showing
  call polympo_initialize()

  call polympo_checkPrecisionForRealKind(APP_RKIND)
  setMeshOption = 1 !create a hard coded planar test mesh
  setMPOption = 1 !create some random test MPs that based on the mesh option you give
  mpMesh = polympo_createMPMesh(setMeshOption,setMPOption) !creates test mesh
 
  !These are hard coded test mesh values 
  nverts = 19 !todo use getNumVtx from the Mesh object
  numCompsVel = 2 !todo use getNumComps from velocity fields
  numCompsCoords = 3
  numMPs = 49 !todo use getNumMPs from the MaterialPoints object
  numElms = 10

  allocate(Mesharray(numCompsVel,nverts))
  allocate(MParray(numCompsVel,numMPs))
  allocate(MPElmID(numMPs))
  allocate(MPPositions(numCompsCoords,numMPs))
  allocate(xArray(nverts))
  allocate(yArray(nverts))
  allocate(zArray(nverts))

  ! MPMesh structure
  call polympo_getMPPositions(mpMesh, numCompsCoords, numMPs, c_loc(MPPositions))
  do i = 1,numMPs 
    call assert(abs(MPPositions(3,i) - 1.1) .lt. test_epsilon, "Assert zPositions for MP array Fail")
  end do
  call polympo_getMeshNumVtxs(mpMesh, nvertsGet)
  call assert(nverts.eq.nvertsGet,"num. verts mismatch")
  call polympo_getMeshNumElms(mpMesh, nElmsGet)
  call assert(numElms.eq.nElmsGet,"num. elms mismatch")

  ! set MP Fields
  do i = 1,numCompsVel
    do j = 1,numMPs 
        MParray(i,j) = (i-1)*numMPs + j
    end do
  end do
  call polympo_setMPVel(mpMesh, numCompsVel, numMPs, c_loc(MParray))
  
  ! check MP Fields
  MParray = -1
  call polympo_getMPVel(mpMesh, numCompsVel, numMPs, c_loc(MParray))
  do i = 1,numCompsVel
    do j = 1,numMPs 
        call assert((MParray(i,j) .eq. (i-1)*numMPs+j), "Assert MPVel Fail")
    end do
  end do

  ! set mesh Fields
  do i = 1,numCompsVel
    do j = 1,nverts 
        Mesharray(i,j) = (i-1)*nverts + j
    end do
  end do
  call polympo_setMeshVtxOnSurfVeloIncr(mpMesh, numCompsVel, nverts, c_loc(Mesharray))
  call polympo_setMeshVtxOnSurfDispIncr(mpMesh, numCompsVel, nverts, c_loc(Mesharray))

  ! check mesh Fields vtxBased
  Mesharray = -1
  call polympo_getMeshVtxOnSurfVeloIncr(mpMesh, numCompsVel, nverts, c_loc(Mesharray))
  do i = 1,numCompsVel
    do j = 1,nverts 
        call assert((Mesharray(i,j) .eq. (i-1)*nverts+j), "Assert MeshOnSurfVeloIncr Fail")
    end do
  end do
  Mesharray = -1
  call polympo_getMeshVtxOnSurfDispIncr(mpMesh, numCompsVel, nverts, c_loc(Mesharray))
  do i = 1,numCompsVel
    do j = 1,nverts 
        call assert((Mesharray(i,j) .eq. (i-1)*nverts+j), "Assert MeshOnSurfDispIncr Fail")
    end do
  end do

  value1 = 1337
  value2 = 42
  xArray = value1
  yArray = value2
  zArray = value1 + value2 
  call polympo_setMeshVtxCoords(mpMesh, nverts, c_loc(xArray), c_loc(yArray), c_loc(zArray))
  xArray = -1
  yArray = -1
  zArray = -1 
  call polympo_getMeshVtxCoords(mpMesh, nverts, c_loc(xArray), c_loc(yArray), c_loc(zArray))
  call assert(all(xArray .eq. value1), "Assert xArray == value1 Failed!")
  call assert(all(yArray .eq. value2), "Assert yArray == value2 Failed!")
  call assert(all(zArray .eq. value1 + value2), "Assert zArray == value1 + value2 Failed!")

  !use xArray and yArrray to hold u and v components 
  do i = 1, nverts
    xArray(i) = i + value1
    yArray(i) = value2 - i
  end do 
  call polympo_setMeshVtxVel(mpMesh, nverts, c_loc(xArray),c_loc(yArray))
  xArray = -1
  yArray = -1
  call polympo_getMeshVtxVel(mpMesh, nverts, c_loc(xArray),c_loc(yArray))
  do i = 1, nverts
    call assert((xArray(i) .eq. i+value1), "Assert MeshVel u-component Velocity Fail")
    call assert((yArray(i) .eq. value2-i), "Assert MeshVel v-component Velocity Fail")
  end do
  
  call polympo_setMeshVtxMass(mpMesh, nverts, c_loc(xArray))
  xArray = -1 
  call polympo_getMeshVtxMass(mpMesh, nverts, c_loc(xArray))
  do i = 1, nverts
    call assert((xArray(i) .eq. i+value1), "Assert MeshVtxMass Mass Fail")
  end do

  ! check mesh Fields vtxBased
  deallocate(xArray)
  allocate(xArray(numElms))
  
  do i = 1, numElms
    xArray(i) = i + value1
  end do 
  
  call polympo_setMeshElmMass(mpMesh, numElms, c_loc(xArray))
  xArray = -1 
  call polympo_getMeshElmMass(mpMesh, numElms, c_loc(xArray))
  do i = 1, numElms
    call assert((xArray(i) .eq. i+value1), "Assert MeshVtxMass Mass Fail")
  end do

  deallocate(MParray)
  deallocate(Mesharray)
  deallocate(xArray)
  deallocate(yArray)
  deallocate(zArray)

  ! test elmcenter
  allocate(xArray(numElms))
  allocate(yArray(numElms))
  allocate(zArray(numElms))
  xArray = value1
  yArray = value2
  zArray = value1 + value2 
  call polympo_setMeshElmCenter(mpMesh, numElms, c_loc(xArray), c_loc(yArray), c_loc(zArray))
  xArray = -1
  yArray = -1
  zArray = -1
  call polympo_getMeshElmCenter(mpMesh, numElms, c_loc(xArray), c_loc(yArray), c_loc(zArray))
  call assert(all(xArray .eq. value1), "Assert xArray == value1 Failed!")
  call assert(all(yArray .eq. value2), "Assert yArray == value2 Failed!")
  call assert(all(zArray .eq. value1 + value2), "Assert zArray == value1 + value2 Failed!")
  
  deallocate(xArray)
  deallocate(yArray)
  deallocate(zArray)
  
  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  stop
end
