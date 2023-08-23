subroutine assert(condition,message)
  implicit none
  logical :: condition
  character(*) :: message
  if (condition .neqv. .true.) then
    write (*,*) message
    call exit(1)
  endif
end subroutine

program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'
    
  integer, parameter :: APP_RKIND = selected_real_kind(15)
  integer :: nverts, numComps, coordDegree, numMPs, numElms
  integer :: i, j
  integer :: setMeshOption, setMPOption
  integer :: mpi_comm_handle = MPI_COMM_WORLD
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
  setMeshOption = 1
  setMPOption = 1    
  mpMesh = polympo_createMPMesh(setMeshOption,setMPOption) !creates test mesh
  
  nverts = 19 !todo use getNumVtx from the Mesh object
  numComps = 2 !todo use getNumComps from velocity fields
  coordDegree = 3
  numMPs = 49 !todo use getNumMPs from the MaterialPoints object
  numElms = 10

  allocate(Mesharray(numComps,nverts))
  allocate(MParray(numComps,numMPs))
  allocate(MPElmID(numMPs))
  allocate(MPPositions(coordDegree,numMPs))
  allocate(xArray(nverts))
  allocate(yArray(nverts))
  allocate(zArray(nverts))

  call polympo_getMPPositions(mpMesh, coordDegree, numMPs, c_loc(MPPositions))
  do i = 1,numMPs 
    call assert(abs(MPPositions(3,i) - 1.1) .lt. 0.000001, "Assert zPositions for MP array Fail")
  end do

  do i = 1,numMPs 
    MPElmID(i) = mod(i, numElms)
  end do
  call polympo_setMPCurElmID(mpMesh, numMPs, c_loc(MPElmID))
  call polympo_getMPCurElmID(mpMesh, numMPs, c_loc(MPElmID))
  do i = 1,numMPs 
    call assert((MPElmID(i) .eq. mod(i, numElms)) , "Assert MPElmID Fail")
  end do

  value1 = 42
  MParray = value1
  call polympo_setMPVelArray(mpMesh, numMPs, c_loc(MParray))
  
  MParray = 1
  call polympo_getMPVelArray(mpMesh, numMPs, c_loc(MParray))
  call assert(all(MParray .eq. value1), "Assert MParray == value1 Failed!")

  value2 = 24
  MParray = value2
  call polympo_setMPVelArray(mpMesh, numMPs, c_loc(MParray))

  MParray = 1
  call polympo_getMPVelArray(mpMesh, numMPs, c_loc(MParray))
  call assert(all(MParray .eq. value2), "Assert MParray == value2 Failed!")

  call polympo_startMeshFill(mpMesh)
  do i = 1,numComps
    do j = 1,nverts 
        Mesharray(i,j) = (i-1)*numComps + j
    end do
  end do
  call polympo_setMeshOnSurfVeloIncrArray(mpMesh, numComps, nverts, c_loc(Mesharray))
  call polympo_setMeshOnSurfDispIncrArray(mpMesh, numComps, nverts, c_loc(Mesharray))

  Mesharray = 1
  call polympo_getMeshOnSurfVeloIncrArray(mpMesh, numComps, nverts, c_loc(Mesharray))
  do i = 1,numComps
    do j = 1,nverts 
        call assert((Mesharray(i,j) .eq. (i-1)*numComps+j), "Assert 2d array Fail")
    end do
  end do
  Mesharray = 1
  call polympo_getMeshOnSurfDispIncrArray(mpMesh, numComps, nverts, c_loc(Mesharray))
  do i = 1,numComps
    do j = 1,nverts 
        call assert((Mesharray(i,j) .eq. (i-1)*numComps+j), "Assert 2d array Fail")
    end do
  end do

  value1 = 1337
  value2 = 42
  xArray = value1
  yArray = value2
  zArray = value1 + value2 
  call polympo_setMeshVtxCoords(mpMesh, nverts, c_loc(xArray), c_loc(yArray), c_loc(zArray))
  xArray = 1
  yArray = 1
  zArray = 1 
  call polympo_getMeshVtxCoords(mpMesh, nverts, c_loc(xArray), c_loc(yArray), c_loc(zArray))
  call assert(all(xArray .eq. value1), "Assert xArray == value1 Failed!")
  call assert(all(yArray .eq. value2), "Assert yArray == value2 Failed!")
  call assert(all(zArray .eq. value1 + value2), "Assert zArray == value1 + value2 Failed!")

  call polympo_endMeshFill(mpMesh)

  deallocate(MParray)
  deallocate(Mesharray)
  deallocate(xArray)
  deallocate(yArray)
  deallocate(zArray)

  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  stop
end
