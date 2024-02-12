module testCreateRebuildMPs
    use :: polympo
    use :: readMPAS
    use iso_c_binding
    implicit none

contains
function epsilonDiff(a,b) result(isSame)
    implicit none
    real(kind=MPAS_RKIND) :: a,b,delta
    parameter (delta=1.0e-8)
    logical :: isSame
    if (abs(a-b) < delta) then
      isSame = .true.
    else
      isSame = .false.
    endif
end function

subroutine createMPsTest(mpMesh, nCells, numMPs, mp2Elm, isMPActive, mpPosition)
    implicit none
    type(c_ptr):: mpMesh
    integer :: nCells, numMPs, i
    integer, parameter :: nDims = 3
    real(kind=MPAS_RKIND) :: ptOne = 0.1_MPAS_RKIND
    integer, parameter :: MP_ACTIVE = 1
    integer, parameter :: MP_INACTIVE = 0
    integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive
    real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition
  
    isMPActive = MP_ACTIVE !no inactive MPs and some changed below
    isMPActive(4) = MP_INACTIVE !first/1-st MP is indexed 1 and 4-th MP is inactive
  
    !mp2Elm = [2,3,2,0,3,4,5,6,...]
    mp2Elm(1) = 2
    mp2Elm(2) = 3
    mp2Elm(3) = 2
    !mp2Elm(4) is not needed/used since 4-th MP is inactive
    do i = 5,numMPs
      mp2Elm(i) = i-2 !i=5 leads to mp2Elm(5)=3 (5-th MP in 3-rd element)
                      !i=numMPs leads to mp2Elm(numMPs=nCells+2)=numMPs-2=nCells
    end do
  
    allocate(mpsPerElm(nCells))
    mpsPerElm = 1 !all elements have 1 MP and some changed below
    mpsPerElm(1) = 0 !1st element has 0 MPs
    mpsPerElm(2) = 2 !2nd element has 2 MPs
    mpsPerElm(3) = 2 !3rd element has 2 MPs 

    call polympo_createMPs(mpMesh,nCells,numMPs,c_loc(mpsPerElm),c_loc(mp2Elm),c_loc(isMPActive))

    !set mp positions
    do i = 1,numMPs
      mpPosition(1,i) = i+ptOne
      mpPosition(2,i) = numMPs+i+ptOne
      mpPosition(3,i) = (2*numMPs)+i+ptOne
    end do

    call polympo_setMPPositions(mpMesh,nDims,numMPs,c_loc(mpPosition))
    mpPosition = 42
    call polympo_getMPPositions(mpMesh,nDims,numMPs,c_loc(mpPosition))
    do i = 1,numMPs
      if(isMPActive(i) .eq. MP_ACTIVE) then
        call assert(epsilonDiff(mpPosition(1,i),i+ptOne), "x position of MP does not match")
        call assert(epsilonDiff(mpPosition(2,i),numMPs+i+ptOne), "y position of MP does not match")
        call assert(epsilonDiff(mpPosition(3,i),(2*numMPs)+i+ptOne), "z position of MP does not match")
      endif
    end do
    
    mp2Elm = -99 !override values and then use get function below
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2Elm))
    call assert(mp2Elm(1) .eq. 2, "wrong element ID for MP 1")
    call assert(mp2Elm(2) .eq. 3, "wrong element ID for MP 2")
    call assert(mp2Elm(3) .eq. 2, "wrong element ID for MP 3")
    !mp2Elm(4) is not needed/used since 4-th MP is inactive
    do i = 5,numMPs
      call assert(mp2Elm(i) .eq. i-2, "wrong element ID for i'th MP")
    end do
    !test end

    !deallocate MP variables
    deallocate(mpsPerElm)
end subroutine

subroutine rebuildMPsTests(mpMesh, numMPs, mp2Elm, isMPActive, mpPosition)
    use :: polympo
    use iso_c_binding
    implicit none
    type(c_ptr):: mpMesh
    integer :: numMPs, i, j, numMPsLarger
    integer, dimension(:), pointer :: mp2Elm, addedMPMask, isMPActive, mp2ElmFromPMPO
    integer, dimension(:), pointer :: mp2ElmLarger, addedMPMaskLarger, mp2ElmFromPMPOLarger, isMPActiveLarger
    real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition, mpPositionFromPMPO
    integer, parameter :: nDims = 3
    integer, parameter :: MP_ACTIVE = 1
    integer, parameter :: MP_INACTIVE = 0
    integer, parameter :: MP_DELETE_ELM_ID = -1

    ! TEST: adding MP
    ! Necessary pre-conditions for test
    call assert(numMPs >= 4, "not enough MPs for test")
    call assert(isMPActive(4) == MP_INACTIVE, "mp2Elm = 4 is active")
    ! PREPARE DATA
    allocate(addedMPMask(numMPs))
    addedMPMask = MP_INACTIVE
    isMPActive(4) = MP_ACTIVE
    mp2Elm(4) = 7
    addedMPMask(4) = MP_ACTIVE
    mpPosition(1,4) = 1.2_MPAS_RKIND
    mpPosition(2,4) = 2.2_MPAS_RKIND
    mpPosition(3,4) = 3.2_MPAS_RKIND
    ! Rebuild MPs
    call polympo_startRebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    call polympo_setMPPositions(mpMesh,nDims,numMPs,c_loc(mpPosition))
    call polympo_finishRebuildMPs(mpMesh)
    ! Test values
    allocate(mp2ElmFromPMPO(numMPs))
    mp2ElmFromPMPO = MP_DELETE_ELM_ID
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    allocate(mpPositionFromPMPO(nDims,numMPs))
    call polympo_getMPPositions(mpMesh,nDims,numMPs,c_loc(mpPositionFromPMPO))

    do i = 1, numMPs !check all values match
        if (isMPActive(i) == MP_ACTIVE) then
            call assert(mp2Elm(i) .eq. mp2ElmFromPMPO(i), "wrong element ID for i'th MP after rebuild")
            do j = 1, nDims
                call assert(mpPosition(j,i) == mpPositionFromPMPO(j,i), "mpPosition not set after rebuild")
            end do
        endif
    end do

    ! TEST: deleting two MPs
    ! Necessary pre-conditions for test
    call assert(numMPs >= 4, "not enough MPs for test")
    call assert(isMPActive(1) == MP_ACTIVE, "mp2Elm = 1 not active")
    call assert(isMPActive(4) == MP_ACTIVE, "mp2Elm = 4 not active")
    ! PREPARE DATA
    isMPActive(1) = MP_INACTIVE
    isMPActive(4) = MP_INACTIVE
    mp2Elm(1) = MP_DELETE_ELM_ID
    mp2Elm(4) = MP_DELETE_ELM_ID
    addedMPMask = MP_INACTIVE
    ! Rebuild MPs
    call polympo_startRebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    call polympo_finishRebuildMPs(mpMesh)
    ! Test values
    mp2ElmFromPMPO = MP_DELETE_ELM_ID
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    do i = 1, numMPs !check all values match
        if (isMPActive(i) == MP_ACTIVE) then
            call assert(mp2Elm(i) .eq. mp2ElmFromPMPO(i), "wrong element ID for i'th MP after rebuild")
        endif
    end do

    ! TEST: adding 1, delete 1 and add in same index, removing 1
    ! Necessary pre-conditions for test
    call assert(numMPs >= 3, "not enough MPs for test")
    call assert(isMPActive(1) == MP_INACTIVE, "mp2Elm = 1 not active")
    call assert(isMPActive(2) == MP_ACTIVE, "mp2Elm = 2 is active")
    call assert(isMPActive(3) == MP_ACTIVE, "mp2Elm = 3 not active")
    ! PREPARE DATA
    isMPActive(1) = MP_ACTIVE
    isMPActive(2) = MP_ACTIVE
    isMPActive(3) = MP_INACTIVE
    addedMPMask(1) = MP_ACTIVE
    addedMPMask(2) = MP_ACTIVE
    mp2Elm(1) = 7 !ADDED
    mp2Elm(2) = 7 !REPLACED
    mp2Elm(3) = MP_DELETE_ELM_ID !DELETED
    ! Rebuild MPs
    call polympo_startRebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    call polympo_finishRebuildMPs(mpMesh)
    ! Test values
    mp2ElmFromPMPO = MP_DELETE_ELM_ID
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    do i = 1, numMPs !check all values match
        if (isMPActive(i) == MP_ACTIVE) then
            call assert(mp2Elm(i) .eq. mp2ElmFromPMPO(i), "wrong element ID for i'th MP after rebuild")
        endif
    end do

    ! TEST: increasing numMPs to numMPsLarger
    ! TEST: adding 1 MP (within numMPs), removing 1 MP (within numMPs), adding 1 (after numMPs, but within numMPsLarger)
    ! Necessary pre-conditions for test
    call assert(numMPs >= 5, "not enough MPs for test")
    call assert(isMPActive(4) == MP_INACTIVE, "mp2Elm = 4 not active")
    call assert(isMPActive(5) == MP_ACTIVE, "mp2Elm = 5 is active")
    ! PREPARE DATA
    numMPsLarger = numMPs+10
    allocate(mp2ElmLarger(numMPsLarger))
    allocate(isMPActiveLarger(numMPsLarger))
    allocate(addedMPMaskLarger(numMPsLarger))
    
    mp2ElmLarger = MP_DELETE_ELM_ID
    isMPActiveLarger = MP_INACTIVE
    addedMPMaskLarger = MP_INACTIVE
    do i = 1, numMPs
        mp2ElmLarger(i) = mp2Elm(i)
        isMPActiveLarger(i) = isMPActive(i)
    end do

    isMPActiveLarger(4) = MP_ACTIVE ! within numMPs
    isMPActiveLarger(5) = MP_INACTIVE ! within numMPs
    isMPActiveLarger(numMPsLarger-2) = MP_ACTIVE ! within numMPsLarger
    mp2ElmLarger(4) = 7
    mp2ElmLarger(5) = MP_DELETE_ELM_ID
    mp2ElmLarger(numMPsLarger-2) =  7
    addedMPMaskLarger(4) = MP_ACTIVE
    addedMPMaskLarger(numMPsLarger-2) = MP_ACTIVE
    ! Rebuild MPs
    call polympo_startRebuildMPs(mpMesh,numMPsLarger,c_loc(mp2ElmLarger),c_loc(addedMPMaskLarger))
    call polympo_finishRebuildMPs(mpMesh)
    ! Test values
    allocate(mp2ElmFromPMPOLarger(numMPsLarger))
    mp2ElmFromPMPOLarger = MP_DELETE_ELM_ID
    call polympo_getMPCurElmID(mpMesh,numMPsLarger,c_loc(mp2ElmFromPMPOLarger))
    do i = 1, numMPs !check all values match
        if (isMPActiveLarger(i) == MP_ACTIVE) then
            call assert(mp2ElmLarger(i) .eq. mp2ElmFromPMPOLarger(i), "wrong element ID for i'th MP after rebuild")
        endif
    end do
    ! Cleanup
    deallocate(addedMPMask)
    deallocate(mp2ElmFromPMPO)
    deallocate(mp2ElmLarger)
    deallocate(isMPActiveLarger)
    deallocate(addedMPMaskLarger)
    deallocate(mp2ElmFromPMPOLarger)
end subroutine
end module testCreateRebuildMPs

program main
  use :: polympo
  use :: readMPAS
  use :: testCreateRebuildMPs
  use :: iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: ierr, self
  integer :: argc, i, arglen
  integer :: setMeshOption, setMPOption
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  character (len=2048) :: filename
  type(c_ptr) :: mpMesh
  integer :: numMPs
  integer, dimension(:), pointer :: mp2Elm, isMPActive
  integer, parameter :: MP_ACTIVE = 1
  integer, parameter :: nDims = 3
  integer, parameter :: MP_INACTIVE = 0
  character (len=64) :: onSphere, stringYes = "YES"
  real(kind=MPAS_RKIND) :: sphereRadius
  integer :: maxEdges, vertexDegree, nCells, nVertices
  integer, dimension(:), pointer :: nEdgesOnCell
  real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  integer, dimension(:,:), pointer :: cellsOnVertex, verticesOnCell, cellsOnCell
  real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle)
  call polympo_initialize()

  argc = command_argument_count()
  if(argc == 1) then
    call get_command_argument(1, filename)
  else
    write(0, *) "Usage: ./testFortranCreateMPs <path to the nc file>"
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

  !check for allocation
  call assert(nCells .ge. 3, "This test requires a mesh with at least three cells")
  numMPs = nCells+2;
  allocate(mp2Elm(numMPs))
  allocate(isMPActive(numMPs))
  allocate(mpPosition(nDims,numMPs))
  
  call createMPsTest(mpMesh, nCells, numMPs, mp2Elm, isMPActive, mpPosition) 
  call rebuildMPsTests(mpMesh, numMPs, mp2Elm, isMPActive, mpPosition)
  
  deallocate(mp2Elm)
  deallocate(isMPActive)
  deallocate(mpPosition)

  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)
 
  !deallocate other mesh variables 
  deallocate(nEdgesOnCell)
  deallocate(xVertex)
  deallocate(yVertex)
  deallocate(zVertex)
  deallocate(verticesOnCell)
  deallocate(cellsOnCell)

  stop
end program
