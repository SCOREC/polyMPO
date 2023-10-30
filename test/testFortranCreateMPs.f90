program main
  use :: polympo
  use :: readMPAS
  use :: iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: ierr, self
  integer :: argc, i, arglen
  integer :: setMeshOption, setMPOption
  integer, parameter :: nDims = 3
  integer, parameter :: MP_ACTIVE = 1
  integer, parameter :: MP_INACTIVE = 0
  real(kind=MPAS_RKIND) :: ptOne = 0.100000000000000000
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  character (len=2048) :: filename
  type(c_ptr) :: mpMesh
  integer :: numMPs
  integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive
  real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition
  character (len=64) :: onSphere, stringYes = "YES"
  real(kind=MPAS_RKIND) :: sphereRadius
  integer :: maxEdges, vertexDegree, nCells, nVertices
  integer, dimension(:), pointer :: nEdgesOnCell
  real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell

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
                        verticesOnCell, cellsOnCell)
  call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
                        verticesOnCell, cellsOnCell)

    !test on new createMPs
    call assert(nCells .ge. 3, "This test requires a mesh with at least three cells")
    numMPs = nCells+2;
    allocate(mpsPerElm(nCells))
    allocate(mp2Elm(numMPs))
    allocate(isMPActive(numMPs))
    
    isMPActive = MP_ACTIVE !no inactive MPs and some changed below
    isMPActive(4) = MP_INACTIVE !first/1-st MP is indexed 1 and 4-th MP is inactive
   
    mpsPerElm = 1 !all elements have 1 MP and some changed below
    mpsPerElm(1) = 0 !1st element has 0 MPs
    mpsPerElm(2) = 2 !2nd element has 2 MPs
    mpsPerElm(3) = 2 !3rd element has 2 MPs 

    ! mp2Elm = [2,3,2,0,3,4,5,6,...]
    mp2Elm(1) = 2
    mp2Elm(2) = 3
    mp2Elm(3) = 2
    !mp2Elm(4) is not needed/used since 4-th MP is inactive
    do i = 5,numMPs
      mp2Elm(i) = i-2 !i=5 leads to mp2Elm(5)=3 (5-th MP in 3-rd element)
                      !i=numMPs leads to mp2Elm(numMPs=nCells+2)=numMPs-2=nCells
    end do
    call polympo_createMPs(mpMesh,nCells,numMPs,c_loc(mpsPerElm),c_loc(mp2Elm),c_loc(isMPActive))

    !set mp positions
    allocate(mpPosition(nDims,numMPs))
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
    deallocate(mpPosition)
    deallocate(mpsPerElm)
    deallocate(mp2Elm)
    deallocate(isMPActive)

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
