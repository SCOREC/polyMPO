module advectionTests
  contains
  include "calculateDisplacement.f90"

  subroutine runAdvectionTest(mpMesh, numPush, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)
    use :: polympo
    use :: readMPAS
    use :: iso_c_binding
    implicit none

    type(c_ptr) :: mpMesh
    integer :: i, numPush, nVertices
    real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
    integer, dimension(:), pointer :: nEdgesOnCell
    integer, dimension(:,:), pointer :: verticesOnCell
    real(kind=MPAS_RKIND) :: sphereRadius


    do i = 1, numPush
      call calcSurfDispIncr(mpMesh, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)
      call polympo_push(mpMesh)
    end do

  end subroutine

  subroutine runAdvectionTest2(mpMesh, numPush, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)
    use :: polympo
    use :: readMPAS
    use :: iso_c_binding
    implicit none

    type(c_ptr) :: mpMesh
    integer :: i, numPush, nVertices
    real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
    integer, dimension(:), pointer :: nEdgesOnCell
    integer, dimension(:,:), pointer :: verticesOnCell
    real(kind=MPAS_RKIND) :: sphereRadius

    PRINT *, "Foward: "
    do i = 1, numPush
      call calcSurfDispIncr(mpMesh, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)
      call polympo_push(mpMesh)
    end do

    PRINT *, "Backward: "
    call calcSurfDispIncr(mpMesh, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius, -numPush)
    call polympo_push(mpMesh)
   
  end subroutine



  subroutine runReconstructionTest(mpMesh, numMPs, numPush, nCells, nVertices, mp2Elm, &
                                  latVertex, lonVertex, nEdgesOnCell, verticesOnCell, sphereRadius)
    use :: polympo
    use :: readMPAS
    use :: iso_c_binding
    implicit none

    type(c_ptr) :: mpMesh
    integer :: i, j, k, vID, numMPs, numPush, nVertices, nCells
    real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpMass, mpVel
    real(kind=MPAS_RKIND), dimension(:), pointer :: meshVtxMass, meshElmMass
    real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
    integer, dimension(:), pointer :: nEdgesOnCell
    integer, dimension(:,:), pointer :: verticesOnCell
    integer, dimension(:), pointer :: mp2Elm
    real(kind=MPAS_RKIND) :: sphereRadius
    real(kind=MPAS_RKIND) :: TEST_VAL = 1.1_MPAS_RKIND
    real(kind=MPAS_RKIND) :: TOLERANCE = 0.0001_MPAS_RKIND

    allocate(mpMass(1,numMPs))
    allocate(mpVel(2,numMPs))
    allocate(meshVtxMass(nVertices))
    allocate(meshElmMass(nCells))

    mpMass = TEST_VAL
    mpVel = TEST_VAL

    call polympo_setMPMass(mpMesh,1,numMPs,c_loc(mpMass))
    call polympo_setMPVel(mpMesh,2,numMPs,c_loc(mpVel))

    ! ! Test vtx reconstruction
    
    ! call polympo_setReconstructionOfMass(mpMesh,0,polympo_getMeshFVtxType())
    ! call polympo_applyReconstruction(mpMesh)
    ! call polympo_getMeshVtxMass(mpMesh,nVertices,c_loc(meshVtxMass))

    ! do i = 1, nVertices
    !   call assert(meshVtxMass(i) < TEST_VAL+TOLERANCE .and. meshVtxMass(i) > TEST_VAL-TOLERANCE, "Error: wrong vtx mass")
    ! end do
    
    ! Test push reconstruction

    do j = 1, numPush
      call calcSurfDispIncr(mpMesh, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)
      call polympo_setReconstructionOfMass(mpMesh,0,polympo_getMeshFElmType())
      call polympo_setReconstructionOfMass(mpMesh,0,polympo_getMeshFVtxType())
      call polympo_setReconstructionOfVel(mpMesh,0,polympo_getMeshFVtxType())
      call polympo_push(mpMesh)
      ! call polympo_getMeshElmMass(mpMesh,nCells,c_loc(meshElmMass))
      ! call polympo_getMeshVtxMass(mpMesh,nVertices,c_loc(meshVtxMass))
      ! call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2Elm))

      ! do i = 1, numMPs
      !   vID = verticesOnCell(1,mp2Elm(i))
      !   call assert(meshVtxMass(vID) < TEST_VAL+TOLERANCE .and. meshVtxMass(vID) > TEST_VAL-TOLERANCE, "Error: wrong vtx mass")

      !   call assert(meshElmMass(mp2Elm(i)) < TEST_VAL+TOLERANCE &
      !         .and. meshElmMass(mp2Elm(i)) > TEST_VAL-TOLERANCE, "Error: wrong elm mass")
      ! end do
    end do

    deallocate(mpMass)
    deallocate(mpVel)
    deallocate(meshVtxMass)
    deallocate(meshElmMass)
  end subroutine

  subroutine runApiTest(mpMesh, numMPs, nVertices, nCells, numPush, mpLatLon, mpPosition, xVertex, yVertex, zVertex, latVertex)
    use :: polympo
    use :: readMPAS
    use :: iso_c_binding
    implicit none
    type(c_ptr) :: mpMesh
    integer :: i, j, k, numMPs, nCells, numPush, nVertices, nCompsDisp
    real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition, mpLatLon, mpMass, mpVel, dispIncr
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    real(kind=MPAS_RKIND), dimension(:), pointer :: meshVtxMass, meshElmMass, meshVtxVel
    real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex
    real(kind=MPAS_RKIND) :: TEST_VAL = 1.1_MPAS_RKIND

    nCompsDisp = 2
    allocate(dispIncr(nCompsDisp,nVertices))
    allocate(mpMass(1,numMPs))
    allocate(mpVel(2,numMPs))
    allocate(meshVtxMass(nVertices))
    allocate(meshElmMass(nCells))
    allocate(meshVtxVel(nVertices))

    dispIncr = TEST_VAL
    mpMass = TEST_VAL
    mpVel = TEST_VAL
    meshVtxMass = TEST_VAL
    meshElmMass = TEST_VAL
    meshVtxVel = TEST_VAL
    
    do j = 1, numPush
      call polympo_setMPPositions(mpMesh,3,numMPs,c_loc(mpPosition))
      call polympo_setMeshVtxOnSurfDispIncr(mpMesh,nCompsDisp,nVertices,c_loc(dispIncr))
      call polympo_setMPMass(mpMesh,1,numMPs,c_loc(mpMass))
      call polympo_setMPVel(mpMesh,2,numMPs,c_loc(mpVel))
      call polympo_setMeshVtxCoords(mpMesh, nVertices, c_loc(xVertex), c_loc(yVertex), c_loc(zVertex))
      call polympo_setMeshVtxRotLat(mpMesh,nVertices,c_loc(latVertex))

      call polympo_getMPMass(mpMesh, 1, numMPs, c_loc(mpMass))
      call polympo_getMPVel(mpMesh, 2, numMPs, c_loc(mpVel))
      call polympo_getMeshElmMass(mpMesh,nCells,c_loc(meshElmMass))
      call polympo_getMeshVtxMass(mpMesh,nVertices,c_loc(meshVtxMass))
      call polympo_getMeshVtxVel(mpMesh,nVertices, c_loc(xVertex),c_loc(yVertex))
    end do

    deallocate(dispIncr)
    deallocate(mpMass)
    deallocate(mpVel)
    deallocate(meshVtxMass)
    deallocate(meshElmMass)
    deallocate(meshVtxVel)

  end subroutine

end module
!---------------------------------------------------------------------------
!> todo add a discription
!---------------------------------------------------------------------------
program main
  use :: polympo
  use :: readMPAS
  use :: iso_c_binding
  use :: advectionTests
  implicit none
  include 'mpif.h'

  !integer, parameter :: APP_RKIND = selected_real_kind(15)
  type(c_ptr) :: mpMesh
  integer :: ierr, self
  integer :: argc, i, j, arglen, k, m, mpsScaleFactorPerVtx, localNumMPs
  integer :: setMeshOption, setMPOption
  integer :: maxEdges, vertexDegree, nCells, nVertices
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  real(kind=MPAS_RKIND) :: xc, yc, zc, xMP, yMP, zMP, radius, lon
  real(kind=MPAS_RKIND) :: pi = 4.0_MPAS_RKIND*atan(1.0_MPAS_RKIND)
  character (len=2048) :: filename, input
  character (len=64) :: onSphere
  real(kind=MPAS_RKIND) :: sphereRadius
  integer, dimension(:), pointer :: nEdgesOnCell
  real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: xCell, yCell, zCell
  integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
  integer :: numMPs, numMPsCount, numPush
  integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive, mp2Elm_new
  real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition, mpLatLon, mpPositions_new, mpLatLon_new
  integer, parameter :: MP_ACTIVE = 1
  integer, parameter :: MP_INACTIVE = 0
  integer, parameter :: INVALID_ELM_ID = -1

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle)
  call polympo_initialize()
  call polympo_enableTiming()

  call polympo_checkPrecisionForRealKind(MPAS_RKIND)
  argc = command_argument_count()
  if(argc == 3) then
    call get_command_argument(1, input)
    read(input, '(I7)') mpsScaleFactorPerVtx
    call get_command_argument(2, input)
    read(input, '(I7)') numPush
    call get_command_argument(3, filename)
  else
    write(0, *) "Usage: ./testFortranMPAdvection <mpsScaleFactorPerVtx> <numPush> <path to the nc file>"
  end if

  call readMPASMeshFromNCFile(filename, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
                        xCell, yCell, zCell, &
                        verticesOnCell, cellsOnCell)
  if (onSphere .ne. 'YES') then
    write (*,*) "The mesh is not spherical!"
    call exit(1)
  end if

  setMeshOption = 0 !create an empty mesh
  setMPOption = 0   !create an empty set of MPs
  mpMesh = polympo_createMPMesh(setMeshOption,setMPOption) !creates test mesh
  call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, &
                        xCell, yCell, zCell, &
                        verticesOnCell, cellsOnCell)

  !createMPs
  numMPs = 0
  do i = 1, nCells
    numMPs = numMPs + nEdgesOnCell(i) * mpsScaleFactorPerVtx
  end do

  print *, "Scale Factor", mpsScaleFactorPerVtx
  print *, "NUM MPs", numMPs

  allocate(mpsPerElm(nCells))
  allocate(mp2Elm(numMPs))
  allocate(mp2Elm_new(numMPs))

  allocate(isMPActive(numMPs))
  allocate(mpPosition(3,numMPs))
  allocate(mpPositions_new(3, numMPs))
  allocate(mpLatLon(2,numMPs))
  allocate(mpLatLon_new(2,numMPs))

  isMPActive = MP_ACTIVE !all active MPs and some changed below

  numMPsCount = 0
  do i = 1, nCells
    localNumMPs = nEdgesOnCell(i) * mpsScaleFactorPerVtx
    mp2Elm(numMPsCount+1:numMPsCount+localNumMPs) = i
    mpsPerElm(i) = localNumMPs
    numMPsCount = numMPsCount + localNumMPs
  end do

  call assert(numMPsCount == numMPs, "num mps miscounted")

  numMPsCount = 0
  do i = 1, nCells
    xc = 0.0_MPAS_RKIND
    yc = 0.0_MPAS_RKIND
    zc = 0.0_MPAS_RKIND
    do k = 1, nEdgesOnCell(i)
      j = verticesOnCell(k,i)
      xc = xc + xVertex(j) 
      yc = yc + yVertex(j) 
      zc = zc + zVertex(j) 
    end do
    xc = xc/nEdgesOnCell(i)
    yc = yc/nEdgesOnCell(i)
    zc = zc/nEdgesOnCell(i)

    do k = 1, nEdgesOnCell(i)
      j = verticesOnCell(k,i)
      
      ! note: m=0 not included but should lead to x_mp=xc and m=mpsScaleFactorPerVtx+1 is also not include but should lead to x_mp=xVertex(j)
      ! so'mpsScaleFactorPerVtx+1' segments and  'mpsScaleFactorPerVtx+2' points along line from 'xc' to 'xVertex(j)'
      ! taking only "inner" or interior 'mpsScaleFactorPerVtx' points (i.e., exclude end points of 'xc' and 'xVertex(j)') and same applies to y- and z-coordinates
      do m = 1, mpsScaleFactorPerVtx
        numMPsCount = numMPsCount + 1
        xMP = (mpsScaleFactorPerVtx+1 - m) * xc + m*xVertex(j) / (mpsScaleFactorPerVtx+1) ! linear interpolation
        yMP = (mpsScaleFactorPerVtx+1 - m) * yc + m*yVertex(j) / (mpsScaleFactorPerVtx+1) ! linear interpolation
        zMP = (mpsScaleFactorPerVtx+1 - m) * zc + m*zVertex(j) / (mpsScaleFactorPerVtx+1) ! linear interpolation
        
        ! normalize to project each MP to be on sphere of radius 'sphereRadius'
        radius = sqrt(xMP*xMP + yMP*yMP + zMP*zMP) ! assuming sphere center to be at origin
        xMP = xMP/radius * sphereRadius
        yMP = yMP/radius * sphereRadius
        zMP = zMP/radius * sphereRadius
        mpPosition(1,numMPsCount) = xMP
        mpPosition(2,numMPsCount) = yMP
        mpPosition(3,numMPsCount) = zMP
        mpLatLon(1,numMPsCount) = asin(zMP/sphereRadius)
        lon = atan2(yMP,xMP)
        if (lon .le. 0.0_MPAS_RKIND) then ! lon[0,2pi]
          lon = lon + 2.0_MPAS_RKIND*pi
        endif 
        mpLatLon(2,numMPsCount) = lon
      end do
    end do
  end do

  call assert(numMPsCount == numMPs, "num mps miscounted")

  call polympo_createMPs(mpMesh,nCells,numMPs,c_loc(mpsPerElm),c_loc(mp2Elm),c_loc(isMPActive))
  call polympo_setMPRotLatLon(mpMesh,2,numMPs,c_loc(mpLatLon))
  call polympo_setMPPositions(mpMesh,3,numMPs,c_loc(mpPosition))

  !call runAdvectionTest(mpMesh, numPush, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)

  call runAdvectionTest2(mpMesh, numPush, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)

  call polympo_getMPPositions(mpMesh, 3, numMPs, c_loc(mpPositions_new))
  call polympo_getMPRotLatLon(mpMesh, 2, numMPs, c_loc(mpLatLon_new))
  call polympo_getMPCurElmID(mpMesh, numMPS, c_loc(mp2Elm_new))
  do i = 1, numMPs
    PRINT *, "Difference: ", i, mpLatLon_new(1,i)-mpLatLon(1,i), mpLatLon_new(2,i)-mpLatLon(2,i)
  end do
  
  call runReconstructionTest(mpMesh, numMPs, numPush, nCells, nVertices, mp2Elm, &
                                   latVertex, lonVertex, nEdgesOnCell, verticesOnCell, sphereRadius)

  !call runApiTest(mpMesh, numMPs, nVertices, nCells, numPush, mpLatLon, mpPosition, xVertex, yVertex, zVertex, latVertex)

  call polympo_summarizeTime();

  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  deallocate(nEdgesOnCell)
  deallocate(xVertex)
  deallocate(yVertex)
  deallocate(zVertex)
  deallocate(latVertex)
  deallocate(lonVertex)
  deallocate(xCell)
  deallocate(yCell)
  deallocate(zCell)
  deallocate(verticesOnCell)
  deallocate(cellsOnCell)
  deallocate(mpsPerElm)
  deallocate(mp2Elm)
  deallocate(isMPActive)
  deallocate(mpPosition)
  deallocate(mpLatLon)

  stop

end program
