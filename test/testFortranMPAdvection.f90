!---------------------------------------------------------------------------
!> todo add a discription
!---------------------------------------------------------------------------
program main
  use :: polympo
  use :: readMPAS
  use :: iso_c_binding
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
  character (len=2048) :: filename
  character (len=64) :: onSphere
  real(kind=MPAS_RKIND) :: sphereRadius
  integer, dimension(:), pointer :: nEdgesOnCell
  real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
  integer :: numMPs, numMPsCount, numPush
  integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive
  real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition, mpLatLon
  integer, parameter :: MP_ACTIVE = 1
  integer, parameter :: MP_INACTIVE = 0
  integer, parameter :: INVALID_ELM_ID = -1

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle)
  call polympo_initialize()

  call polympo_checkPrecisionForRealKind(MPAS_RKIND)
  argc = command_argument_count()
  if(argc == 1) then
    call get_command_argument(1, filename)
  else
    write(0, *) "Usage: ./testFortranInterpolatePush <path to the nc file>"
  end if

  call readMPASMeshFromNCFile(filename, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
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
                        verticesOnCell, cellsOnCell)

  !createMPs
  numMPs = 0
  mpsScaleFactorPerVtx = 5
  do i = 1, nCells
    numMPs = numMPs + nEdgesOnCell(i) * mpsScaleFactorPerVtx
  end do
  allocate(mpsPerElm(nCells))
  allocate(mp2Elm(numMPs))
  allocate(isMPActive(numMPs))
  allocate(mpPosition(3,numMPs))
  allocate(mpLatLon(2,numMPs))

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

  numPush = 5
  do i = 1, numPush
    call calcSurfDispIncr(mpMesh, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)
    call polympo_push(mpMesh)
    ! TODO: add timer 
  end do

  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  deallocate(nEdgesOnCell)
  deallocate(xVertex)
  deallocate(yVertex)
  deallocate(zVertex)
  deallocate(latVertex)
  deallocate(lonVertex)
  deallocate(verticesOnCell)
  deallocate(cellsOnCell)
  deallocate(mpsPerElm)
  deallocate(mp2Elm)
  deallocate(isMPActive)
  deallocate(mpPosition)
  deallocate(mpLatLon)

  stop

  contains
  include "calculateDisplacement.f90"

end program
