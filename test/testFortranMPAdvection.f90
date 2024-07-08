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
  integer :: argc, i, j, arglen, k, mpsPerEdge, localNumMPs
  integer :: setMeshOption, setMPOption
  integer :: maxEdges, vertexDegree, nCells, nVertices
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  real(kind=MPAS_RKIND) :: xc, yc, zc, radiusX, radiusY, radiusZ, lon
  real(kind=MPAS_RKIND) :: pi = 4.0_MPAS_RKIND*atan(1.0_MPAS_RKIND)
  character (len=2048) :: filename
  character (len=64) :: onSphere
  real(kind=MPAS_RKIND) :: sphereRadius, xComputed, yComputed, zComputed, latComputed, lonComputed
  integer, dimension(:), pointer :: nEdgesOnCell
  real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
  integer :: numMPs 
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
  mpsPerEdge = 5
  numMPs = maxEdges * mpsPerEdge
  allocate(mpsPerElm(nCells))
  allocate(mp2Elm(numMPs))
  allocate(isMPActive(numMPs))
  allocate(mpPosition(3,numMPs))
  allocate(mpLatLon(2,numMPs))

  isMPActive = MP_ACTIVE !all active MPs and some changed below

  numMPs = 1
  do i = 1, nCells
    localNumMPs = nEdgesOnCell(i) * mpsPerEdge
    mp2Elm(numMPs:numMPs+localNumMPs) = i
    mpsPerElm(i) = localNumMPs
    numMPs = numMPs + localNumMPs
  end do

  do i = 1, numMPs
    xc = 0.0_MPAS_RKIND
    yc = 0.0_MPAS_RKIND
    zc = 0.0_MPAS_RKIND
    do k = 1, nEdgesOnCell(i)
      j = verticesOnCell(k,i)
      xc = xc + xVertex(j) 
      yc = yc + yVertex(j) 
      zc = zc + zVertex(j) 
      xComputed = sphereRadius*cos(latVertex(j))*cos(lonVertex(j))
      yComputed = sphereRadius*cos(latVertex(j))*sin(lonVertex(j))
      zComputed = sphereRadius*sin(latVertex(j))
      latComputed = asin(zVertex(j)/sphereRadius)
      lonComputed = atan2(yVertex(j),xVertex(j))
      if (lonComputed .le. 0.0_MPAS_RKIND) then ! lon[0,2pi]
        lonComputed = lonComputed + 2.0_MPAS_RKIND*pi
      endif

    end do
    xc = xc/nEdgesOnCell(i)
    yc = yc/nEdgesOnCell(i)
    zc = zc/nEdgesOnCell(i)

    do k = 1, nEdgesOnCell(i)
      j = verticesOnCell(k,i)
      radiusX = (mpsPerEdge+1) * (k - xc) / (xVertex(j) - xc)  ! linear interpolation
      radiusY = (mpsPerEdge+1) * (k - yc) / (yVertex(j) - yc)  ! linear interpolation
      radiusZ = (mpsPerEdge+1) * (k - zc) / (zVertex(j) - zc)  ! linear interpolation
      
      xc = xc/radiusX * sphereRadius
      yc = yc/radiusY * sphereRadius
      zc = zc/radiusZ * sphereRadius
      mpPosition(1,i) = xc
      mpPosition(2,i) = yc
      mpPosition(3,i) = zc
      mpLatLon(1,i) = asin(zc/sphereRadius)
      lon = atan2(yc,xc)
      if (lon .le. 0.0_MPAS_RKIND) then ! lon[0,2pi]
        lon = lon + 2.0_MPAS_RKIND*pi
      endif 
      mpLatLon(2,i) = lon
    end do
  end do

  call polympo_createMPs(mpMesh,nCells,numMPs,c_loc(mpsPerElm),c_loc(mp2Elm),c_loc(isMPActive))
  call polympo_setMPRotLatLon(mpMesh,2,numMPs,c_loc(mpLatLon))
  call polympo_setMPPositions(mpMesh,3,numMPs,c_loc(mpPosition))

  call calculateSurfaceDisplacement(mpMesh, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius)
  call polympo_push(mpMesh) ! TODO: preform multiple times configurable (beta)
  ! TODO: add timer 
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
