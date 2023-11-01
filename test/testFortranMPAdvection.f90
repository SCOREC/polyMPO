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

  !integer, parameter :: APP_RKIND = selected_real_kind(15)
  integer :: ierr, self
  integer :: argc, i, j, arglen, k
  integer :: setMeshOption, setMPOption
  integer :: maxEdges, vertexDegree, nCells, nVertices
  integer :: nCompsDisp
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  real(kind=MPAS_RKIND) :: xc, yc, zc, radius, maxlon, minlon, deltaLon, lon
  real(kind=MPAS_RKIND) :: pi = 4*atan(1.0)
  character (len=2048) :: filename
  real(kind=MPAS_RKIND), dimension(:,:), pointer :: dispIncr
  character (len=64) :: onSphere
  real(kind=MPAS_RKIND) :: sphereRadius, xComputed, yComputed, zComputed, latComputed, lonComputed
  integer, dimension(:), pointer :: nEdgesOnCell
  real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
  integer :: numMPs 
  integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive
  real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition, mpLatLon

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
  if (onSphere .ne. 'YES') then
    write (*,*) "The mesh is not spherical!"
    call exit(1)
  end if
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
  allocate(mpPosition(3,numMPs))
  allocate(mpLatLon(2,numMPs))
  isMPActive = 1 !no inactive MPs and some changed below
  mpsPerElm = 1 !all elements have 1 MP and some changed below
  do i = 1,numMPs
    mp2Elm(i) = i
  end do
  call polympo_createMPs(mpMesh,nCells,numMPs,c_loc(mpsPerElm),c_loc(mp2Elm),c_loc(isMPActive))
  do i = 1, nCells
    if (.true.) then
      xc = 0.0_MPAS_RKIND
      yc = 0.0_MPAS_RKIND
      zc = 0.0_MPAS_RKIND
      do k = 1, nEdgesOnCell(i)
        j = verticesOnCell(k,i)
        xc = xc + xVertex(j) 
        yc = yc + yVertex(j) 
        zc = zc + zVertex(j) 
        !xComputed = sphereRadius*cos(lonVertex(j))*cos(latVertex(j))
        !yComputed = sphereRadius*sin(lonVertex(j))*cos(latVertex(j))
        !zComputed = sphereRadius*sin(latVertex(j))
        xComputed = sphereRadius*cos(latVertex(j))*cos(lonVertex(j))
        yComputed = sphereRadius*cos(latVertex(j))*sin(lonVertex(j))
        zComputed = sphereRadius*sin(latVertex(j))
        write(*,*)  xVertex(j), xComputed
        write(*,*)  yVertex(j), yComputed
        write(*,*)  zVertex(j), zComputed
        latComputed = asin(zVertex(j)/sphereRadius)
        lonComputed = atan2(yVertex(j),xVertex(j))
        if (lonComputed .le. 0.0) then ! lon[0,2pi]
          lonComputed = lonComputed + 2*pi
        endif
        write(*,*)  latVertex(j), latComputed
        write(*,*)  lonVertex(j), lonComputed

      end do
      xc = xc/nEdgesOnCell(i)
      yc = yc/nEdgesOnCell(i)
      zc = zc/nEdgesOnCell(i)
      ! normalize
      radius = sqrt(xc*xc + yc*yc + zc*zc)! assuming sphere center to be at origin
      xc = xc/radius * sphereRadius
      yc = yc/radius * sphereRadius
      zc = zc/radius * sphereRadius
      mpPosition(1,i) = xc
      mpPosition(2,i) = yc
      mpPosition(3,i) = zc
      mpLatLon(1,i) = asin(zc/sphereRadius)
      lon = atan2(yc,xc)
      if (lon .le. 0.0) then ! lon[0,2pi]
        lon = lon + 2*pi
      endif 
      mpLatLon(2,i) = lon
    endif
  end do
  ! check first element/cell for delta
  maxlon = minval(lonVertex)
  minlon = maxval(lonVertex)
  do i = 1, nEdgesOnCell(1)
    j = verticesOnCell(i,1)
    if(maxlon .lt. lonVertex(j)) then
      maxlon = lonVertex(j)
    endif
    if(minlon .gt. lonVertex(j)) then
      minlon = lonVertex(j)
    endif
  end do
  !todo setLatLonMPPositions
  !todo setXYZMPPositions
  
  deltaLon = maxlon - minlon
  !todo getLatLonMPPositions
  do i = 1,numMPs
    ! R*cos(latMP) *deltaLon
    dispIncr(1,i) = 0!sphereRadius*cos(
    dispIncr(2,i) = deltaLon
  end do
  call polympo_setMeshOnSurfDispIncr(mpMesh, nCompsDisp, nVertices, c_loc(dispIncr))
  call polympo_push(mpMesh)
 
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
  deallocate(dispIncr)
  deallocate(mpsPerElm)
  deallocate(mp2Elm)
  deallocate(isMPActive)
  deallocate(mpPosition)
  deallocate(mpLatLon)

  stop
end program
