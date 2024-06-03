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
  integer :: argc, i, j, arglen, k
  integer :: setMeshOption, setMPOption
  integer :: maxEdges, vertexDegree, nCells, nVertices
  integer :: nCompsDisp
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  real(kind=MPAS_RKIND) :: pi = 4.0_MPAS_RKIND*atan(1.0_MPAS_RKIND)
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
  real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpMass, mpVel
  real(kind=MPAS_RKIND), dimension(:), pointer :: meshMass, meshVelU, meshVelV
  logical :: inBound
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
    write(0, *) "Usage: ./testFortranMPReconstruction <path to the nc file>"
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
  sphereRadius = 1
  mpMesh = polympo_createMPMesh(setMeshOption,setMPOption) !creates test mesh
  call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, &
                        verticesOnCell, cellsOnCell)
 
  nCompsDisp = 2
  allocate(dispIncr(nCompsDisp,nVertices))
  !createMPs
  numMPs = nCells
  allocate(mpsPerElm(nCells))
  allocate(mp2Elm(numMPs))
  allocate(isMPActive(numMPs))
  allocate(mpPosition(3,numMPs))
  allocate(mpLatLon(2,numMPs))
  allocate(mpMass(1,numMPs))
  allocate(mpVel(2,numMPs))
  allocate(meshMass(nVertices))
  allocate(meshVelU(nVertices))
  allocate(meshVelV(nVertices))

  isMPActive = MP_ACTIVE !all active MPs and some changed below
  mpsPerElm = 1 !all elements have 1 MP and some changed below
  mpMass = 1.1
  mpVel = 1.1
  do i = 1,numMPs
    mp2Elm(i) = i
    j = verticesOnCell(1,i)
    mpLatLon(1,j) = latVertex(j)
    mpLatLon(2,j) = lonVertex(j) 
    mpPosition(1,j) = xVertex(j)
    mpPosition(2,j) = yVertex(j)
    mpPosition(3,j) = zVertex(j)
    !write(*,*) i, mpVel(1,i), mpVel(2,i)
  end do
  !write(*,*) mpVel(1,1), mpVel(2,1)
  !write(*,*) mpVel(1,2), mpVel(2,2)
  !write(*,*) mpVel(1,3), mpVel(2,3)

  call polympo_createMPs(mpMesh,nCells,numMPs,c_loc(mpsPerElm),c_loc(mp2Elm),c_loc(isMPActive))
  call polympo_setMPRotLatLon(mpMesh,2,numMPs,c_loc(mpLatLon))
  call polympo_setMPPositions(mpMesh,3,numMPs,c_loc(mpPosition))
 
  !write(*,*) mpMass 
  !write(*,*) mpVel
  call polympo_setMPMass(mpMesh,1,numMPs,c_loc(mpMass))
  call polympo_setMPVel(mpMesh,2,numMPs,c_loc(mpVel))
  !mpMass = 0
  !mpVel = -1
  !call polympo_getMPMass(mpMesh,1,numMPs,c_loc(mpMass))
  !write(*,*) mpMass 
  !call polympo_getMPVel(mpMesh,2,numMPs,c_loc(mpVel))
  !write(*,*) mpVel 

  call polympo_setReconstructionOfVel(mpMesh,0,polympo_getMeshFVtxType())
  call polympo_applyReconstruction(mpMesh)

  meshVelU = -1
  meshVelV = -1
  !call polympo_getMeshVtxVel(mpMesh,nCells,c_loc(meshMass)) todo
  call polympo_getMeshVtxVel(mpMesh,nVertices,c_loc(meshVelU),c_loc(meshVelV))
  do i = 1,numMPs
    j = verticesOnCell(1,i)
    write(*,*) mpVel(1,i), meshVelU(j)
    write(*,*) mpVel(1,i), meshVelV(j)
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
  deallocate(dispIncr)
  deallocate(mpsPerElm)
  deallocate(mp2Elm)
  deallocate(isMPActive)
  deallocate(mpPosition)
  deallocate(mpLatLon)
  deallocate(mpMass)
  deallocate(mpVel)
  deallocate(meshMass)
  deallocate(meshVelU)
  deallocate(meshVelV)

  stop
end program
