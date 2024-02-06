program main
    use :: polympo
    use :: readMPAS
    use :: iso_c_binding
    implicit none
    include 'mpif.h'

    integer :: ierr, self
    integer :: mpi_comm_handle = MPI_COMM_WORLD
    type(c_ptr) :: mpMesh
    integer :: setMeshOption, setMPOption
    character (len=64) :: onSphere
    integer :: maxEdges, vertexDegree, nCells, nVertices
    real(kind=MPAS_RKIND) :: sphereRadius
    integer, dimension(:), pointer :: nEdgesOnCell
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
    integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive
    integer, parameter :: MP_ACTIVE = 1
    integer :: numMPs

    ! Initialize
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_handle, self, ierr)
    call polympo_setMPICommunicator(mpi_comm_handle)
    call polympo_initialize()

    ! Create Mesh
    setMeshOption = 0 !create an empty mesh
    setMPOption = 0   !create an empty set of MPs
    mpMesh = polympo_createMPMesh(setMeshOption, setMPOption)

    nCells = 1
    nVertices = 0
    maxEdges = 0
    vertexDegree = 0
    sphereRadius = 0

    allocate(xVertex(nVertices))
    allocate(yVertex(nVertices))
    allocate(zVertex(nVertices))
    allocate(latVertex(nVertices))
    allocate(nEdgesOnCell(nCells))
    allocate(cellsOnCell(maxEdges, nCells))
    allocate(verticesOnCell(maxEdges, nCells))

    ! call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
    !                     nCells, nVertices, nEdgesOnCell, &
    !                     onSphere, sphereRadius, &
    !                     xVertex, yVertex, zVertex, &
    !                     latVertex, &
    !                     verticesOnCell, cellsOnCell)


    numMPs = 1
    allocate(mpsPerElm(nCells))
    allocate(mp2Elm(numMPs))
    allocate(isMPActive(numMPs))
    mpsPerElm = 1
    mp2Elm = 1
    isMPActive = MP_ACTIVE
    ! call polympo_createMPs(mpMesh, nCells, numMPs, c_loc(mpsPerElm), c_loc(mp2Elm), c_loc(isMPActive))


    ! Clean Up
    deallocate(nEdgesOnCell)
    deallocate(xVertex)
    deallocate(yVertex)
    deallocate(zVertex)
    deallocate(latVertex)
    deallocate(cellsOnCell)
    deallocate(verticesOnCell)

    call polympo_deleteMPMesh(mpMesh)
    call polympo_finalize()
    call mpi_finalize(ierr)
end program