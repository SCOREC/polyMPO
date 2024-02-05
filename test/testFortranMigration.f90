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

    ! Initialize
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_handle, self, ierr)
    call polympo_setMPICommunicator(mpi_comm_handle)
    call polympo_initialize()

    ! Create Mesh
    setMeshOption = 0 !create an empty mesh
    setMPOption = 0   !create an empty set of MPs
    mpMesh = polympo_createMPMesh(setMeshOption, setMPOption)

    ! call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
    !                     nCells, nVertices, nEdgesOnCell, &
    !                     onSphere, sphereRadius, &
    !                     xVertex, yVertex, zVertex, &
    !                     latVertex, &
    !                     verticesOnCell, cellsOnCell)


    ! Clean Up
    deallocate(nEdgesOnCell)
    deallocate(xVertex)
    deallocate(yVertex)
    deallocate(zVertex)
    deallocate(verticesOnCell)
    deallocate(cellsOnCell)

    call polympo_deleteMPMesh(mpMesh)
    call polympo_finalize()
    call mpi_finalize(ierr)
end program