program main
    use :: polympo
    use :: iso_c_binding
    implicit none
    include 'mpif.h'
    integer :: ierr, self
    integer :: mpi_comm_handle = MPI_COMM_WORLD
    type(c_ptr) :: mpMesh
    integer :: setMeshOption, setMPOption

    ! Initialize
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_handle, self, ierr)
    call polympo_setMPICommunicator(mpi_comm_handle)
    call polympo_initialize()

    ! Create Mesh
    setMeshOption = 1 !create a hard coded planar test mesh
    setMPOption = 0   !create an empty set of MPs
    mpMesh = polympo_createMPMesh(setMeshOption, setMPOption)

    ! Clean Up
    call polympo_deleteMPMesh(mpMesh)
    call polympo_finalize()
    call mpi_finalize(ierr)
end program