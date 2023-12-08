module MYDATA_QUEUE
    include "queue.f90"
end module MYDATA_QUEUE

program main
    use MYDATA_QUEUE
    use :: polympo
    use :: iso_c_binding
    implicit none
    include 'mpif.h'

    type (QUEUE_STRUCT), pointer :: queue
    logical :: success
    integer :: setMeshOption, setMPOption
    integer :: ierr, self
    integer :: mpi_comm_handle = MPI_COMM_WORLD
    type(c_ptr) :: mpMesh

    ! Initialize
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_handle, self, ierr)
    call polympo_setMPICommunicator(mpi_comm_handle)
    call polympo_initialize()
    ! Create Queue
    call queue_create(queue, 5)
    call queue_append_data( queue, 122, success )
    call queue_append_data( queue, 233, success )
    call queue_append_data( queue, 344, success )
    call queue_append_data( queue, 455, success )
    call queue_append_data( queue, 566, success )

    print *, queue_retrieve_data(queue)
    print *, queue_retrieve_data(queue)
    print *, queue_retrieve_data(queue)
    ! Create Mesh
    setMeshOption = 1 !create a hard coded planar test mesh
    setMPOption = 0   !create an empty set of MPs
    mpMesh = polympo_createMPMesh(setMeshOption, setMPOption)
    ! Clean Up
    call polympo_deleteMPMesh(mpMesh)
    call queue_destroy(queue)
    call polympo_finalize()
    call mpi_finalize(ierr)
end program