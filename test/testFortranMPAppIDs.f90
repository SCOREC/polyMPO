module MYDATA_QUEUE
    include "queue.f90"

    subroutine create_app_ids(queue)
        type (QUEUE_STRUCT), pointer :: queue
        logical :: success

        call queue_create(queue, 5)
        call queue_append_data( queue, 122, success )
        call queue_append_data( queue, 233, success )
        call queue_append_data( queue, 344, success )
        call queue_append_data( queue, 455, success )
        call queue_append_data( queue, 566, success )
    end subroutine

    integer function GetAppID(queueIn) result(id)
        use :: iso_c_binding
        type(c_ptr), value :: queueIn
        type (QUEUE_STRUCT), pointer :: queue
        call c_f_pointer(queueIn, queue)
        id = queue_retrieve_data(queue)
    end function
end module MYDATA_QUEUE

program main
    use MYDATA_QUEUE
    use :: polympo
    use :: iso_c_binding
    implicit none
    include 'mpif.h'

    interface
        subroutine testAppIDPointer(mpMesh, numMPs, mp2Elm) &
        bind(C, NAME='testAppIDPointer')
            use :: iso_c_binding
            type(c_ptr), value :: mpMesh
            integer(c_int), value :: numMPs
            type(c_ptr), value :: mp2Elm
        end subroutine
    end interface

    type (QUEUE_STRUCT), pointer :: queue
    integer :: setMeshOption, setMPOption
    integer :: ierr, self
    integer :: mpi_comm_handle = MPI_COMM_WORLD
    integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive
    type(c_ptr) :: mpMesh
    integer :: nCells, numMPs, appID
    integer, parameter :: MP_ACTIVE = 1

    ! Initialize
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_handle, self, ierr)
    call polympo_setMPICommunicator(mpi_comm_handle)
    call polympo_initialize()
    ! Create Queue
    call create_app_ids(queue)

    ! Create Mesh
    setMeshOption = 1 !create a hard coded planar test mesh
    setMPOption = 0   !create an empty set of MPs
    mpMesh = polympo_createMPMesh(setMeshOption, setMPOption)

    numMPs = 1
    allocate(mpsPerElm(nCells))
    allocate(mp2Elm(numMPs))
    allocate(isMPActive(numMPs))
    nCells = polympo_getMeshNumElms(mpMesh)
    mpsPerElm = 1
    mp2Elm = 1
    isMPActive = MP_ACTIVE
    call polympo_createMPs(mpMesh, nCells, numMPs, c_loc(mpsPerElm), c_loc(mp2Elm), c_loc(isMPActive))
    
    ! Set function and opaque data structure(list/queue) used to retrieve appIDS
    call polympo_setAppIDFunc(mpMesh, c_funloc(GetAppID), c_loc(queue));
    
    ! Test to emulate using appID function during migration
    call testAppIDPointer(mpMesh, numMPs, c_loc(mp2Elm))

    ! Clean Up
    call polympo_deleteMPMesh(mpMesh)
    call queue_destroy(queue)
    call polympo_finalize()
    call mpi_finalize(ierr)

end program