module MYDATA_QUEUE
    include "queue.f90"
end module MYDATA_QUEUE

program main
    use MYDATA_QUEUE
    type (QUEUE_STRUCT), pointer :: queue
    logical :: success
    integer :: result

    call queue_create(queue, 5)
    call queue_append_data( queue, 122, success )
    call queue_append_data( queue, 233, success )
    call queue_append_data( queue, 344, success )
    call queue_append_data( queue, 455, success )
    call queue_append_data( queue, 566, success )

    print *, queue_retrieve_data(queue)
    print *, queue_retrieve_data(queue)
    print *, queue_retrieve_data(queue)

    call queue_destroy(queue)
end program