! Adapted from : https://github.com/mapmeld/fortran-machine/blob/5e5e6bbdab353149047971a3f872238a2a4b716e/flibs-0.9/flibs/src/datastructures/queues.f90

! queues.f90 --
!     Include file for defining queues with a fixed capacity
!     Queues as implemented here are simply arrays where
!     data are inserted at the end and retrieved from the
!     top.
!
type QUEUE_STRUCT
logical                                 :: full
integer                                 :: start
integer                                 :: end
integer, dimension(:), pointer          :: data
end type QUEUE_STRUCT

!
! Define the subroutines and functions
!
contains

! queue_create --
!     Create and initialise a queue
! Arguments:
!     queue      Pointer to new queue
!     capacity   The number of data that can be stored
! Note:
!     There is no check that the capacity is positive!
!
subroutine queue_create( queue, capacity )
type(QUEUE_STRUCT), pointer        :: queue
integer                            :: capacity

allocate( queue )
allocate( queue%data(0:capacity) )

queue%full  = .false.
queue%start = 1
queue%end   = 0
end subroutine queue_create

! queue_destroy --
!     Destroy a queue
! Arguments:
!     queue       Pointer to the queue to be destroyed
!
subroutine queue_destroy( queue )
type(QUEUE_STRUCT), pointer  :: queue

deallocate( queue%data )
deallocate( queue )
end subroutine queue_destroy

! queue_empty --
!     Check if the queue is empty
! Arguments:
!     queue       Pointer to the queue
! Result:
!     logical indicating if the queue is
!     empty or not
!
logical function queue_empty( queue )
type(QUEUE_STRUCT), intent(in)  :: queue

queue_empty = (.not. queue%full) .and. &
    queue%end .eq. queue%start - 1

end function queue_empty

! queue_retrieve_data
!     Return the data stored at the top,
!     and remove them from the queue
! Arguments:
!     queue      Queue to be examined
! Result:
!     Data stored at the top, afterwards
!     removed
!
function queue_retrieve_data( queue ) result(data)
type(QUEUE_STRUCT)             :: queue
integer                        :: data

if (queue_empty(queue)) then
    print *, "ERROR: Queue is empty"
    call exit
endif

data = queue%data(queue%start)

if ( .not. queue_empty(queue) ) then
    queue%start = queue%start + 1
    if ( queue%start .gt. size(queue%data) ) then
        queue%start = 1
    endif
    queue%full = .false.
endif

end function queue_retrieve_data

! queue_append_data
!     Append data to the end of the queue
! Arguments:
!     queue      Queue to which to add the data
!     data       The data to be added
!
subroutine queue_append_data( queue, data)
type(QUEUE_STRUCT)           :: queue
integer, intent(in) :: data

if (queue%full) then
    print *, "ERROR: Queue is full"
    call exit
endif

queue%end = queue%end + 1
if ( queue%end .gt. size(queue%data) ) then
    queue%end = 1
endif
if ( queue%start .eq. queue%end+1 ) then
    queue%full = .true.
endif
if ( queue%end   .eq. size(queue%data) .and. &
        queue%start .eq. 1 ) then
    queue%full = .true.
endif
queue%data(queue%end) = data

end subroutine queue_append_data