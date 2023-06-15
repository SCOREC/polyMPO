subroutine assert(condition,message)
  implicit none
  logical :: condition
  character(*) :: message
  if (condition .neqv. .true.) then
    write (*,*) message
    call exit(1)
  endif
end subroutine

program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: nverts = 19 !matches the number of verts in the test mesh
  integer :: numComps = 2 !matches the Velocity coordinates mesh array
  integer :: numMPs = 51 !matches the MPs created in createMpMesh
  integer :: i, j
  real(c_double), dimension(:), allocatable :: MParray
  real(c_double), dimension(:), allocatable :: Mesharray
  integer :: ierr, self
  type(c_ptr) :: mpMesh
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)
  call polympo_initialize()
  mpMesh = polympo_createMpMesh() !creates test mesh
  allocate(Mesharray(nverts*numComps))
  Mesharray = 42
  allocate(MParray(numMPs*numComps))
  MParray = 42
  call polympo_setMPVelArray(mpMesh, numMPs, MParray);
  call polympo_setMeshVelArray(mpMesh, nverts, Mesharray);
  call polympo_getMPVelArray(mpMesh, numMPs, MParray);
  call polympo_getMeshVelArray(mpMesh, nverts, Mesharray);
  
  do, i=1,numMPs-1
    do, j=1,numComps-1
      call assert(MParray(i*numComps+j)==42, "Assert MParray == 42 Failed!")
    enddo
  enddo
  MParray = 24
  call polympo_setMPVelArray(mpMesh, numMPs, MParray);
  call polympo_getMPVelArray(mpMesh, numMPs, MParray);
  do, i=1,numMPs-1
    do, j=1,numComps-1
      call assert(MParray(i*numComps+j)==24, "Assert MParray == 24 Failed!")
    enddo
  enddo
  do, i=1,nverts-1
    do, j=1,numComps-1
      call assert(Mesharray(i*numComps+j)==42, "Assert Mesharray == 42 Failed!")
    enddo
  enddo
  Mesharray = 24
  call polympo_setMeshVelArray(mpMesh, nverts, Mesharray);
  call polympo_getMeshVelArray(mpMesh, nverts, Mesharray);
  do, i=1,nverts-1
    do, j=1,numComps-1
      call assert(Mesharray(i*numComps+j)==24, "Assert Mesharray == 24 Failed!")
    enddo
  enddo
  deallocate(MParray)
  deallocate(Mesharray)
  call polympo_deleteMpMesh(mpMesh)
  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
