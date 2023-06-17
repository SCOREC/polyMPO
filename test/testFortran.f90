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
  call polympo_setCommunicator(MPI_COMM_WORLD) !this is not supported yet! only for showing
  call polympo_initialize()
  mpMesh = polympo_createMpMesh() !creates test mesh

  allocate(Mesharray(nverts*numComps))
  allocate(MParray(numMPs*numComps))

  MParray = 42
  call polympo_setMPVelArray(mpMesh, numMPs, MParray);
  Mesharray = 42
  call polympo_setMeshVelArray(mpMesh, nverts, Mesharray);

  MParray = 1
  call polympo_getMPVelArray(mpMesh, numMPs, MParray);
  call assert(all(MParray .eq. 42), "Assert MParray == 42 Failed!")
  Mesharray = 1
  call polympo_getMeshVelArray(mpMesh, nverts, Mesharray);
  call assert(all(Mesharray .eq. 42), "Assert Mesharray == 42 Failed!")

  MParray = 24
  call polympo_setMPVelArray(mpMesh, numMPs, MParray);
  MParray = 1
  call polympo_getMPVelArray(mpMesh, numMPs, MParray);
  call assert(all(MParray .eq. 24), "Assert MParray == 24 Failed!")

  Mesharray = 24
  call polympo_setMeshVelArray(mpMesh, nverts, Mesharray);
  Mesharray = 1
  call polympo_getMeshVelArray(mpMesh, nverts, Mesharray);
  call assert(all(Mesharray .eq. 24), "Assert Mesharray == 24 Failed!")

  deallocate(MParray)
  deallocate(Mesharray)
  call polympo_deleteMpMesh(mpMesh)
  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
