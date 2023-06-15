program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: nverts = 19 !matches the number of verts in the test mesh
  integer :: numComps = 2 !matches the Velocity coordinates mesh array
  integer :: numMPs = 51 !matches the MPs created in createMpMesh
  real(c_double), dimension(:), allocatable :: MParray
  real(c_double), dimension(:), allocatable :: Mesharray
  integer :: ierr, self
  type(c_ptr) :: mpMesh
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)
  call polympo_initialize()
  mpMesh = polympo_createMpMesh() !creates test mesh
  allocate(Mesharray(nverts*numComps))
  Mesharray = 1
  write (*,*) Mesharray
  allocate(MParray(numMPs*numComps))
  MParray = 1
  write (*,*) MParray
  call polympo_setMPVelArray(mpMesh, numMPs, MParray);
  call polympo_setMeshVelArray(mpMesh, nverts, Mesharray);
  call polympo_getMPVelArray(mpMesh, numMPs, MParray);
  call polympo_getMeshVelArray(mpMesh, nverts, Mesharray);
  write (*,*) MParray
  deallocate(MParray)
  write (*,*) Mesharray
  deallocate(Mesharray)
  call polympo_deleteMpMesh(mpMesh)
  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
