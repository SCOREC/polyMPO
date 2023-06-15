program main
  use polympo
  use iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: nverts = 19 !matches the number of verts in the test mesh
  integer :: numComps = 3 !matches the XYZ coordinates mesh array
  real(c_double), dimension(:), allocatable :: array
  integer :: ierr, self
  type(c_ptr) :: mpMesh
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)
  call polympo_initialize()
  mpMesh = polympo_createMpMesh() !creates test mesh
  allocate(array(nverts*numComps))
  array = 1
  write (*,*) array
  call polympo_setMPVelArray(mpMesh, nverts, array);
  call polympo_getMeshVelArray(mpMesh, nverts, array);
  !call polympo_getMeshCurPosXYZArray(mpMesh, nverts, array);
  write (*,*) array
  deallocate(array)
  call polympo_deleteMpMesh(mpMesh)
  call polympo_finalize()
  call mpi_finalize(ierr)
  stop
end
