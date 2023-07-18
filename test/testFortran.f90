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

  integer(c_int) :: nverts 
  integer :: numComps
  integer :: numMPs 
  integer :: i, j
  integer :: setMeshOption, setMPOption
  real(c_double) :: value1, value2
  real(c_double), dimension(:), allocatable :: MParray
  real(c_double), dimension(:,:), allocatable, target :: MP2dArray
  real(c_double), dimension(:), allocatable :: Mesharray
  integer :: ierr, self
  type(c_ptr) :: mpMesh

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, self, ierr)

  call polympo_setMPICommunicator(MPI_COMM_WORLD) !this is not supported yet! only for showing
  call polympo_initialize()

  setMeshOption = 1
  setMPOption = 1    
  mpMesh = polympo_createMPMesh(setMeshOption,setMPOption) !creates test mesh
  
  nverts = 19 !todo use getNumVtx from the Mesh object
  numComps = 2 !todo use getNumComps from velocity fields
  numMPs = 49 !todo use getNumMPs from the MaterialPoints object

  allocate(Mesharray(nverts*numComps))
  allocate(MParray(numMPs*numComps))
  allocate(MP2dArray(numMPs,numComps))

  do i = 1,numMPs
    do j = 1,numComps
      MP2dArray(i,j) = (i-1)*numComps+(j-1)
    end do
  end do
  call polympo_setMP2dVelArray(mpMesh, numMPs, numComps, c_loc(MP2dArray))

  value1 = 42
  MParray = value1
  call polympo_setMPVelArray(mpMesh, numMPs, MParray)

  Mesharray = value1
  call polympo_setMeshVelArray(mpMesh, nverts, Mesharray)

  MParray = 1
  call polympo_getMPVelArray(mpMesh, numMPs, MParray)
  call assert(all(MParray .eq. value1), "Assert MParray == value1 Failed!")

  Mesharray = 1
  call polympo_getMeshVelArray(mpMesh, nverts, Mesharray)
  call assert(all(Mesharray .eq. value1), "Assert Mesharray == value1 Failed!")

  value2 = 24
  MParray = value2
  call polympo_setMPVelArray(mpMesh, numMPs, MParray)

  MParray = 1
  call polympo_getMPVelArray(mpMesh, numMPs, MParray)
  call assert(all(MParray .eq. value2), "Assert MParray == value2 Failed!")

  Mesharray = value2
  call polympo_setMeshVelArray(mpMesh, nverts, Mesharray)

  Mesharray = 1
  call polympo_getMeshVelArray(mpMesh, nverts, Mesharray)
  call assert(all(Mesharray .eq. value2), "Assert Mesharray == value2 Failed!")

  deallocate(MParray)
  deallocate(Mesharray)

  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  stop
end
