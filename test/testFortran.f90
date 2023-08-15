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
    
  integer, parameter :: APP_RKIND = selected_real_kind(15)
  integer :: nverts 
  integer :: numComps
  integer :: numMPs 
  integer :: i, j
  integer :: setMeshOption, setMPOption
  integer :: mpi_comm_handle = MPI_COMM_WORLD
  real(kind=APP_RKIND) :: value1, value2
  real(kind=APP_RKIND), dimension(:,:), pointer :: MParray
  real(kind=APP_RKIND), dimension(:,:), pointer :: Mesharray
  integer :: ierr, self
  type(c_ptr) :: mpMesh

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle) !this is not supported yet! only for showing
  call polympo_initialize()

  call polympo_checkPrecisionForRealKind(APP_RKIND)
  setMeshOption = 1
  setMPOption = 1    
  mpMesh = polympo_createMPMesh(setMeshOption,setMPOption) !creates test mesh
  
  nverts = 19 !todo use getNumVtx from the Mesh object
  numComps = 2 !todo use getNumComps from velocity fields
  numMPs = 49 !todo use getNumMPs from the MaterialPoints object

  allocate(Mesharray(numComps,nverts))
  allocate(MParray(numComps,numMPs))

  value1 = 42
  MParray = value1
  call polympo_setMPVelArray(mpMesh, numMPs, c_loc(MParray))
  
  MParray = 1
  call polympo_getMPVelArray(mpMesh, numMPs, c_loc(MParray))
  call assert(all(MParray .eq. value1), "Assert MParray == value1 Failed!")

  value2 = 24
  MParray = value2
  call polympo_setMPVelArray(mpMesh, numMPs, c_loc(MParray))

  MParray = 1
  call polympo_getMPVelArray(mpMesh, numMPs, c_loc(MParray))
  call assert(all(MParray .eq. value2), "Assert MParray == value2 Failed!")

  do i = 1,numComps
    do j = 1,nverts 
        Mesharray(i,j) = (i-1)*numComps + j
    end do
  end do
  call polympo_setMeshOnSurfVeloIncrArray(mpMesh, numComps, nverts, c_loc(Mesharray))
  call polympo_setMeshOnSurfDispIncrArray(mpMesh, numComps, nverts, c_loc(Mesharray))

  Mesharray = 1
  call polympo_getMeshOnSurfVeloIncrArray(mpMesh, numComps, nverts, c_loc(Mesharray))
  do i = 1,numComps
    do j = 1,nverts 
        call assert((Mesharray(i,j) .eq. (i-1)*numComps+j), "Assert 2d array Fail")
    end do
  end do
  Mesharray = 1
  call polympo_getMeshOnSurfDispIncrArray(mpMesh, numComps, nverts, c_loc(Mesharray))
  do i = 1,numComps
    do j = 1,nverts 
        call assert((Mesharray(i,j) .eq. (i-1)*numComps+j), "Assert 2d array Fail")
    end do
  end do

  deallocate(MParray)
  deallocate(Mesharray)

  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  stop
end
