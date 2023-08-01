program main
  use :: polympo
  use :: iso_c_binding
  implicit none
  include 'mpif.h'

  integer :: ierr, self
  integer :: argc, i, arglen
  integer :: setMeshOption, setMPOption
  integer(c_int) :: mpi_comm_handle = MPI_COMM_WORLD
  character (len=2048) :: filename
  integer(c_int) :: maxEdges, vertexDegree, nCells, nVertices
  real(c_double) :: sphereRadius
  integer(c_int), dimension(:), pointer :: nEdgesOnCell
  real(c_double), dimension(:), pointer :: xVertex, yVertex, zVertex
  integer(c_int), dimension(:,:), pointer :: verticesOnCell, cellsOnCell
  type(c_ptr) :: mpMesh

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setMPICommunicator(mpi_comm_handle)
  call polympo_initialize()

  argc = command_argument_count()
  if(argc == 1) then
    call get_command_argument(1, filename)
  else
    write(0, *) "Usage: ./testFortranReadMPAS <path to the nc file>"
  end if

  setMeshOption = 1 !create a test mesh
  setMPOption = 1   !create a test set of MPs
  mpMesh = polympo_createMPMesh(setMeshOption, setMPOption)

  call polympo_readMPASMesh(trim(filename), maxEdges, vertexDegree, &
                            nCells, nVertices, nEdgesOnCell, sphereRadius, &
                            xVertex, yVertex, zVertex, &
                            verticesOnCell, cellsOnCell)

  !check on maxEdges and vertexDegree
  call polympo_checkMeshMaxSettings(mpMesh,maxEdges,vertexDegree)

  !set MeshType GeomType sphereRadius
  call polympo_setMeshType(mpMesh,1); !-1=unrecognized,0=general,1=CVT
  call polympo_setMeshGeomType(mpMesh,1); !-1=unrecognized,0=planar,1=spherical
  call polympo_setMeshSphereRadius(mpMesh,sphereRadius);

  !set nCells nVertices
  call polympo_setMeshNumVtxs(mpMesh,nVertices)
  call polympo_setMeshNumElms(mpMesh,nCells)

  !set VtxCoords and connectivities
  call polympo_setMeshVtxCoords(mpMesh,nVertices,c_loc(xVertex),c_loc(yVertex),c_loc(zVertex))
  call polympo_setMeshElm2VtxConn(mpMesh,nCells,maxEdges,c_loc(verticesOnCell))
  call polympo_setMeshNumEdgesPerElm(mpMesh,nCells,c_loc(nEdgesOnCell))
  call polympo_setMeshElm2ElmConn(mpMesh,nCells,maxEdges,c_loc(cellsOnCell))

  !todo check the value using get functions. 
  
  !unloadMPASMesh to deallocated
  deallocate(nEdgesOnCell)
  deallocate(xVertex)
  deallocate(yVertex)
  deallocate(zVertex)
  deallocate(verticesOnCell)
  deallocate(cellsOnCell)

  call polympo_deleteMPMesh(mpMesh)
  call polympo_finalize()

  call mpi_finalize(ierr)

  stop

contains
subroutine polympo_readMPASMesh(filename, maxEdges, vertexDegree, &
                                nCells, nVertices, nEdgesOnCell, sphereRadius, &
                                xVertex, yVertex, zVertex, &
                                verticesOnCell, cellsOnCell)
    use :: netcdf
    use :: iso_c_binding
    implicit none
    
    character (len=*), intent(in) :: filename
    integer(c_int), intent(inout) :: maxEdges, vertexDegree, &
                                     nCells, nVertices
    real(c_double) :: sphereRadius
    integer(c_int), dimension(:), pointer :: nEdgesOnCell
    real(c_double), dimension(:), pointer :: xVertex, yVertex, zVertex
    integer(c_int), dimension(:,:), pointer :: verticesOnCell, cellsOnCell

    integer :: ncid, status, nCellsID, nVerticesID, maxEdgesID, vertexDegreeID, &
               nEdgesOnCellID, xVertexID, yVertexID, zVertexID, &
               verticesOnCellID, cellsOnCellID
    
    status = nf90_open(path=trim(filename), mode=nf90_nowrite, ncid=ncid)
    if (status /= nf90_noerr) then
        write(0,*) 'polympo_readMPASMesh: Error occured when opening MPAS grid:'//filename
        write(0,*) trim(nf90_strerror(status))
        stop 
    end if

   status = nf90_inq_dimid(ncid, 'nCells', nCellsID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error when getting dimid of 'nCells'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_dimid(ncid, 'nVertices', nVerticesID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error when getting dimid of 'nVertices'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, nCellsID, len = nCells)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire dimension of 'nCellsID'"
        write(0, *) trim(nf90_strerror(status))
        stop
end if

    status = nf90_inq_dimid(ncid, 'maxEdges', maxEdgesID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error when getting dimid of 'maxEdgesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_dimid(ncid, 'vertexDegree', vertexDegreeID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error when getting dimid of 'vertexDegreeID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, nVerticesID, len = nVertices)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire dimension of 'nVerticesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, maxEdgesID, len = maxEdges)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire dimension of 'maxEdgesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, vertexDegreeID, len = vertexDegree)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire dimension of 'vertexDegreeID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    allocate(xVertex(nVertices))
    allocate(yVertex(nVertices))
    allocate(zVertex(nVertices))
    allocate(nEdgesOnCell(nCells))
    allocate(verticesOnCell(maxEdges, nCells))
    allocate(cellsOnCell(maxEdges, nCells))

    status = nf90_get_att(ncid, nf90_global, 'sphere_radius', sphereRadius)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get attribute 'sphereRadius'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if
    
    status = nf90_inq_varid(ncid, 'xVertex', xVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire varid of 'xVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'yVertex', yVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire varid of 'yVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'zVertex', zVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire varid of 'zVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'nEdgesOnCell', nEdgesOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire varid of 'nEdgesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'verticesOnCell', verticesOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire varid of 'verticesOnCellID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'cellsOnCell', cellsOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire varid of 'cellsOnCellID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, xVertexID, xVertex)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get var of 'xVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if
    status = nf90_get_var(ncid, yVertexID, yVertex)

    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get var of 'yVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, zVertexID, zVertex)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get var of 'zVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, nEdgesOnCellID, nEdgesOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get var of 'nEdgesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, verticesOnCellID, verticesOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get var of 'verticesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, cellsOnCellID, cellsOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get var of 'cellsOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_close(ncid)   
end subroutine polympo_readMPASMesh
end
