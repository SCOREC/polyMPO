program main
  use :: polympo
  use :: iso_c_binding
  use :: mpi_f08
  implicit none

  integer :: ierr, self
  TYPE(MPI_Comm) :: mpi_comm_handle = MPI_COMM_WORLD
  character (len=128) :: filename = "./grid_full.nc"
  integer(c_int) :: maxEdges, vertexDegree, nCells, nVertices
  integer(c_int), dimension(:), pointer :: nEdgesOnCell
  real(c_double), dimension(:), pointer :: xVertex, yVertex, zVertex
  integer(c_int), dimension(:,:), pointer :: verticesOnCell, cellsOnVertex, cellsOnCell
  type(c_ptr) :: mpMesh

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_handle, self, ierr)

  call polympo_setCommunicator(mpi_comm_handle%MPI_VAL)
  call polympo_initialize()

  mpMesh = polympo_createMpMesh() !creates test mesh
  call polympo_readMPASMesh(filename, maxEdges, vertexDegree, &
                            nCells, nVertices, nEdgesOnCell, &
                            xVertex, yVertex, zVertex, &
                            verticesOnCell, cellsOnVertex, cellsOnCell)

  !check on maxEdges and vertexDegree
  call polympo_checkMeshSetting(mpMesh,maxEdges,vertexDegree)
  !set nCells nVertices
  call polympo_setNumVtxs(mpMesh,nVertices)
  call polympo_setNumElms(mpMesh,nCells)
  !todo 1d array 2d array
  call polympo_setMeshVtxCoords(mpMesh,nVertices,c_loc(xVertex),c_loc(yVertex),c_loc(zVertex))
  
  !call polympo_setMeshElm2VtxConn(mpMesh,nCells,maxEdges,c_loc(verticesOnCell))
  !call polympo_setMeshVtx2ElmConn(mpMesh,nVertices,vertexDegree,c_loc(cellsOnVertex))
  call polympo_setMeshElm2ElmConn(mpMesh,nCells,maxEdges,c_loc(cellsOnCell))

  !todo how to check the value 
  call polympo_printMesh(mpMesh)
  !unloadMPASMesh to deallocated
  deallocate(nEdgesOnCell)
  deallocate(xVertex)
  deallocate(yVertex)
  deallocate(zVertex)
  deallocate(verticesOnCell)
  deallocate(cellsOnVertex)
  deallocate(cellsOnCell)
  call polympo_finalize()
  call mpi_finalize(ierr)
  stop

contains
subroutine polympo_readMPASMesh(filename, maxEdges, vertexDegree, &
                                nCells, nVertices, nEdgesOnCell, &
                                xVertex, yVertex, zVertex, &
                                verticesOnCell, cellsOnVertex, cellsOnCell)
    use :: netcdf
    use :: iso_c_binding
    implicit none
    
    character (len=*), intent(in) :: filename
    integer(c_int), intent(inout) :: maxEdges, vertexDegree, nCells, nVertices
    integer(c_int), dimension(:), pointer :: nEdgesOnCell
    real(c_double), dimension(:), pointer :: xVertex, yVertex, zVertex
    integer(c_int), dimension(:,:), pointer :: verticesOnCell, cellsOnVertex, cellsOnCell

    integer :: ncid, status, nCellsID, nVerticesID, maxEdgesID, vertexDegreeID, &
               nEdgesOnCellID, xVertexID, yVertexID, zVertexID, &
               verticesOnCellID, cellsOnVertexID, cellsOnCellID
    
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

    status = nf90_inq_dimid(ncid, 'maxEdges', maxEdgesID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error when getting dimid of 'maxEdges'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, nCellsID, len = nCells)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire dimension of 'nCellsID'"
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

    allocate(xVertex(nVertices))
    allocate(yVertex(nVertices))
    allocate(zVertex(nVertices))
    allocate(nEdgesOnCell(nCells))
    allocate(verticesOnCell(maxEdges, nCells))
    allocate(cellsOnVertex(maxEdges, nCells))
    allocate(cellsOnCell(maxEdges, nCells))

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

    status = nf90_inq_varid(ncid, 'verticesOnCell', cellsOnVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on inquire varid of 'cellsOnVertexID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'verticesOnCell', cellsOnCellID)
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

    status = nf90_get_var(ncid, cellsOnVertexID, cellsOnVertex)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get var of 'cellsOnVertex'"
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
