subroutine assert(condition,message)
  implicit none
  logical :: condition
  character(*) :: message
  if (condition .neqv. .true.) then
    write (*,*) message
    call exit(1)
  endif
end subroutine

module readMPAS
    use :: polympo
    use iso_c_binding
    implicit none
    integer, parameter :: MPAS_RKIND = selected_real_kind(12)
    
contains

!---------------------------------------------------------------------------
!> @brief get the MP positions array from a polympo array
!> @param mpmesh(in/out) MPMesh object to fill, allocated by users
!> @param filename(in) the .nc file want to read
!---------------------------------------------------------------------------
subroutine loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, &
                        xCell, yCell, zCell, &
                        verticesOnCell, cellsOnCell)
    use :: netcdf
    use :: iso_c_binding
    implicit none
   
    type(c_ptr), intent(inout) :: mpMesh
    character (len=64), intent(in) :: onSphere
    character (len=64) :: stringYes = "YES"
    integer, intent(in) :: maxEdges, vertexDegree, nCells, nVertices
    real(kind=MPAS_RKIND), intent(in) :: sphereRadius
    integer, dimension(:), pointer :: nEdgesOnCell
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
    real(kind=MPAS_RKIND), dimension(:), pointer :: xCell, yCell, zCell
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell

    call polympo_checkPrecisionForRealKind(MPAS_RKIND)
    !check on maxEdges and vertexDegree
    call polympo_checkMeshMaxSettings(mpMesh,maxEdges,vertexDegree)

    call polympo_startMeshFill(mpMesh)
    !set MeshType GeomType sphereRadius
    call polympo_setMeshTypeCVTPoly(mpMesh)
    if (onSphere == stringYes) then
        call polympo_setMeshGeomTypeSpherical(mpMesh)
    else
        call polympo_setMeshGeomTypePlanar(mpMesh)
    end if
        call polympo_setMeshSphereRadius(mpMesh,sphereRadius)

    !set nCells nVertices
    call polympo_setMeshNumVtxs(mpMesh,nVertices)
    call polympo_setMeshNumElms(mpMesh,nCells)

    !set connectivities
    call polympo_setMeshElm2VtxConn(mpMesh,maxEdges,nCells,c_loc(verticesOnCell))
    call polympo_setMeshElm2ElmConn(mpMesh,maxEdges,nCells,c_loc(cellsOnCell))
    call polympo_setMeshNumEdgesPerElm(mpMesh,nCells,c_loc(nEdgesOnCell))
    
    !end mesh structure fill
    call polympo_endMeshFill(mpMesh)

    !set vtxCoords which is a mesh field 
    call polympo_setMeshVtxCoords(mpMesh,nVertices,c_loc(xVertex),c_loc(yVertex),c_loc(zVertex))
    call polympo_setMeshVtxRotLat(mpMesh,nVertices,c_loc(latVertex))
    
    !set mesh element center
    call polympo_setMeshElmCenter(mpMesh,nCells,c_loc(xCell),c_loc(yCell),c_loc(zCell))
end subroutine

subroutine readMPASMeshFromNCFile(filename, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
                        xCell, yCell, zCell, &
                        verticesOnCell, cellsOnCell)
    use :: netcdf
    use :: iso_c_binding
    implicit none
    
    character (len=*), intent(inout) :: filename
    character (len=*), intent(inout) :: onSphere
    character (len=64) :: stringYes = "YES"
    integer, intent(inout) :: maxEdges, vertexDegree, &
                                     nCells, nVertices
    real(kind=MPAS_RKIND) :: sphereRadius
    integer, dimension(:), pointer :: nEdgesOnCell
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
    real(kind=MPAS_RKIND), dimension(:), pointer :: xCell, yCell, zCell
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell

    integer :: ncid, status, nCellsID, nVerticesID, maxEdgesID, vertexDegreeID, &
               nEdgesOnCellID, xVertexID, yVertexID, zVertexID, &
               latVertexID, lonVertexID, &
               xCellID, yCellID, zCellID, &
               verticesOnCellID, cellsOnCellID
    
    status = nf90_open(path=trim(filename), mode=nf90_nowrite, ncid=ncid)
    if (status /= nf90_noerr) then
        write(0,*) 'readMPASMeshFromNCFile: Error occured when opening MPAS grid:'//filename
        write(0,*) trim(nf90_strerror(status))
        stop 
    end if

   status = nf90_inq_dimid(ncid, 'nCells', nCellsID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error when getting dimid of 'nCells'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_dimid(ncid, 'nVertices', nVerticesID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error when getting dimid of 'nVertices'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, nCellsID, len = nCells)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire dimension of 'nCellsID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_dimid(ncid, 'maxEdges', maxEdgesID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error when getting dimid of 'maxEdgesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_dimid(ncid, 'vertexDegree', vertexDegreeID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error when getting dimid of 'vertexDegreeID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, nVerticesID, len = nVertices)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire dimension of 'nVerticesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, maxEdgesID, len = maxEdges)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire dimension of 'maxEdgesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, vertexDegreeID, len = vertexDegree)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire dimension of 'vertexDegreeID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    allocate(xVertex(nVertices))
    allocate(yVertex(nVertices))
    allocate(zVertex(nVertices))
    allocate(latVertex(nVertices))
    allocate(lonVertex(nVertices))
    allocate(xCell(nCells))
    allocate(yCell(nCells))
    allocate(zCell(nCells))
    allocate(nEdgesOnCell(nCells))
    allocate(verticesOnCell(maxEdges, nCells))
    allocate(cellsOnCell(maxEdges, nCells))
   
    status = nf90_get_att(ncid, nf90_global, "on_a_sphere", onSphere)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get attribute 'on_a_sphere'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if
    if (onSphere == stringYes) then
        status = nf90_get_att(ncid, nf90_global, 'sphere_radius', sphereRadius)
        if (status /= nf90_noerr) then
            write(0, *) "readMPASMeshFromNCFile: Error on get attribute 'sphere_radius'"
            write(0, *) trim(nf90_strerror(status))
            stop
        end if
    else
        sphereRadius = 0.0
    end if
    
    status = nf90_inq_varid(ncid, 'xVertex', xVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'xVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'yVertex', yVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'yVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'zVertex', zVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'zVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'latVertex', latVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'latVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'lonVertex', lonVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'lonVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'xCell', xCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'xCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'yCell', yCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'yCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'zCell', zCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'zCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'nEdgesOnCell', nEdgesOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'nEdgesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'verticesOnCell', verticesOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'verticesOnCellID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'cellsOnCell', cellsOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on inquire varid of 'cellsOnCellID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, xVertexID, xVertex)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'xVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, yVertexID, yVertex)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'yVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, zVertexID, zVertex)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'zVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, latVertexID, latVertex)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'latVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, lonVertexID, lonVertex)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'lonVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, xCellID, xCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'xCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, yCellID, yCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'yCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, zCellID, zCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'zCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, nEdgesOnCellID, nEdgesOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'nEdgesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, verticesOnCellID, verticesOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'verticesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if
    
    status = nf90_get_var(ncid, cellsOnCellID, cellsOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMeshFromNCFile: Error on get var of 'cellsOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_close(ncid)   
end subroutine readMPASMeshFromNCFile
subroutine setWithMPASMeshByFortran(mpMesh, fileName, n) bind(C, name="setWithMPASMeshByFortran")
    use :: netcdf
    use :: iso_c_binding
    implicit none
    type(c_ptr):: mpMesh
    character (kind=c_char), dimension(*), intent(in) :: fileName
    integer :: maxEdges, vertexDegree, nCells, nVertices
    integer, value :: n
    character (len=n), target :: fileNameFortran
    character (len=64) :: onSphere, stringYes = "YES"
    real(kind=MPAS_RKIND) :: sphereRadius
    integer, dimension(:), pointer :: nEdgesOnCell
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
    real(kind=MPAS_RKIND), dimension(:), pointer :: xCell, yCell, zCell
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
 
    fileNameFortran = transfer(fileName(1:n), fileNameFortran) 
    mpMesh = polympo_createMPMesh(0, 0)

    call readMPASMeshFromNCFile(fileNameFortran, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
                        xCell, yCell, zCell, &
                        verticesOnCell, cellsOnCell)
    call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, &
                        xCell, yCell, zCell, &
                        verticesOnCell, cellsOnCell)
    
    deallocate(nEdgesOnCell)
    deallocate(xVertex)
    deallocate(yVertex)
    deallocate(zVertex)
    deallocate(xCell)
    deallocate(yCell)
    deallocate(zCell)
    deallocate(verticesOnCell)
    deallocate(cellsOnCell)
end subroutine
end module readMPAS
