module readMPAS
    use :: polympo
    use iso_c_binding
    implicit none
    integer, parameter :: MPAS_RKIND = selected_real_kind(12)
    
contains
subroutine polympo_setWithMPASMesh(mpMesh, filename)
    use :: netcdf
    use :: iso_c_binding
    implicit none
   
    type(c_ptr) :: mpMesh
    character (len=*), intent(in) :: filename
    character (len=64) :: onSphere, stringYes = "YES"
    integer :: maxEdges, vertexDegree, nCells, nVertices
    real(kind=MPAS_RKIND) :: sphereRadius
    integer, dimension(:), pointer :: nEdgesOnCell
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
    
    call polympo_readMPASMesh(trim(filename), maxEdges, vertexDegree, &
                              nCells, nVertices, nEdgesOnCell, &
                              onSphere, sphereRadius, &
                              xVertex, yVertex, zVertex, &
                              verticesOnCell, cellsOnCell)
    
    call polympo_checkPrecisionForRealKind(MPAS_RKIND)
    !check on maxEdges and vertexDegree
    call polympo_checkMeshMaxSettings(mpMesh,maxEdges,vertexDegree)

    !set MeshType GeomType sphereRadius
    call polympo_setMeshTypeCVTPoly(mpMesh);
    if (onSphere == stringYes) then
        call polympo_setMeshGeomTypeSpherical(mpMesh);
    else
        call polympo_setMeshGeomTypePlanar(mpMesh);
    end if
        call polympo_setMeshSphereRadius(mpMesh,sphereRadius);

    !set nCells nVertices
    call polympo_setMeshNumVtxs(mpMesh,nVertices)
    call polympo_setMeshNumElms(mpMesh,nCells)

    !set VtxCoords and connectivities
    call polympo_setMeshVtxCoords(mpMesh,nVertices,c_loc(xVertex),c_loc(yVertex),c_loc(zVertex))
    call polympo_setMeshElm2VtxConn(mpMesh,maxEdges,nCells,c_loc(verticesOnCell))
    call polympo_setMeshElm2ElmConn(mpMesh,maxEdges,nCells,c_loc(cellsOnCell))
    call polympo_setMeshNumEdgesPerElm(mpMesh,nCells,c_loc(nEdgesOnCell))

    !unloadMPASMesh to deallocated
    deallocate(nEdgesOnCell)
    deallocate(xVertex)
    deallocate(yVertex)
    deallocate(zVertex)
    deallocate(verticesOnCell)
    deallocate(cellsOnCell)
end subroutine

subroutine polympo_readMPASMesh(filename, maxEdges, vertexDegree, &
                                nCells, nVertices, nEdgesOnCell, &
                                onSphere, sphereRadius, &
                                xVertex, yVertex, zVertex, &
                                verticesOnCell, cellsOnCell)
    use :: netcdf
    use :: iso_c_binding
    implicit none
    
    character (len=*), intent(in) :: filename
    character (len=*), intent(inout) :: onSphere
    character (len=64) :: stringYes = "YES"
    integer, intent(inout) :: maxEdges, vertexDegree, &
                                     nCells, nVertices
    real(kind=MPAS_RKIND) :: sphereRadius
    integer, dimension(:), pointer :: nEdgesOnCell
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell

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
   
    status = nf90_get_att(ncid, nf90_global, "on_a_sphere", onSphere)
    if (status /= nf90_noerr) then
        write(0, *) "polympo_readMPASMesh: Error on get attribute 'on_a_sphere'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if
    if (onSphere == stringYes) then
        status = nf90_get_att(ncid, nf90_global, 'sphere_radius', sphereRadius)
        if (status /= nf90_noerr) then
            write(0, *) "polympo_readMPASMesh: Error on get attribute 'sphere_radius'"
            write(0, *) trim(nf90_strerror(status))
            stop
        end if
    else
        sphereRadius = 0.0
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
subroutine polympo_setWithMPASMeshByFortran(mpMesh, fileName, n) bind(C, name="polympo_setWithMPASMeshByFortran")
    use :: netcdf
    use :: iso_c_binding
    implicit none
    type(c_ptr):: mpMesh
    character (kind=c_char), dimension(*), intent(in) :: fileName
    integer, value :: n
    character (len=n), target :: fileNameFortran
 
    fileNameFortran = transfer(fileName(1:n), fileNameFortran) 
    mpMesh = polympo_createMPMesh(0, 0)

    call polympo_setWithMPASMesh(mpMesh, fileNameFortran)
end subroutine
end module readMPAS
