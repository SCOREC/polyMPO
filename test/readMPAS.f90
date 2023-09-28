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
subroutine loadMPASMesh(mpMesh, filename)
    use :: netcdf
    use :: iso_c_binding
    implicit none
   
    type(c_ptr) :: mpMesh
    character (len=*), intent(in) :: filename
    character (len=64) :: onSphere, stringYes = "YES"
    integer :: i
    integer :: maxEdges, vertexDegree, nCells, nVertices
    integer :: nDims !FIXME - parameter???
    integer :: numMPs
    real(kind=MPAS_RKIND) :: sphereRadius
    integer, dimension(:), pointer :: nEdgesOnCell
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
    integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive
    real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition

    nDims = 3
    
    call readMPASMesh(trim(filename), maxEdges, vertexDegree, &
                              nCells, nVertices, nEdgesOnCell, &
                              onSphere, sphereRadius, &
                              xVertex, yVertex, zVertex, &
                              verticesOnCell, cellsOnCell)
    
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

    !test on new createMPs
    call assert(nCells .ge. 3, "This test requires a mesh with at least three cells")
    numMPs = nCells+2;
    allocate(mpsPerElm(nCells))
    allocate(mp2Elm(numMPs))
    allocate(isMPActive(numMPs))
    
    isMPActive = 1 !no inactive MPs and some changed below
    isMPActive(4) = 0 !first/1-st MP is indexed 1 and 4-th MP is inactive
   
    mpsPerElm = 1 !all elements have 1 MP and some changed below
    mpsPerElm(1) = 0 !1st element has 0 MPs
    mpsPerElm(2) = 2 !2nd element has 2 MPs
    mpsPerElm(3) = 2 !3rd element has 2 MPs 

    ! mp2Elm = [2,3,2,0,3,4,5,6,...]
    mp2Elm(1) = 2
    mp2Elm(2) = 3
    mp2Elm(3) = 2
    !mp2Elm(4) is not needed/used since 4-th MP is inactive
    do i = 5,numMPs
      mp2Elm(i) = i-2 !i=5 leads to mp2Elm(5)=3 (5-th MP in 3-rd element)
                      !i=numMPs leads to mp2Elm(numMPs=nCells+2)=numMPs-2=nCells
    end do
    call polympo_createMPs(mpMesh,nCells,numMPs,c_loc(mpsPerElm),c_loc(mp2Elm),c_loc(isMPActive))

    !set mp positions
    allocate(mpPosition(nDims,numMPs))
    do i = 1,numMPs
      mpPosition(1,i) = i+0.1
      mpPosition(2,i) = numMPs+i+0.1
      mpPosition(3,i) = (2*numMPs)+i+0.1
    end do

    call polympo_setMPPositions(mpMesh,nDims,numMPs,c_loc(mpPosition))
    
    mp2Elm = -99 !override values and then use get function below
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2Elm))
    call assert(mp2Elm(1) .eq. 2, "wrong element ID for MP 1")
    call assert(mp2Elm(2) .eq. 3, "wrong element ID for MP 2")
    call assert(mp2Elm(3) .eq. 2, "wrong element ID for MP 3")
    !mp2Elm(4) is not needed/used since 4-th MP is inactive
    do i = 5,numMPs
      call assert(mp2Elm(i) .eq. i-2, "wrong element ID for i'th MP")
    end do
    !test end

    deallocate(mpPosition)
    deallocate(mpsPerElm)
    deallocate(mp2Elm)
    deallocate(isMPActive)

    !set vtxCoords which is a mesh field 
    call polympo_setMeshVtxCoords(mpMesh,nVertices,c_loc(xVertex),c_loc(yVertex),c_loc(zVertex))
    !unloadMPASMesh to deallocated
    deallocate(nEdgesOnCell)
    deallocate(xVertex)
    deallocate(yVertex)
    deallocate(zVertex)
    deallocate(verticesOnCell)
    deallocate(cellsOnCell)
end subroutine

subroutine readMPASMesh(filename, maxEdges, vertexDegree, &
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
        write(0,*) 'readMPASMesh: Error occured when opening MPAS grid:'//filename
        write(0,*) trim(nf90_strerror(status))
        stop 
    end if

   status = nf90_inq_dimid(ncid, 'nCells', nCellsID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error when getting dimid of 'nCells'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_dimid(ncid, 'nVertices', nVerticesID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error when getting dimid of 'nVertices'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, nCellsID, len = nCells)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire dimension of 'nCellsID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_dimid(ncid, 'maxEdges', maxEdgesID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error when getting dimid of 'maxEdgesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_dimid(ncid, 'vertexDegree', vertexDegreeID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error when getting dimid of 'vertexDegreeID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, nVerticesID, len = nVertices)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire dimension of 'nVerticesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, maxEdgesID, len = maxEdges)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire dimension of 'maxEdgesID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inquire_dimension(ncid, vertexDegreeID, len = vertexDegree)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire dimension of 'vertexDegreeID'"
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
        write(0, *) "readMPASMesh: Error on get attribute 'on_a_sphere'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if
    if (onSphere == stringYes) then
        status = nf90_get_att(ncid, nf90_global, 'sphere_radius', sphereRadius)
        if (status /= nf90_noerr) then
            write(0, *) "readMPASMesh: Error on get attribute 'sphere_radius'"
            write(0, *) trim(nf90_strerror(status))
            stop
        end if
    else
        sphereRadius = 0.0
    end if
    
    status = nf90_inq_varid(ncid, 'xVertex', xVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire varid of 'xVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'yVertex', yVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire varid of 'yVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'zVertex', zVertexID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire varid of 'zVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'nEdgesOnCell', nEdgesOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire varid of 'nEdgesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'verticesOnCell', verticesOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire varid of 'verticesOnCellID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_inq_varid(ncid, 'cellsOnCell', cellsOnCellID)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on inquire varid of 'cellsOnCellID'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, xVertexID, xVertex)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on get var of 'xVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if
    status = nf90_get_var(ncid, yVertexID, yVertex)

    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on get var of 'yVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, zVertexID, zVertex)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on get var of 'zVertex'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, nEdgesOnCellID, nEdgesOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on get var of 'nEdgesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_get_var(ncid, verticesOnCellID, verticesOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on get var of 'verticesOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if
    
    status = nf90_get_var(ncid, cellsOnCellID, cellsOnCell)
    if (status /= nf90_noerr) then
        write(0, *) "readMPASMesh: Error on get var of 'cellsOnCell'"
        write(0, *) trim(nf90_strerror(status))
        stop
    end if

    status = nf90_close(ncid)   
end subroutine readMPASMesh
subroutine setWithMPASMeshByFortran(mpMesh, fileName, n) bind(C, name="setWithMPASMeshByFortran")
    use :: netcdf
    use :: iso_c_binding
    implicit none
    type(c_ptr):: mpMesh
    character (kind=c_char), dimension(*), intent(in) :: fileName
    integer, value :: n
    character (len=n), target :: fileNameFortran
 
    fileNameFortran = transfer(fileName(1:n), fileNameFortran) 
    mpMesh = polympo_createMPMesh(0, 0)

    call loadMPASMesh(mpMesh, fileNameFortran)
end subroutine
end module readMPAS
