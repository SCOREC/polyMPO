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

function epsilonDiff(a,b) result(isSame)
  implicit none
  real(kind=MPAS_RKIND) :: a,b,delta
  parameter (delta=1.0e-8)
  logical :: isSame
  if (abs(a-b) < delta) then
    isSame = .true.
  else
    isSame = .false.
  endif
end function


subroutine rebuildTests(mpMesh, numMPs, mp2Elm, isMPActive)
    use :: polympo
    use iso_c_binding
    implicit none
    type(c_ptr):: mpMesh
    integer :: numMPs, i, numMPsLarger
    integer, dimension(:), pointer :: mp2Elm, addedMPMask, isMPActive, mp2ElmFromPMPO
    integer, dimension(:), pointer :: mp2ElmLarger, addedMPMaskLarger, mp2ElmFromPMPOLarger, isMPActiveLarger
    integer, parameter :: MP_ACTIVE = 1
    integer, parameter :: MP_INACTIVE = 0
    integer, parameter :: MP_DELETE_ELM_ID = -1

    ! TEST: adding MP
    ! Necessary pre-conditions for test
    call assert(numMPs >= 4, "not enough MPs for test")
    call assert(isMPActive(4) == MP_INACTIVE, "mp2Elm = 4 is active")
    ! PREPARE DATA
    allocate(addedMPMask(numMPs))
    addedMPMask = MP_INACTIVE
    isMPActive(4) = MP_ACTIVE
    mp2Elm(4) = 7
    addedMPMask(4) = MP_ACTIVE
    ! Rebuild MPs
    call polympo_rebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    ! Test values
    allocate(mp2ElmFromPMPO(numMPs))
    mp2ElmFromPMPO = MP_DELETE_ELM_ID
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    do i = 1, numMPs
        if (isMPActive(i) == MP_ACTIVE) then
            call assert(mp2Elm(i) .eq. mp2ElmFromPMPO(i), "wrong element ID for i'th MP after rebuild")
        endif
    end do

    ! TEST: deleting two MPs
    ! Necessary pre-conditions for test
    call assert(numMPs >= 4, "not enough MPs for test")
    call assert(isMPActive(1) == MP_ACTIVE, "mp2Elm = 1 not active")
    call assert(isMPActive(4) == MP_ACTIVE, "mp2Elm = 4 not active")
    ! PREPARE DATA
    isMPActive(1) = MP_INACTIVE
    isMPActive(4) = MP_INACTIVE
    mp2Elm(1) = MP_DELETE_ELM_ID
    mp2Elm(4) = MP_DELETE_ELM_ID
    addedMPMask = MP_INACTIVE
    ! Rebuild MPs
    call polympo_rebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    ! Test values
    mp2ElmFromPMPO = MP_DELETE_ELM_ID
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    do i = 1, numMPs
        if (isMPActive(i) == MP_ACTIVE) then
            call assert(mp2Elm(i) .eq. mp2ElmFromPMPO(i), "wrong element ID for i'th MP after rebuild")
        endif
    end do

    ! TEST: adding 1, delete 1 and add in same index, removing 1
    ! Necessary pre-conditions for test
    call assert(numMPs >= 3, "not enough MPs for test")
    call assert(isMPActive(1) == MP_INACTIVE, "mp2Elm = 1 not active")
    call assert(isMPActive(2) == MP_ACTIVE, "mp2Elm = 2 is active")
    call assert(isMPActive(3) == MP_ACTIVE, "mp2Elm = 3 not active")
    ! PREPARE DATA
    isMPActive(1) = MP_ACTIVE
    isMPActive(2) = MP_ACTIVE
    isMPActive(3) = MP_INACTIVE
    addedMPMask(1) = MP_ACTIVE
    addedMPMask(2) = MP_ACTIVE
    mp2Elm(1) = 7 !ADDED
    mp2Elm(2) = 7 !REPLACED
    mp2Elm(3) = MP_DELETE_ELM_ID !DELETED
    ! Rebuild MPs
    call polympo_rebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    ! Test values
    mp2ElmFromPMPO = MP_DELETE_ELM_ID
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    do i = 1, numMPs
        if (isMPActive(i) == MP_ACTIVE) then
            call assert(mp2Elm(i) .eq. mp2ElmFromPMPO(i), "wrong element ID for i'th MP after rebuild")
        endif
    end do

    ! TEST: increasing numMPs to numMPsLarger
    ! TEST: adding 1 MP (within numMPs), removing 1 MP (within numMPs), adding 1 (after numMPs, but within numMPsLarger)
    ! Necessary pre-conditions for test
    call assert(numMPs >= 5, "not enough MPs for test")
    call assert(isMPActive(4) == MP_INACTIVE, "mp2Elm = 4 not active")
    call assert(isMPActive(5) == MP_ACTIVE, "mp2Elm = 5 is active")
    ! PREPARE DATA
    numMPsLarger = numMPs+10
    allocate(mp2ElmLarger(numMPsLarger))
    allocate(isMPActiveLarger(numMPsLarger))
    allocate(addedMPMaskLarger(numMPsLarger))
    
    mp2ElmLarger = MP_DELETE_ELM_ID
    isMPActiveLarger = MP_INACTIVE
    addedMPMaskLarger = MP_INACTIVE
    do i = 1, numMPs
        mp2ElmLarger(i) = mp2Elm(i)
        isMPActiveLarger(i) = isMPActive(i)
    end do

    isMPActiveLarger(4) = MP_ACTIVE ! within numMPs
    isMPActiveLarger(5) = MP_INACTIVE ! within numMPs
    isMPActiveLarger(numMPsLarger-2) = MP_ACTIVE ! within numMPsLarger
    mp2ElmLarger(4) = 7
    mp2ElmLarger(5) = MP_DELETE_ELM_ID
    mp2ElmLarger(numMPsLarger-2) =  7
    addedMPMaskLarger(4) = MP_ACTIVE
    addedMPMaskLarger(numMPsLarger-2) = MP_ACTIVE
    ! Rebuild MPs
    call polympo_rebuildMPs(mpMesh,numMPsLarger,c_loc(mp2ElmLarger),c_loc(addedMPMaskLarger))
    ! Test values
    allocate(mp2ElmFromPMPOLarger(numMPsLarger))
    mp2ElmFromPMPOLarger = MP_DELETE_ELM_ID
    call polympo_getMPCurElmID(mpMesh,numMPsLarger,c_loc(mp2ElmFromPMPOLarger))
    do i = 1, numMPs
        if (isMPActiveLarger(i) == MP_ACTIVE) then
            call assert(mp2ElmLarger(i) .eq. mp2ElmFromPMPOLarger(i), "wrong element ID for i'th MP after rebuild")
        endif
    end do
    ! Cleanup
    deallocate(addedMPMask)
    deallocate(mp2ElmFromPMPO)
    deallocate(mp2ElmLarger)
    deallocate(isMPActiveLarger)
    deallocate(addedMPMaskLarger)
    deallocate(mp2ElmFromPMPOLarger)
end subroutine

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
end subroutine

subroutine readMPASMeshFromNCFile(filename, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
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
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell

    integer :: ncid, status, nCellsID, nVerticesID, maxEdgesID, vertexDegreeID, &
               nEdgesOnCellID, xVertexID, yVertexID, zVertexID, &
               latVertexID, lonVertexID, &
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
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
 
    fileNameFortran = transfer(fileName(1:n), fileNameFortran) 
    mpMesh = polympo_createMPMesh(0, 0)

    call readMPASMeshFromNCFile(fileNameFortran, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, lonVertex, &
                        verticesOnCell, cellsOnCell)
    call loadMPASMeshInPolyMPO(mpMesh, maxEdges, vertexDegree, &
                        nCells, nVertices, nEdgesOnCell, &
                        onSphere, sphereRadius, &
                        xVertex, yVertex, zVertex, &
                        latVertex, &
                        verticesOnCell, cellsOnCell)
    
    deallocate(nEdgesOnCell)
    deallocate(xVertex)
    deallocate(yVertex)
    deallocate(zVertex)
    deallocate(verticesOnCell)
    deallocate(cellsOnCell)
end subroutine
end module readMPAS
