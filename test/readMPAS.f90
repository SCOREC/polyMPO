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


subroutine rebuildTests(mpMesh, numMPs, mp2Elm)
    use :: polympo
    use iso_c_binding
    implicit none
    type(c_ptr):: mpMesh
    integer :: numMPs, i, MPACTIVE, MPINACTIVE, MPDELETE_ELM_ID
    integer, dimension(:), pointer :: mp2Elm, addedMPMask, mp2ElmFromPMPO, mp2ElmLarger, addedMPMaskLarger

    MPACTIVE = 1
    MPINACTIVE = 0
    MPDELETE_ELM_ID = -1

    !TODO: introduce MPISACTIVE

    allocate(addedMPMask(numMPs))
    addedMPMask = MPINACTIVE

    ! TEST: adding MP

    ! PREPARE DATA
    do i = 1, numMPs
        if (mp2Elm(i) == MPINACTIVE) then
            addedMPMask(i) = MPACTIVE
            mp2Elm(i) = 1
            exit
        endif
    end do
    ! Rebuild MPs
    call polympo_rebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    ! Test values
    allocate(mp2ElmFromPMPO(numMPs))
    mp2ElmFromPMPO = -1
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    do i = 1, numMPs
        call assert(mp2Elm(i) .eq. mp2ElmFromPMPO(i), "wrong element ID for i'th MP after rebuild")
    end do

    ! TEST: deleting two MPs

    ! PREPARE DATA
    call assert(numMPs > 4, "not enough MPs for test")
    call assert(mp2Elm(1) /= MPDELETE_ELM_ID, "mp2Elm = 1 not active")
    call assert(mp2Elm(4) /= MPDELETE_ELM_ID, "mp2Elm = 4 not active")
    mp2Elm(1) = MPDELETE_ELM_ID
    mp2Elm(4) = MPDELETE_ELM_ID
    addedMPMask = MPINACTIVE
    ! Rebuild MPs
    call polympo_rebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    ! Test values
    mp2ElmFromPMPO = -1
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    do i = 1, numMPs
        if (mp2Elm(i) /= MPDELETE_ELM_ID) then
            call assert(mp2Elm(i) .eq. mp2ElmFromPMPO(i), "wrong element ID for i'th MP after rebuild")
        endif
    end do

    ! TEST: adding 1, replacing 1, removing 1

    ! PREPARE DATA
    addedMPMask(1) = MPACTIVE
    addedMPMask(2) = MPACTIVE
    mp2Elm(1) = 7 !ADDED
    mp2Elm(2) = 7 !REPLACED
    mp2Elm(3) = MPDELETE_ELM_ID !DELETED
    ! Rebuild MPs
    call polympo_rebuildMPs(mpMesh,numMPs,c_loc(mp2Elm),c_loc(addedMPMask))
    ! Test values
    mp2ElmFromPMPO = -1
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2ElmFromPMPO))
    call assert(mp2ElmFromPMPO(1) == 7, "MP = 1 not added")
    call assert(mp2ElmFromPMPO(2) == 7, "MP = 2 not replaced")
    call assert(mp2ElmFromPMPO(3) == MPINACTIVE, "MP = 3 not removed")

    ! TEST: increasing numMP and adding 1 and adding 1 in newMP and deleting 1

    ! PREPARE DATA
    allocate(mp2ElmLarger(numMPs + 10))
    allocate(addedMPMaskLarger(numMPs + 10))

    do i = 1, numMPs
        mp2ElmLarger(i) = mp2Elm(i)
    end do
    addedMPMaskLarger = 0

    mp2ElmLarger(4) = 7
    mp2ElmLarger(5) = MPDELETE_ELM_ID
    mp2ElmLarger(numMPs+8) =  7
    addedMPMaskLarger(4) = MPACTIVE
    addedMPMaskLarger(numMPs+8) =  MPACTIVE
    ! Rebuild MPs
    call polympo_rebuildMPs(mpMesh,numMPs+10,c_loc(mp2ElmLarger),c_loc(addedMPMaskLarger))
    ! Test values
    mp2ElmLarger = -1
    call polympo_getMPCurElmID(mpMesh,numMPs+10,c_loc(mp2ElmLarger))
    call assert(mp2ElmLarger(4) == 7, "MP = 4 not added")
    call assert(mp2ElmLarger(5) == MPINACTIVE, "MP = 5 not deleted")
    call assert(mp2ElmLarger(numMPs+8) == 7, "MP = numMPs+8 not added")
    ! Cleanup
    deallocate(addedMPMask)
    deallocate(mp2ElmFromPMPO)
    deallocate(mp2ElmLarger)
    deallocate(addedMPMaskLarger)
end subroutine

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
    integer, parameter :: nDims = 3
    integer, parameter :: MP_ACTIVE = 1
    integer, parameter :: MP_INACTIVE = 0
    integer :: numMPs
    real(kind=MPAS_RKIND) :: ptOne = 0.100000000000000000
    real(kind=MPAS_RKIND) :: sphereRadius
    integer, dimension(:), pointer :: nEdgesOnCell
    real(kind=MPAS_RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
    integer, dimension(:), pointer :: mpsPerElm, mp2Elm, isMPActive
    real(kind=MPAS_RKIND), dimension(:,:), pointer :: mpPosition

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
    
    isMPActive = MP_ACTIVE !no inactive MPs and some changed below
    isMPActive(4) = MP_INACTIVE !first/1-st MP is indexed 1 and 4-th MP is inactive
   
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
      mpPosition(1,i) = i+ptOne
      mpPosition(2,i) = numMPs+i+ptOne
      mpPosition(3,i) = (2*numMPs)+i+ptOne
    end do

    call polympo_setMPPositions(mpMesh,nDims,numMPs,c_loc(mpPosition))
    mpPosition = 42
    call polympo_getMPPositions(mpMesh,nDims,numMPs,c_loc(mpPosition))
    do i = 1,numMPs
      if(isMPActive(i) .eq. MP_ACTIVE) then
        call assert(epsilonDiff(mpPosition(1,i),i+ptOne), "x position of MP does not match")
        call assert(epsilonDiff(mpPosition(2,i),numMPs+i+ptOne), "y position of MP does not match")
        call assert(epsilonDiff(mpPosition(3,i),(2*numMPs)+i+ptOne), "z position of MP does not match")
      endif
    end do
    
    mp2Elm = -99 !override values and then use get function below
    call polympo_getMPCurElmID(mpMesh,numMPs,c_loc(mp2Elm))
    call assert(mp2Elm(1) .eq. 2, "wrong element ID for MP 1")
    call assert(mp2Elm(2) .eq. 3, "wrong element ID for MP 2")
    call assert(mp2Elm(3) .eq. 2, "wrong element ID for MP 3")
    !mp2Elm(4) is not needed/used since 4-th MP is inactive
    do i = 5,numMPs
      call assert(mp2Elm(i) .eq. i-2, "wrong element ID for i'th MP")
    end do
    
    call rebuildTests(mpMesh, numMPs, mp2Elm)
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
