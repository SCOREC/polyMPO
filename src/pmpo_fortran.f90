!---------------------------------------------------------------------------
!> @file pmo_fortran.f90
!> @brief PolyMPO FORTRAN interface using iso_c_binding
!---------------------------------------------------------------------------
module polympo
  use :: iso_c_binding
  public
  interface
  !---------------------------------------------------------------------------
  !> @brief initialize polympo, call this before any other polympo api
  !> @remark the user must initialize MPI prior to this call
  !---------------------------------------------------------------------------
  subroutine polympo_initialize() bind(C, NAME='polympo_initialize')
    use :: iso_c_binding
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief finalize polympo, no polympo apis may be called after this
  !> @remark the user must not finalize MPI until after this call
  !---------------------------------------------------------------------------
  subroutine polympo_finalize() bind(C, NAME='polympo_finalize')
    use :: iso_c_binding
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief create MPMesh object
  !> @param testMeshOption >0 init a planar test Mesh
  !>                       <0 init a blank Mesh
  !> @param testMPOption   >0 init a set of test MPs
  !>                       <0 init a blank set of MPs
  !> @return polympo_createMPMesh(out) MPMesh object
  !---------------------------------------------------------------------------
  function polympo_createMPMesh(setMeshOption, setMPOption) bind(C, NAME='polympo_createMPMesh')
    use :: iso_c_binding
    type(c_ptr) polympo_createMPMesh
    integer(c_int), value :: setMeshOption, setMPOption
  end function
  !---------------------------------------------------------------------------
  !> @brief delete MPMesh object
  !> @param mpmesh(in/out) MPMesh object
  !---------------------------------------------------------------------------
  subroutine polympo_deleteMPMesh(mpMesh) bind(C, NAME='polympo_deleteMPMesh')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the MPI communicator used by polympo
  !> @param comm(in) MPI communicator
  !---------------------------------------------------------------------------
  subroutine polympo_setMPICommunicator(comm) &
             bind(C, NAME='polympo_setMPICommunicator')
    use :: iso_c_binding
    integer(c_int), value :: comm    
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief create the material points
  !> @brief the fields associated with the MPs are NOT initialized
  !> @param mpmesh(in/out) MPMesh object
  !> @param numElms(in) total number of mesh elements
  !> @param numMPs(in) total number of MPs, total = number of active + number of inactive
  !> @param mpsPerElm(in) number of MPs for each mesh element
  !> @param mp2Elm(in) element ID for each MP
  !> @param isMPActive(in) set to 1 if the MP is active, 0 otherwise
  !---------------------------------------------------------------------------
  subroutine polympo_createMPs(mpMesh, numElms, numMPs, mpsPerElm, mp2Elm, isMPActive) &
             bind(C, NAME='polympo_createMPs')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numElms
    integer(c_int), value :: numMPs
    type(c_ptr), intent(in), value :: mpsPerElm
    type(c_ptr), intent(in), value :: mp2Elm
    type(c_ptr), intent(in), value :: isMPActive
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the current element ID MP array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param numMPs(in) length of array, number of the MPs
  !> @param array(in/out) output MP element ID 1D array (numMPs), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPCurElmID(mpMesh, numMPs, array) &
             bind(C, NAME='polympo_getMPCurElmID')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numMPs
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the MP positions array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 3
  !> @param numMPs(in) number of the MPs
  !> @param array(in/out) output MP current position 2D array (3,numMPs),
  !>                      allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPPositions(mpMesh, nComps, numMPs, array) &
             bind(C, NAME='polympo_getMPPositions')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the velocity MP array from a host array
  !> @warning THIS IS NOT SUPPORTED YET 
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) half length of array
  !> @param array(in) input MP velocity 1D array (numMPs*2)
  !---------------------------------------------------------------------------
  subroutine polympo_setMPVel(mpMesh, n, array) &
             bind(C, NAME='polympo_setMPVel')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the velocity MP array from a polympo array
  !> @warning THIS IS NOT SUPPORTED YET 
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) half length of array
  !> @param array(in/out) output MP velocity 1D array (numMPs*2), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPVel(mpMesh, n, array) &
             bind(C, NAME='polympo_getMPVel')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief Enable the setting of mesh topology (number of entities and entity adjacencies). 
  !>        By default, modifying the mesh topology without calling this function first will result
  !>        in a runtime failure.
  !> @param mpMesh(in/out) the MPMesh is valid/created
  !---------------------------------------------------------------------------
  subroutine polympo_startMeshFill(mpMesh) &
             bind(C, NAME='polympo_startMeshFill')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief Disable the modification of mesh topology (number of entities and entity adjacencies). 
  !>        Once the topology is set, call this function to prevent unexpected/accidental modifications.
  !> @param mpMesh(in/out) the MPMesh is valid/created
  !---------------------------------------------------------------------------
  subroutine polympo_endMeshFill(mpMesh) &
             bind(C, NAME='polympo_endMeshFill')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief check the Mesh is valid/runable in polympo
  !> @param mpMesh(in/out) the MPMesh is valid/created
  !> @param maxEdges(in) the maxEdges per cell of your mesh 
  !> @param vertexDegree(in) the max vertexDegree of a vertex 
  !---------------------------------------------------------------------------
  subroutine polympo_checkMeshMaxSettings(mpMesh,maxEdges,vertexDegree) &
             bind(C, NAME='polympo_checkMeshMaxSettings')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: maxEdges, vertexDegree
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh type to general polygonal
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshTypeGeneralPoly(mpMesh) &
             bind(C, NAME='polympo_setMeshTypeGeneralPoly')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh type to CVT polygonal
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshTypeCVTPoly(mpMesh) &
             bind(C, NAME='polympo_setMeshTypeCVTPoly')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh geometry type to planar
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshGeomTypePlanar(mpMesh) &
             bind(C, NAME='polympo_setMeshGeomTypePlanar')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh geometry type to spherical
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshGeomTypeSpherical(mpMesh) &
             bind(C, NAME='polympo_setMeshGeomTypeSpherical')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh sphere radius
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshSphereRadius(mpMesh,sphereRadius) &
             bind(C, NAME='polympo_setMeshSphereRadius')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    real(c_double), value :: sphereRadius
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the number of vetices of the mesh
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !> @param numVtxs(in) the number of vertices need to set
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshNumVtxs(mpMesh,numVtxs) &
             bind(C, NAME='polympo_setMeshNumVtxs')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numVtxs
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the number of vertices from the mesh holding by polyMPO
  !> @param mpMesh(in) mpMesh object
  !> @param numVtxs(return)) the number of vertices
  !---------------------------------------------------------------------------
  function polympo_getMeshNumVtxs(mpMesh) result(numVtxs) &
            bind(C, NAME = 'polympo_getMeshNumVtxs')
    use :: iso_c_binding
    type(c_ptr), intent(in), value :: mpMesh
    integer(c_int) :: numVtxs
  end function
  !---------------------------------------------------------------------------
  !> @brief set the number of elements of the mesh
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !> @param numElms(in) the number of elements
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshNumElms(mpMesh,numElms) &
             bind(C, NAME='polympo_setMeshNumElms')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numElms
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the number of elements from the mesh holding by polyMPO
  !> @param mpMesh(in) mpMesh object
  !> @param numVtxs(return)) the number of elements
  !---------------------------------------------------------------------------
  function polympo_getMeshNumElms(mpMesh) result(numElms) &
            bind(C, NAME = 'polympo_getMeshNumElms')
    use :: iso_c_binding
    type(c_ptr), intent(in), value :: mpMesh
    integer(c_int) :: numElms
  end function
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh element to vertices connectivity
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpmesh(in/out) MPMesh object
  !> @param maxEdges,nCells(in) length of array in each direction
  !> @param verticesOnCell(in) element to vertices connectivity 2D array 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshElm2VtxConn(mpMesh, maxEdges, nCells, verticesOnCell) &
             bind(C, NAME='polympo_setMeshElm2VtxConn')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: maxEdges, nCells
    type(c_ptr), intent(in), value :: verticesOnCell
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh element to elements connectivity
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpmesh(in/out) MPMesh object
  !> @param maxEdges,nCells(in) length of array in each direction 
  !> @param cellsOnCell(in) element to elements connectivity 2D array 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshElm2ElmConn(mpMesh, maxEdges, nCells, cellsOnCell) &
             bind(C, NAME='polympo_setMeshElm2ElmConn')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: maxEdges, nCells
    type(c_ptr), intent(in), value :: cellsOnCell
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh number of edges per element
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpmesh(in/out) MPMesh object
  !> @param nCells(in) length of array (numElms)
  !> @param nEdgesOnCell(in) number of edges per element
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshNumEdgesPerElm(mpMesh, nCells, nEdgesOnCell) &
             bind(C, NAME='polympo_setMeshNumEdgesPerElm')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nCells
    type(c_ptr), intent(in), value :: nEdgesOnCell
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh vertices coordinates
  !> @param mpmesh(in/out) MPMesh object
  !> @param nVertices(in) length of array in 
  !> @param x/y/zArray(in) the 1D arrays of vertices coordinates
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshVtxCoords(mpMesh, nVertices, xArray, yArray, zArray) &
             bind(C, NAME='polympo_setMeshVtxCoords')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nVertices
    type(c_ptr), intent(in), value :: xArray, yArray, zArray
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the polympo mesh vertices coordinates
  !> @param mpmesh(in/out) MPMesh object
  !> @param nVertices(in) length of array in, use for assertion
  !> @param x/y/zArray(in/out) the 1D arrays of vertices coordinates
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshVtxCoords(mpMesh, nVertices, xArray, yArray, zArray) &
             bind(C, NAME='polympo_getMeshVtxCoords')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nVertices
    type(c_ptr), value :: xArray, yArray, zArray
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh vertices latitude and longitude
  !> @param mpmesh(in/out) MPMesh object
  !> @param nVertices(in) length of array in 
  !> @param latitude/longitude(in) the 1D arrays of vertices lat/lon
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshVtxLatLon(mpMesh, nVertices, latitude, longitude) &
             bind(C, NAME='polympo_setMeshVtxLatLon')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nVertices
    type(c_ptr), intent(in), value :: latitude, longitude
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the polympo mesh vertices latitude and longitude
  !> @param mpmesh(in/out) MPMesh object
  !> @param nVertices(in) length of array in, use for assertion
  !> @param latitude/longitude(in/out) the 1D arrays of vertices lat/lon
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshLatLon(mpMesh, nVertices, latitude, longitude) &
             bind(C, NAME='polympo_getMeshLatLon')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nVertices
    type(c_ptr), value :: latitude, longitude
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the spherical velocity increment mesh array 
  !>        from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 2
  !> @param nVertices(in) numVertices
  !> @param array(in) input mesh velocity 2D array (2,numVtx)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshOnSurfVeloIncr(mpMesh, nComps, nVertices, array) &
             bind(C, NAME='polympo_setMeshOnSurfVeloIncr')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, nVertices
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the spherical velocity increment mesh array 
  !>        from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 2
  !> @param nVertices(in) numVertices
  !> @param array(in/out) output mesh spherical velocity increment
  !>        2D array (2,numVtx), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshOnSurfVeloIncr(mpMesh, nComps, nVertices, array) &
             bind(C, NAME='polympo_getMeshOnSurfVeloIncr')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, nVertices
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the spherical displacement increment mesh array 
  !>        from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 2
  !> @param nVertices(in) numVertices
  !> @param array(in) input mesh velocity 2D array (2,numVtx)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshOnSurfDispIncr(mpMesh, nComps, nVertices, array) &
             bind(C, NAME='polympo_setMeshOnSurfDispIncr')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, nVertices
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the spherical displacement increment mesh array
  !>        from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 2
  !> @param nVertices(in) numVertices
  !> @param array(in/out) output mesh spherical displacement increment 
  !>        2D array (2,numVtx), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshOnSurfDispIncr(mpMesh, nComps, nVertices, array) &
             bind(C, NAME='polympo_getMeshOnSurfDispIncr')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, nVertices
    type(c_ptr), value :: array
  end subroutine
  end interface
  contains
  !---------------------------------------------------------------------------
  !> @brief check to make sure the real value is same precision with c_double
  !> @param APP_RKIND(in) the real value kind use in fortran
  !---------------------------------------------------------------------------
  subroutine polympo_checkPrecisionForRealKind(APP_RKIND)
    use :: iso_c_binding
    implicit none
    integer, intent(in) :: APP_RKIND
    if (APP_RKIND .ne. c_double) then
        write(0, *) "polyMPO: Precision does not match!"
        call exit(1)
    end if
  end subroutine
end module
