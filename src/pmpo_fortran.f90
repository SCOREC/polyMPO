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
  subroutine polympo_initialize() bind(C, NAME='polympo_initialize_f')
    use :: iso_c_binding
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief finalize polympo, no polympo apis may be called after this
  !> @remark the user must not finalize MPI until after this call
  !---------------------------------------------------------------------------
  subroutine polympo_finalize() bind(C, NAME='polympo_finalize_f')
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
  function polympo_createMPMesh(setMeshOption, setMPOption) bind(C, NAME='polympo_createMPMesh_f')
    use :: iso_c_binding
    type(c_ptr) polympo_createMPMesh
    integer(c_int), value :: setMeshOption, setMPOption
  end function
  !---------------------------------------------------------------------------
  !> @brief delete MPMesh object
  !> @param mpmesh(in/out) MPMesh object
  !---------------------------------------------------------------------------
  subroutine polympo_deleteMPMesh(mpMesh) bind(C, NAME='polympo_deleteMPMesh_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the MPI communicator used by polympo
  !> @param comm(in) MPI communicator
  !---------------------------------------------------------------------------
  subroutine polympo_setMPICommunicator(comm) &
             bind(C, NAME='polympo_setMPICommunicator_f')
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
             bind(C, NAME='polympo_createMPs_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numElms
    integer(c_int), value :: numMPs
    type(c_ptr), intent(in), value :: mpsPerElm
    type(c_ptr), intent(in), value :: mp2Elm
    type(c_ptr), intent(in), value :: isMPActive
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief move MPs to a new element, add new MPs, or delete MPs
  !> @brief the fields associated with the MPs are NOT initialized
  !> @param mpmesh(in/out) MPMesh object
  !> @param numMPs(in) total number of MPs, total = number of active + number of inactive
  !> @param allMP2Elm(in) the target element for each MP (length of numMPs)
  !> @param addedMPMask(in) set to 1 for each new MP, 0 otherwise (length of numMPs)
  !---------------------------------------------------------------------------
  subroutine polympo_startRebuildMPs(mpMesh, numMPs, allMP2Elm, addedMPMask) &
    bind(C, NAME='polympo_startRebuildMPs_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numMPs
    type(c_ptr), intent(in), value :: allMP2Elm
    type(c_ptr), intent(in), value :: addedMPMask
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief called after startRebuild()
  !> @brief called after initializing MP fields
  !> @param mpmesh(in/out) MPMesh object
  !---------------------------------------------------------------------------
  subroutine polympo_finishRebuildMPs(mpMesh) &
    bind(C, NAME='polympo_finishRebuildMPs_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief Stores pointer to appID data structure and a function to retrieve them used in migration
  !> @param mpmesh(in/out) MPMesh object
  !> @param getNext(in) Pointer to function that returns next App IDs
  !> @param appIDs(in) Pointer to opaque data application data structure (that may contain all available app IDs)
  !---------------------------------------------------------------------------
  subroutine polympo_setAppIDFunc(mpMesh, getNext, appIDs) &
    bind(C, NAME='polympo_setAppIDFunc_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    type(c_funptr), value :: getNext
    type(c_ptr), value :: appIDs
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the current element ID MP array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param numMPs(in) length of array, number of the MPs
  !> @param array(in/out) output MP element ID 1D array (numMPs), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPCurElmID(mpMesh, numMPs, array) &
             bind(C, NAME='polympo_getMPCurElmID_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numMPs
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the mp lat lon is rotational or normal
  !> @param mpmesh(in/out) MPMesh object
  !> @param isRotateFlag(in) Flag>0 = True, otherwise False.
  !---------------------------------------------------------------------------
  subroutine polympo_setMPLatLonRotatedFlag(mpMesh, isRotateFlag) &
             bind(C, NAME='polympo_setMPLatLonRotatedFlag_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: isRotateFlag
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the MP positions array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 3
  !> @param numMPs(in) number of the MPs
  !> @param array(in) MP current position 2D array (3,numMPs), allocated by user on host
  !---------------------------------------------------------------------------
  subroutine polympo_setMPPositions(mpMesh, nComps, numMPs, array) &
             bind(C, NAME='polympo_setMPPositions_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
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
             bind(C, NAME='polympo_getMPPositions_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the MP latitude and longtitude array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 2
  !> @param numMPs(in) number of the MPs
  !> @param array(in)  input MP current lat and lon 2D array (2,numMPs),
  !>                   allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_setMPRotLatLon(mpMesh, nComps, numMPs, array) &
             bind(C, NAME='polympo_setMPRotLatLon_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the MP latitude and longtitude array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 2
  !> @param numMPs(in) number of the MPs
  !> @param array(in/out) output MP current lat and lon 2D array (2,numMPs),
  !>                      allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPRotLatLon(mpMesh, nComps, numMPs, array) &
             bind(C, NAME='polympo_getMPRotLatLon_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mass MP array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 1
  !> @param numMPs(in) number of the MPs
  !> @param array(in) input MP Mass 1D array (numMPs)
  !---------------------------------------------------------------------------
  subroutine polympo_setMPMass(mpMesh, nComps, numMPs, array) &
             bind(C, NAME='polympo_setMPMass_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the Mass MP array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 1
  !> @param numMPs(in) number of the MPs
  !> @param array(in/out) output MP Mass 1D array (numMPs), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPMass(mpMesh, nComps, numMPs, array) &
             bind(C, NAME='polympo_getMPMass_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the velocity MP array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 2
  !> @param numMPs(in) number of the MPs
  !> @param array(in) input MP velocity 1D array (numMPs*2)
  !---------------------------------------------------------------------------
  subroutine polympo_setMPVel(mpMesh, nComps, numMPs, array) &
             bind(C, NAME='polympo_setMPVel_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the velocity MP array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nComps(in) number of components, should always be 2
  !> @param numMPs(in) number of the MPs
  !> @param array(in/out) output MP velocity 1D array (numMPs*2), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPVel(mpMesh, nComps, numMPs, array) &
             bind(C, NAME='polympo_getMPVel_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, numMPs
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief Enable the setting of mesh topology (number of entities and entity adjacencies). 
  !>        By default, modifying the mesh topology without calling this function first will result
  !>        in a runtime failure.
  !> @param mpMesh(in/out) the MPMesh is valid/created
  !---------------------------------------------------------------------------
  subroutine polympo_startMeshFill(mpMesh) &
             bind(C, NAME='polympo_startMeshFill_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief Disable the modification of mesh topology (number of entities and entity adjacencies). 
  !>        Once the topology is set, call this function to prevent unexpected/accidental modifications.
  !> @param mpMesh(in/out) the MPMesh is valid/created
  !---------------------------------------------------------------------------
  subroutine polympo_endMeshFill(mpMesh) &
             bind(C, NAME='polympo_endMeshFill_f')
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
             bind(C, NAME='polympo_checkMeshMaxSettings_f')
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
             bind(C, NAME='polympo_setMeshTypeGeneralPoly_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh type to CVT polygonal
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshTypeCVTPoly(mpMesh) &
             bind(C, NAME='polympo_setMeshTypeCVTPoly_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh geometry type to planar
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshGeomTypePlanar(mpMesh) &
             bind(C, NAME='polympo_setMeshGeomTypePlanar_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh geometry type to spherical
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshGeomTypeSpherical(mpMesh) &
             bind(C, NAME='polympo_setMeshGeomTypeSpherical_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh sphere radius
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshSphereRadius(mpMesh,sphereRadius) &
             bind(C, NAME='polympo_setMeshSphereRadius_f')
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
             bind(C, NAME='polympo_setMeshNumVtxs_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numVtxs
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the number of mesh vertices
  !> @param mpMesh(in) mpMesh object
  !> @param numVtxs(out) the number of vertices
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshNumVtxs(mpMesh, numVtx) &
            bind(C, NAME = 'polympo_getMeshNumVtxs_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), intent(inout) :: numVtx
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the number of elements of the mesh
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpMesh(in/out) mpMesh object 
  !> @param numElms(in) the number of elements
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshNumElms(mpMesh,numElms) &
             bind(C, NAME='polympo_setMeshNumElms_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numElms
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the number of mesh elements
  !> @param mpMesh(in) mpMesh object
  !> @param numVtxs(out) the number of elements
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshNumElms(mpMesh, numElm) &
            bind(C, NAME = 'polympo_getMeshNumElms_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), intent(inout) :: numElm
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh number of edges per element
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpmesh(in/out) MPMesh object
  !> @param nCells(in) length of array (numElms)
  !> @param nEdgesOnCell(in) number of edges per element
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshNumEdgesPerElm(mpMesh, nCells, nEdgesOnCell) &
             bind(C, NAME='polympo_setMeshNumEdgesPerElm_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nCells
    type(c_ptr), intent(in), value :: nEdgesOnCell
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh element to vertices connectivity
  !>        modifies mesh topology polympo_startMeshFill required
  !> @param mpmesh(in/out) MPMesh object
  !> @param maxEdges,nCells(in) length of array in each direction
  !> @param verticesOnCell(in) element to vertices connectivity 2D array 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshElm2VtxConn(mpMesh, maxEdges, nCells, verticesOnCell) &
             bind(C, NAME='polympo_setMeshElm2VtxConn_f')
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
             bind(C, NAME='polympo_setMeshElm2ElmConn_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: maxEdges, nCells
    type(c_ptr), intent(in), value :: cellsOnCell
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh vertices coordinates
  !> @param mpmesh(in/out) MPMesh object
  !> @param nVertices(in) length of array in 
  !> @param x/y/zArray(in) the 1D arrays of vertices coordinates
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshVtxCoords(mpMesh, nVertices, xArray, yArray, zArray) &
             bind(C, NAME='polympo_setMeshVtxCoords_f')
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
             bind(C, NAME='polympo_getMeshVtxCoords_f')
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
  subroutine polympo_setMeshVtxRotLat(mpMesh, nVertices, latitude) &
             bind(C, NAME='polympo_setMeshVtxRotLat_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nVertices
    type(c_ptr), intent(in), value :: latitude
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the polympo mesh vertices latitude and longitude
  !> @param mpmesh(in/out) MPMesh object
  !> @param nVertices(in) length of array in, use for assertion
  !> @param latitude/longitude(in/out) the 1D arrays of vertices lat/lon
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshVtxRotLat(mpMesh, nVertices, latitude) &
             bind(C, NAME='polympo_getMeshVtxRotLat_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nVertices
    type(c_ptr), value :: latitude
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the vertices velocity from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param nVertices(in) numVertices
  !> @param uVel(in) vertices u-component velocity 1D array (numVtx)
  !> @param vVel(in) vertices v-component velocity 1D array (numVtx)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshVel(mpMesh, nVertices, uVel, vVel) &
             bind(C, NAME='polympo_setMeshVel_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nVertices
    type(c_ptr), intent(in), value :: uVel, vVel
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the vertices velocity from polyMPO
  !> @param mpmesh(in/out) MPMesh object
  !> @param nVertices(in) numVertices
  !> @param uVel(in/out) output vertices u-component velocity
  !>        1D array (numVtx), allocated by user
  !> @param vVel(in/out) output vertices v-component velocity
  !>        1D array (numVtx), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshVel(mpMesh, nVertices, uVel, vVel) &
             bind(C, NAME='polympo_getMeshVel_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nVertices
    type(c_ptr), value :: uVel, vVel
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
             bind(C, NAME='polympo_setMeshOnSurfVeloIncr_f')
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
             bind(C, NAME='polympo_getMeshOnSurfVeloIncr_f')
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
             bind(C, NAME='polympo_setMeshOnSurfDispIncr_f')
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
             bind(C, NAME='polympo_getMeshOnSurfDispIncr_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: nComps, nVertices
    type(c_ptr), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief calculate the MPs from given mesh vertices rotational latitude
  !>        longitude, update the MP slices
  !>        MPs MUST have rotated flag set to True(>0)
  !> @param mpmesh(in/out) MPMesh object
  !---------------------------------------------------------------------------
  subroutine polympo_push(mpMesh) &
             bind(C, NAME='polympo_push_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Reconstruction of MP fields options
  !> @param mpmesh(in/out) MPMesh object
  !> @param reconsOption(in) 1 = True, otherwise False.
  !> TODO: not support yet!
  !---------------------------------------------------------------------------
  subroutine polympo_setReconstructionOption(mpMesh, reconsOption) &
             bind(C, NAME='polympo_setReconstructionOption_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: reconsOption
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief directly call the reconstruct of the MP fields to mesh fields
  !> @param mpmesh(in/out) MPMesh object
  !---------------------------------------------------------------------------
  subroutine polympo_reconstruct(mpMesh) &
             bind(C, NAME='polympo_reconstruct_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief start the reconstruction of MP Mass to Mesh Vertices
  !> @param mpmesh(in/out) MPMesh object
  !> @param order Order of the reconstruction field
  !> @param meshEntType Mesh entity type
  !> TODO not support yet!
  !---------------------------------------------------------------------------
  subroutine polympo_reconstructMass(mpMesh, order, meshEntType) &
             bind(C, NAME='polympo_reconstructMass_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: order, meshEntType
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief start the reconstruction of MP Velocity to Mesh Vertices
  !> @param mpmesh(in/out) MPMesh object
  !> @param order Order of the reconstruction field
  !> @param meshEntType Mesh entity type
  !> TODO not support yet!
  !---------------------------------------------------------------------------
  subroutine polympo_reconstructVel(mpMesh, order, meshEntType) &
             bind(C, NAME='polympo_reconstructVel_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: order, meshEntType
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief start the reconstruction of MP Strain Rate to Mesh Vertices
  !> @param mpmesh(in/out) MPMesh object
  !> @param order Order of the reconstruction field
  !> @param meshEntType Mesh entity type
  !> TODO not support yet!
  !---------------------------------------------------------------------------
  subroutine polympo_reconstructStrainRate(mpMesh, order, meshEntType) &
             bind(C, NAME='polympo_reconstructStrainRate_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: order, meshEntType
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief start the reconstruction of MP Stress to Mesh Vertices
  !> @param mpmesh(in/out) MPMesh object
  !> @param order Order of the reconstruction field
  !> @param meshEntType Mesh entity type
  !> TODO not support yet!
  !---------------------------------------------------------------------------
  subroutine polympo_reconstructStress(mpMesh, order, meshEntType) &
             bind(C, NAME='polympo_reconstructStress_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: order, meshEntType
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
