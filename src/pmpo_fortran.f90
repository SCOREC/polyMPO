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
  !> @brief set the velocity MP array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) half length of array
  !> @param array(in) input MP velocity 1D array (numMPs*2)
  !---------------------------------------------------------------------------
  subroutine polympo_setMPVelArray(mpMesh, n, array) &
             bind(C, NAME='polympo_setMPVelArray_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the velocity MP array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) half length of array
  !> @param array(in/out) output MP velocity 1D array (numMPs*2), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPVelArray(mpMesh, n, array) &
             bind(C, NAME='polympo_getMPVelArray_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    type(c_ptr), value :: array
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
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshTypeGeneralPoly(mpMesh) &
             bind(C, NAME='polympo_setMeshTypeGeneralPoly_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh type to CVT polygonal
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshTypeCVTPoly(mpMesh) &
             bind(C, NAME='polympo_setMeshTypeCVTPoly_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh geometry type to planar
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshGeomTypePlanar(mpMesh) &
             bind(C, NAME='polympo_setMeshGeomTypePlanar_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh geometry type to spherical
  !> @param mpMesh(in/out) mpMesh object 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshGeomTypeSpherical(mpMesh) &
             bind(C, NAME='polympo_setMeshGeomTypeSpherical_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh sphere radius
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
  !> @brief set the number of elements of the mesh
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
  !> @brief set the polympo mesh element to vertices connectivity
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
  !> @brief set the polympo mesh number of edges per element
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
  !> @brief set the velocity mesh array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) half length of array (n = numVtx)
  !> @param array(in) input mesh velocity 1D array (numVtx*2)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshVelArray(mpMesh, n, array) &
             bind(C, NAME='polympo_setMeshVelArray_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the velocity mesh array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) half length of the array
  !> @param array(in/out) output mesh velocity 1D array (numVtx*2), allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshVelArray(mpMesh, n, array) &
             bind(C, NAME='polympo_getMeshVelArray_f')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
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
