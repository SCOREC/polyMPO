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
  !> @todo take mpi communicator from user
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
  !> @return mpmesh(in/out) MPMesh object
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
  !> @brief replicate a new MPMesh object
  !> @param mpmesh(in/out) MPMesh object
  !> @param scaleFactor(in/out) the scaleFactor of the new MPMesh
  !>        new MPMesh = MPMeshObj * scaleFactor
  !---------------------------------------------------------------------------
  function polympo_replicateMPMesh(mpMesh, scaleFactor) bind(C, NAME='polympo_replicateMPMesh')
    use :: iso_c_binding
    type(c_ptr) polympo_replicateMPMesh
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: scaleFactor
  end function
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
  !> @brief set the velocity MP array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) length of array
  !> @param array(in) input MP velocity array
  !---------------------------------------------------------------------------
  subroutine polympo_setMPVelArray(mpMesh, n, array) &
             bind(C, NAME='polympo_setMPVelArray')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    real(c_double), intent(in), dimension(n) :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the velocity MP array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) length of array
  !> @param array(in/out) output MP velocity array, allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMPVelArray(mpMesh, n, array) &
             bind(C, NAME='polympo_getMPVelArray')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    real(c_double), intent(inout), dimension(n) :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief check the Mesh is valid/runable in polympo
  !> @param mpMesh(in/out) the MPMesh is valid/created
  !> @param maxEdges(in) the maxEdges per Cell of your mesh 
  !> @param vertexDegree(in) the max vertexDegree of a vertex 
  !---------------------------------------------------------------------------
  subroutine polympo_checkMeshMaxSettings(mpMesh,maxEdges,vertexDegree) &
             bind(C, NAME='polympo_checkMeshMaxSettings')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: maxEdges, vertexDegree
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh Type
  !> @param mpMesh(in/out) mpMesh Object 
  !> @param numVtxs(in) the Mesh Type: (mesh_unrecognized = -1
  !>                                    mesh_general_polygoal = 0
  !>                                    mesh_CVT_polygonal = 1)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshType(mpMesh,meshType) &
             bind(C, NAME='polympo_setMeshType')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: meshType
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh geometry Type
  !> @param mpMesh(in/out) mpMesh Object 
  !> @param numVtxs(in) the Geom Type: (mesh_unrecognized = -1
  !>                                    mesh_general_polygoal = 0
  !>                                    mesh_CVT_polygonal = 1)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshGeomType(mpMesh,geomType) &
             bind(C, NAME='polympo_setMeshGeomType')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: geomType
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the Mesh sphere radius
  !> @param mpMesh(in/out) mpMesh Object 
  !> @param numVtxs(in) the Geom Type: (mesh_unrecognized = -1
  !>                                    mesh_general_polygoal = 0
  !>                                    mesh_CVT_polygonal = 1)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshSphereRadius(mpMesh,sphereRadius) &
             bind(C, NAME='polympo_setMeshSphereRadius')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    real(c_double), value :: sphereRadius
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the number of vetices of the mesh
  !> @param mpMesh(in/out) mpMesh Object 
  !> @param numVtxs(in) the number of vertices need to set
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshNumVtxs(mpMesh,numVtxs) &
             bind(C, NAME='polympo_setMeshNumVtxs')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numVtxs
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the number of elements of the mesh
  !> @param mpMesh(in/out) mpMesh Object 
  !> @param numElms(in) the number of elements
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshNumElms(mpMesh,numElms) &
             bind(C, NAME='polympo_setMeshNumElms')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: numElms
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh vertices coordinates
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) number of the vertices coordinates 
  !> @param x/y/zArray(in) the arrays of vertices coordinates
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshVtxCoords(mpMesh, n, xArray, yArray, zArray) &
             bind(C, NAME='polympo_setMeshVtxCoords')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    type(c_ptr), intent(in), value :: xArray, yArray, zArray
    !todo: change all the 1d arrays to pointers
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh element to vertices connectivity
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) length of array (numElms*maxEdges)
  !> @param array(in) element to vertices connectivity array (verticesOnCell)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshElm2VtxConn(mpMesh, m, n, array) &
             bind(C, NAME='polympo_setMeshElm2VtxConn')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: m, n
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh element to elements connectivity
  !> @param mpmesh(in/out) MPMesh object
  !> @param m,n(in) length of array (numElms*maxEdges)
  !> @param array(in) element to elements connectivity array (cellsOnCell)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshElm2ElmConn(mpMesh, m, n, array) &
             bind(C, NAME='polympo_setMeshElm2ElmConn')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: m, n
    type(c_ptr), intent(in), value :: array
  end subroutine

  subroutine polympo_setMeshNumEdgesPerElm(mpMesh, m, array) &
             bind(C, NAME='polympo_setMeshNumEdgesPerElm')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: m
    type(c_ptr), intent(in), value :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the velocity mesh array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) length of array
  !> @param array(in) input mesh velocity array
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshVelArray(mpMesh, n, array) &
             bind(C, NAME='polympo_setMeshVelArray')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    real(c_double), intent(in), dimension(n) :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief get the velocity mesh array from a polympo array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) length of array
  !> @param array(in/out) output mesh velocity array, allocated by user
  !---------------------------------------------------------------------------
  subroutine polympo_getMeshVelArray(mpMesh, n, array) &
             bind(C, NAME='polympo_getMeshVelArray')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    real(c_double), intent(inout), dimension(n) :: array
  end subroutine
  end interface
end module
