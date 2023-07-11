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
  function polympo_createMpMesh() bind(C, NAME='polympo_createMpMesh')
    use :: iso_c_binding
    type(c_ptr) polympo_createMpMesh
  end function
  !---------------------------------------------------------------------------
  !> @brief delete MPMesh object
  !> @param mpmesh(in/out) MPMesh object
  !---------------------------------------------------------------------------
  subroutine polympo_deleteMpMesh(mpMesh) bind(C, NAME='polympo_deleteMpMesh')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
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
  !---------------------------------------------------------------------------
  !> @brief set the MPI communicator used by polympo
  !> @param comm(in) MPI communicator
  !---------------------------------------------------------------------------
  subroutine polympo_setCommunicator(comm) &
             bind(C, NAME='polympo_setCommunicator')
    use :: iso_c_binding
    integer(c_int), value :: comm    
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
    real(c_double), intent(in), dimension(:,:) :: xArray, yArray, zArray
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh element to vertices connectivity
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) length of array
  !> @param array(in) element to vertices connectivity array (verticesOnCell)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshElm2VtxConn(mpMesh, n, array) &
             bind(C, NAME='polympo_setMeshElm2VtxConn')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    integer(c_int), intent(in), dimension(:,:) :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh vertex to elements connectivity
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) length of array
  !> @param array(in) vertex to elements connectivity array (cellsOnVertex) 
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshVtx2ElmConn(mpMesh, n, array) &
             bind(C, NAME='polympo_setMeshVtx2ElmConn')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    integer(c_int), intent(in), dimension(:,:) :: array
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief set the polympo mesh element to elements connectivity
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) length of array
  !> @param array(in) element to elements connectivity array (cellsOnCell)
  !---------------------------------------------------------------------------
  subroutine polympo_setMeshElm2ElmConn(mpMesh, n, array) &
             bind(C, NAME='polympo_setMeshElm2ElmConn')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n
    integer(c_int), intent(in), dimension(:,:) :: array
  end subroutine

!=============================================================================
  !---------------------------------------------------------------------------
  !> @brief read MPAS mesh(Version 1.0)
  !> @param filename(in)
  !> 
  !---------------------------------------------------------------------------
  subroutine polympo_readMPASMesh(filename, &
                                  nCells, nVertices, nEdgesOnCell, &
                                  xVertex, yVertex, zVertex, &
                                  verticesOnCell, cellsOnVertex, cellsOnCell)
    use :: netcdf
    implicit none
    
    character (len=*), intent(in) :: filename
    integer, intent(inout) :: nCells, nVertices
    integer, dimension(:), pointer :: nEdgesOnCell
    double precision, dimension(:), pointer :: xVertex, yVertex, zVertex
    integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex, cellsOnCell

    integer :: maxEdge, vertexDegree
    integer :: ncid, nCellsID, nVerticesID, maxEdgeID, vertexDegreeID, &
               nEdgesOnCellID, xVertexID, yVertexID, zVertexID, &
               verticesOnCellID, cellsOnVertexID, cellsOnCellID
    
    status = nf90_open(path=trim(filename), mode=nf90_nowrite, ncid=ncid)
    !if (status /= nf90_noerr) then
    !    write(0,*) 'polympo_readMPASMesh: Error occured when opening MPAS grid:'//filename
    !    write(0,*) trim(nf90_strerror(status))
    !    stop 
    !end if    
    
    
  end subroutine polympo_readMPASMesh
  end interface
end module
