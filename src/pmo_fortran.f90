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
  !> @brief set the 2d velocity MP array from a host array
  !> @param mpmesh(in/out) MPMesh object
  !> @param n(in) number of items
  !> @param m(in) components per item
  !> @param array(in) input MP velocity array
  !---------------------------------------------------------------------------
  subroutine polympo_setMP2dVelArray(mpMesh, n, m, array) &
             bind(C, NAME='polympo_setMP2dVelArray')
    use :: iso_c_binding
    type(c_ptr), value :: mpMesh
    integer(c_int), value :: n, m
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
  end interface
end module
