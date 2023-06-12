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
  !---------------------------------------------------------------------------
  subroutine polympo_initialize() bind(C, NAME='polympo_initialize')
    use :: iso_c_binding
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief finalize polympo, no polympo apis may be called after this
  !---------------------------------------------------------------------------
  subroutine polympo_finalize() bind(C, NAME='polympo_finalize')
    use :: iso_c_binding
  end subroutine
  !---------------------------------------------------------------------------
  !> @brief initialize polympo, call this before any other polympo api
  !---------------------------------------------------------------------------
  function polympo_createMpMesh() bind(C, NAME='polympo_createMpMesh')
    use :: iso_c_binding
    type(c_ptr) polympo_createMpMesh
  end function

  !---------------------------------------------------------------------------
  !> @brief modify the specified array
  !> @param n(in) length of array
  !> @param ranks(in/out) array to modify
  !---------------------------------------------------------------------------
  subroutine polympo_modifyArray(n, array) &
             bind(C, NAME='polympo_modifyArray')
    use :: iso_c_binding
    integer(c_int), value :: n
    real(c_double), intent(inout), dimension(n) :: array
  end subroutine
  end interface
end module
