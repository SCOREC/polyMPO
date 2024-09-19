subroutine calcSurfDispIncr(mpMesh, latVertex, lonVertex, nEdgesOnCell, verticesOnCell, nVertices, sphereRadius, scale)
  use :: polympo
  use :: readMPAS
  use :: iso_c_binding
  implicit none

  real(kind=MPAS_RKIND), dimension(:), pointer :: latVertex, lonVertex
  real(kind=MPAS_RKIND) :: maxlon, minlon, deltaLon, sphereRadius
  integer, dimension(:), pointer :: nEdgesOnCell
  integer, dimension(:,:), pointer :: verticesOnCell
  integer :: i, j, nVertices, nCompsDisp, scale_use
  real(kind=MPAS_RKIND), dimension(:,:), pointer :: dispIncr
  type(c_ptr) :: mpMesh
  INTEGER, INTENT(IN), OPTIONAL :: scale
  
  IF (PRESENT(scale)) THEN
    scale_use=scale
  ELSE
    scale_use=1
  END IF

  nCompsDisp = 2
  allocate(dispIncr(nCompsDisp,nVertices))

  ! check first element/cell for delta
  maxlon = minval(lonVertex)
  minlon = maxval(lonVertex)
  do i = 1, nEdgesOnCell(1)
    j = verticesOnCell(i,1)
    if(maxlon .lt. lonVertex(j)) then
      maxlon = lonVertex(j)
    endif
    if(minlon .gt. lonVertex(j)) then
      minlon = lonVertex(j)
    endif
  end do

  deltaLon = maxlon - minlon

  do i = 1,nVertices
    dispIncr(1,i) = scale_use*sphereRadius*cos(latVertex(i))*deltaLon
    dispIncr(2,i) = 0.0_MPAS_RKIND
  end do
  call polympo_setMeshVtxOnSurfDispIncr(mpMesh,nCompsDisp,nVertices,c_loc(dispIncr))

  deallocate(dispIncr)
end subroutine

