#include "pmt_mesh.hpp"

#ifdef POLYMPI_HAS_NETCDF
#include <netcdf.h>

#define ERRexit(e){ printf("Error: %s\n", nc_strerror(e)); exit(2);}

namespace polyMpmTest{

Mesh Mesh::readMPASMesh(std::string filename){
  int ncid;
  int retval;
  if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
    ERRexit(retval);
  Mesh mesh = readMPASMesh(ncid);
  if ((retval = nc_close(ncid)))
    ERRexit(retval);
  return mesh;
}

Mesh Mesh::readMPASMesh(int ncid){
  int retval,
      nCells, nCellsID,
      nVertices, nVerticesID,
      maxEdges, maxEdgesID,
      vertexDegree, vertexDegreeID;
  size_t temp;

  int xVertexID, yVertexID, zVertexID, lonVertexID, latVertexID,
      verticesOnCellID, cellsOnVertexID, nEdgesOnCellID;
  double *xVertex;
  double *yVertex;
  double *zVertex;     // nVertices
  double *lonVertex;   // 2d y
  double *latVertex;   // 2d x
  int *verticesOnCell; //[maxEdges,nCells]
  int *cellsOnVertex;  //[3,nVertices]
  int *nEdgesOnCell;   //[nCells]
  if ((retval = nc_inq_dimid(ncid, "nCells", &nCellsID)))
    ERRexit(retval);
  if ((retval = nc_inq_dimid(ncid, "nVertices", &nVerticesID)))
    ERRexit(retval);
  if ((retval = nc_inq_dimid(ncid, "maxEdges", &maxEdgesID)))
    ERRexit(retval);
  if ((retval = nc_inq_dimid(ncid, "vertexDegree", &vertexDegreeID)))
    ERRexit(retval);

  if ((retval = nc_inq_dimlen(ncid, nCellsID, &temp)))
    ERRexit(retval);
  nCells = temp;
  if ((retval = nc_inq_dimlen(ncid, nVerticesID, &temp)))
    ERRexit(retval);
  nVertices = temp;
  if ((retval = nc_inq_dimlen(ncid, maxEdgesID, &temp)))
    ERRexit(retval);
  maxEdges = temp;
  if ((retval = nc_inq_dimlen(ncid, vertexDegreeID, &temp)))
    ERRexit(retval);
  vertexDegree = temp;

  if ((retval = nc_inq_varid(ncid, "xVertex", &xVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "yVertex", &yVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "zVertex", &zVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "lonVertex", &lonVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "latVertex", &latVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "verticesOnCell", &verticesOnCellID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "cellsOnVertex", &cellsOnVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "nEdgesOnCell", &nEdgesOnCellID)))
    ERRexit(retval);

  xVertex = new double[nVertices];
  yVertex = new double[nVertices];
  zVertex = new double[nVertices];
  lonVertex = new double[nVertices];
  latVertex = new double[nVertices];
  verticesOnCell = new int[maxEdges * nCells];

  cellsOnVertex = new int[vertexDegree * nVertices]; // vertex dimension is vertexDegree
  nEdgesOnCell = new int[nCells];

  if (maxEdges > maxVtxsPerElm)
  {
    perror("maxEdges out of bound!\n");
    exit(1);
  }
  if (vertexDegree > maxElmsPerVtx)
  {
    perror("vertexDegree is too large!\n");
    exit(1);
  }

  if ((retval = nc_get_var(ncid, xVertexID, xVertex)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, yVertexID, yVertex)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, zVertexID, zVertex)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, lonVertexID, lonVertex)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, latVertexID, latVertex)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, verticesOnCellID, verticesOnCell)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, cellsOnVertexID, cellsOnVertex)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, nEdgesOnCellID, nEdgesOnCell)))
    ERRexit(retval);

  Vec2dView vtxCoords("verticesCoordinates", nVertices);
  IntElm2VtxView vtx2ElmConn("vertexToElementsConnection", nVertices); // 4 = vertexDegree + 1

  Vec2dView::HostMirror h_vtxCoords = Kokkos::create_mirror_view(vtxCoords);
  IntElm2VtxView::HostMirror h_vtx2ElmConn = Kokkos::create_mirror_view(vtx2ElmConn);
  for (int i = 0; i < nVertices; i++)
  {
    h_vtxCoords(i) = Vec2d(xVertex[i], yVertex[i]);
    h_vtx2ElmConn(i, 0) = vertexDegree;
    for (int j = 0; j < vertexDegree; j++)
    {
      h_vtx2ElmConn(i, j + 1) = cellsOnVertex[i * 3 + j];
    }
  }
  Kokkos::deep_copy(vtxCoords, h_vtxCoords);
  Kokkos::deep_copy(vtx2ElmConn, h_vtx2ElmConn);

  IntVtx2ElmView elm2VtxConn("elementToVerticesConnection", nCells);
  IntVtx2ElmView::HostMirror h_elm2VtxConn = Kokkos::create_mirror_view(elm2VtxConn);

  for (int i = 0; i < nCells; i++)
  {
    h_elm2VtxConn(i, 0) = nEdgesOnCell[i];
    for (int j = 0; j < nEdgesOnCell[i]; j++)
    {
      h_elm2VtxConn(i, j + 1) = verticesOnCell[i * maxEdges + j];
    }
  }

  Kokkos::deep_copy(elm2VtxConn, h_elm2VtxConn);

  // delete dynamic allocation
  delete[] xVertex;
  delete[] yVertex;
  delete[] zVertex;
  delete[] lonVertex;
  delete[] latVertex;
  delete[] verticesOnCell;
  delete[] cellsOnVertex;

  return Mesh(nVertices,
              nCells,
              vtxCoords,
              elm2VtxConn,
              vtx2ElmConn);
}
#else

namespace polyMpmTest{
  Mesh Mesh::readMPASMesh(std::string) {
    return readMPASMesh(0);
  }
  Mesh Mesh::readMPASMesh(int){
    fprintf(stderr, "ERROR: readMPASMesh requires compiling with NetCDF enabled\n");
    exit(EXIT_FAILURE);
    return Mesh();
  }
#endif

} // namespace polyMpmTest
