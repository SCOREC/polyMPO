#include "pmt_mesh.hpp"
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
      vertexDegree, vertexDegreeID,
      nEdges, nEdgesID;
  size_t temp;

  int xVertexID, yVertexID, verticesOnCellID, cellsOnVertexID, nEdgesOnCellID, cellsOnEdgeID;
  double *xVertex;
  double *yVertex;
  int *verticesOnCell; //[maxEdges,nCells]
  int *cellsOnVertex;  //[3,nVertices]
  int *nEdgesOnCell;   //[nCells]
  int *cellsOnEdge; //[nEdges*maxEdges2]

  if ((retval = nc_inq_dimid(ncid, "nCells", &nCellsID)))
    ERRexit(retval);
  if ((retval = nc_inq_dimid(ncid, "nVertices", &nVerticesID)))
    ERRexit(retval);
  if ((retval = nc_inq_dimid(ncid, "maxEdges", &maxEdgesID)))
    ERRexit(retval);
  if ((retval = nc_inq_dimid(ncid, "vertexDegree", &vertexDegreeID)))
    ERRexit(retval);
  if ((retval = nc_inq_dimid(ncid, "nEdges", &nEdgesID)))
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
  if ((retval = nc_inq_dimlen(ncid, nEdgesID, &temp)))
    ERRexit(retval);
  nEdges = temp;

  if ((retval = nc_inq_varid(ncid, "xVertex", &xVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "yVertex", &yVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "verticesOnCell", &verticesOnCellID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "cellsOnVertex", &cellsOnVertexID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "nEdgesOnCell", &nEdgesOnCellID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "cellsOnEdge", &cellsOnEdgeID)))
    ERRexit(retval);

  xVertex = new double[nVertices];
  yVertex = new double[nVertices];
  verticesOnCell = new int[maxEdges * nCells];

  cellsOnVertex = new int[vertexDegree * nVertices]; // vertex dimension is vertexDegree
  nEdgesOnCell = new int[nCells];
  cellsOnEdge = new int[2 * nEdges];

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
  if ((retval = nc_get_var(ncid, verticesOnCellID, verticesOnCell)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, cellsOnVertexID, cellsOnVertex)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, nEdgesOnCellID, nEdgesOnCell)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, cellsOnEdgeID, cellsOnEdge)))
    ERRexit(retval);
//TODO: add cellOnVertex[nVertices*VertexDegree] and verticesOnEdge[nEdge*2] to make edgeToEdge


  Vector2View vtxCoords("verticesCoordinates", nVertices);
  IntElm2VtxView vtx2ElmConn("vertexToElementsConnection", nVertices); // 4 = vertexDegree + 1

  Vector2View::HostMirror h_vtxCoords = Kokkos::create_mirror_view(vtxCoords);
  IntElm2VtxView::HostMirror h_vtx2ElmConn = Kokkos::create_mirror_view(vtx2ElmConn);
  for (int i = 0; i < nVertices; i++)
  {
    h_vtxCoords(i) = Vector2(xVertex[i], yVertex[i]);
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
  delete[] verticesOnCell;
  delete[] cellsOnVertex;

  return Mesh();
  /*
  return Mesh(nVertices,
              nCells,
              vtxCoords,
              elm2VtxConn,
              vtx2ElmConn);  
  */
}

} // namespace polyMpmTest
