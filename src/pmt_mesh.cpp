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

  int xVertexID, yVertexID, verticesOnCellID, cellsOnVertexID, nEdgesOnCellID, cellsOnEdgeID, verticesOnEdgeID, cellsOnCellID, edgesOnCellID;
  double *xVertex;
  double *yVertex;
  int *verticesOnCell; //[maxEdges,nCells]
  int *cellsOnVertex;  //[3,nVertices]
  int *nEdgesOnCell;   //[nCells]
  int *cellsOnEdge; //[nEdges*2]
  int *verticesOnEdge; //[nEdges*2]
  int *cellsOnCell; //[nCells*maxEdges]
  int *edgesOnCell;

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
  if ((retval = nc_inq_varid(ncid, "verticesOnEdge", &verticesOnEdgeID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "cellsOnCell", &cellsOnCellID)))
    ERRexit(retval);
  if ((retval = nc_inq_varid(ncid, "edgesOnCell", &edgesOnCellID)))
    ERRexit(retval);

  xVertex = new double[nVertices];
  yVertex = new double[nVertices];
  verticesOnCell = new int[nCells * maxEdges];

  cellsOnVertex = new int[vertexDegree * nVertices]; // vertex dimension is vertexDegree
  nEdgesOnCell = new int[nCells];
  cellsOnEdge = new int[2 * nEdges];
  verticesOnEdge = new int[2 * nEdges];
  cellsOnCell = new int[nCells * maxEdges];
  edgesOnCell = new int[nCells * maxEdges];

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
  if ((retval = nc_get_var(ncid, verticesOnEdgeID, verticesOnEdge)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, cellsOnCellID, cellsOnCell)))
    ERRexit(retval);
  if ((retval = nc_get_var(ncid, edgesOnCellID, edgesOnCell)))
    ERRexit(retval);
//TODO: add verticesOnEdge[nEdge*2] to make edgeToEdge
   /* 
   for(int i=0; i<10; i++){
        //i=cell
        printf("%2d ",i+1);
        int n = nEdgesOnCell[i];
        for (int j = 0; j < n; j++){
            printf("| %4d= %4d,%4d ",verticesOnCell[i*maxEdges+j],verticesOnEdge[(edgesOnCell[i*maxEdges+(j+1)%n]-1)*2],verticesOnEdge[(edgesOnCell[i*maxEdges+(j+1)%n]-1)*2+1]);
        }
        printf("\n   ");
        for (int j = 0; j < n; j++){
            printf("| %4d= %4d,%4d ",cellsOnCell[i*maxEdges+(j+1)%n], cellsOnEdge[(edgesOnCell[i*maxEdges+(j+1)%n]-1)*2], cellsOnEdge[(edgesOnCell[i*maxEdges+(j+1)%n]-1)*2+1]);
        }
        printf("\n");
    }//==*/
/*
    for(int i=0; i<nEdges; i++){
        int v1 = verticesOnEdge[i*2];
        int v2 = verticesOnEdge[i*2+1];
        if(v1 == 1 || v2 == 1)
            printf("%d:%d,%d|%d,%d\n",i,v1,v2,cellsOnEdge[i*2],cellsOnEdge[i*2+1]);
    }
    
    for(int i=0; i<nCells; i++){
        int e1 = edgesOnCell[i*4];
        int e2 = edgesOnCell[i*4+1];
        int e3 = edgesOnCell[i*4+2];
        int e4 = edgesOnCell[i*4+3];
        if(e1 == 7529 || e2 == 7529 || e3 == 7529 || e4 == 7529)
            printf("%d:%d,%d,%d,%d\n",i,e1,e2,e3,e4);
    }
///====*/
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
  IntElm2ElmView elm2ElmConn("elementToElmentConnection", nCells);
  IntElm2ElmView::HostMirror h_elm2ElmConn = Kokkos::create_mirror_view(elm2ElmConn);

  for (int i = 0; i < nCells; i++)
  {
    int n = nEdgesOnCell[i];
    h_elm2VtxConn(i, 0) = n;
    h_elm2ElmConn(i, 0) = n;
    for (int j = 0; j < n; j++)
    {
      h_elm2VtxConn(i, j+1) = verticesOnCell[i * maxEdges + j];
      h_elm2ElmConn(i, j+1) = cellsOnCell[i*maxEdges+(j+1)%n]-1;
    }
  }


  Kokkos::deep_copy(elm2VtxConn, h_elm2VtxConn);
  Kokkos::deep_copy(elm2ElmConn, h_elm2ElmConn);
/*
    printf("<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"%d\">\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",nVertices,nCells);
    for(int i=0; i<nVertices; i++){
        printf("          %f %f 0.0\n",xVertex[i],yVertex[i]);
    }
    printf("        </DataArray>\n      </Points>\n      <Polys>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int i=0; i<nCells; i++){
        printf("          ");
        for(int j=0; j< nEdgesOnCell[i]; j++){
            printf("%d ", h_elm2VtxConn(i,j+1)-1);
        } 
        printf("\n");
    }
    printf("        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    
    int count = 0;
    for(int i=0;i<nCells; i++){
        count += h_elm2VtxConn(i,0);
        printf("          %d\n",count);
    }
    printf("        </DataArray>\n      </Polys>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");
*/

  // delete dynamic allocation
  delete[] xVertex;
  delete[] yVertex;
  delete[] verticesOnCell;
  delete[] cellsOnVertex;

  return Mesh(nVertices,
              nCells,
              vtxCoords,
              elm2VtxConn,
              vtx2ElmConn,
              elm2ElmConn);  
}

} // namespace polyMpmTest
