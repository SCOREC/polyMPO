#ifndef POLYMPMTEST_MESH_H
#define POLYMPMTEST_MESH_H

namespace polyMpmTest{

#define maxVtxsPerElm 8
#define maxElmsPerVtx 5

using IntVtx2ElmView = Kokkos::View<int*[maxVtxsPerElm+1]>;
using IntElm2VtxView = Kokkos::View<int*[maxElmsPerVtx+1]>;

class Mesh {
  private:
    int numVtxs_;
    int numElms_;
    VectorView vtxCoords_;
    IntVtx2ElmView elm2VtxConn_;
    IntElm2VtxView vtx2ElmConn_;

  public:
    Mesh();
    Mesh(/*TODO: variables*/);
    ~Mesh();

    int getNumVertices();
    int getNumElements();
    VectorView getVtxCoords();
    IntVertiView getElm2VtxConn();
    IntElemsPerVertView getVtx2ElmConn();

    void setVtx2ElmConn(IntElm2VtxView vtx2ElmConn);  
};


}

#endif

