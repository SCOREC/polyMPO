#ifndef POLYMPMTEST_MESH_H
#define POLYMPMTEST_MESH_H

namespace polyMpmTest{

#define maxVtxsPerElm 8
#define maxElmsPerVtx 5

using IntVtxView = Kokkos::View<int*[maxVtxsPerElm+1]>;
using IntElmsPerVtxView = Kokkos::View<int*[maxElmsPerVtx+1]>;

class Mesh {
  private:
    VectorView vtxCoords_;
    IntVtxView elm2VtxConn_;
    int numElms_;
    Int ElmsPerVtxView vtx2ElmConn_;

  public:
    Mesh();
    Mesh(/*TODO: variables*/);
    ~Mesh();

    int getNumVertices();
    int getNumElements();
    VectorView getVtxCoords();
    IntVertiView getElm2VtxConn();
    IntElemsPerVertView getVtx2ElmConn();

    void setVtx2ElmConn(IntElmsPerVertView vertex2Elems);  
};


}i

#endif

