#ifndef POLYMPMTEST_MESH_H
#define POLYMPMTEST_MESH_H

#include "utils.hpp"

namespace polyMpmTest{

#define maxVtxsPerElm 8
#define maxElmsPerVtx 5

using IntVtx2ElmView = Kokkos::View<int*[maxVtxsPerElm+1]>;
using IntElm2VtxView = Kokkos::View<int*[maxElmsPerVtx+1]>;

class Mesh {
  private:
    int numVtxs_;
    int numElms_;
    Vector2View vtxCoords_;
    IntVtx2ElmView elm2VtxConn_;
    IntElm2VtxView vtx2ElmConn_;

  public:
    Mesh(){};
    Mesh( int numVtxs,
          int numElms,
          Vector2View vtxCoords,
          IntVtx2ElmView elm2VtxConn,
          IntElm2VtxView vtx2ElmConn):
          numVtxs_(numVtxs),
          numElms_(numElms),
          vtxCoords_(vtxCoords),
          elm2VtxConn_(elm2VtxConn),
          vtx2ElmConn_(vtx2ElmConn){};

    int getNumVertices() { return numVtxs_; }
    int getNumElements() { return numElms_; }
    Vector2View getVtxCoords() { return vtxCoords_; }
    IntVertiView getElm2VtxConn() { return elm2VtxConn_; }
    IntElemsPerVertView getVtx2ElmConn() { return vtx2ElmConn_; }

    void setVtx2ElmConn(IntElm2VtxView vtx2ElmConn) { vtx2ElmConn_ = vtx2ElmConn; }
    void setElm2VtxConn(IntVtx2ElmView elm2VtxConn) { elm2VtxConn_ = elm2VtxConn; }
};

Mesh readMPASMesh(std::string filename);
Mesh readMPASMesh(int ncid);

}

#endif

