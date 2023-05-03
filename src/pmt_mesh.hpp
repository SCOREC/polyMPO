#ifndef POLYMPMTEST_MESH_H
#define POLYMPMTEST_MESH_H

#include "pmt_utils.hpp"

namespace polyMpmTest{

#define maxVtxsPerElm 8
#define maxElmsPerVtx 5

using IntVtx2ElmView = Kokkos::View<int*[maxVtxsPerElm+1]>;
using IntElm2VtxView = Kokkos::View<int*[maxElmsPerVtx+1]>;
using IntElm2ElmView = Kokkos::View<int*[maxVtxsPerElm*2+1]>;

class Mesh {
  private:
    int numVtxs_;
    int numElms_;
    Vector2View vtxCoords_;
    IntVtx2ElmView elm2VtxConn_;
    IntElm2VtxView vtx2ElmConn_;
    IntElm2ElmView elm2ElmConn_;    

    static Mesh readMPASMesh(int ncid);

  public:
    Mesh(){};
    Mesh( int numVtxs,
          int numElms,
          Vector2View vtxCoords,
          IntVtx2ElmView elm2VtxConn,
          IntElm2VtxView vtx2ElmConn,
          IntElm2ElmView elm2ElmConn):
          numVtxs_(numVtxs),
          numElms_(numElms),
          vtxCoords_(vtxCoords),
          elm2VtxConn_(elm2VtxConn),
          vtx2ElmConn_(vtx2ElmConn),
          elm2ElmConn_(elm2ElmConn){};
    static Mesh readMPASMesh(std::string filename);

    int getNumVertices() { return numVtxs_; }
    int getNumElements() { return numElms_; }
    Vector2View getVtxCoords() { return vtxCoords_; }
    IntVtx2ElmView getElm2VtxConn() { return elm2VtxConn_; }
    IntElm2VtxView getVtx2ElmConn() { return vtx2ElmConn_; }
    IntElm2ElmView getElm2ElmConn() { return elm2ElmConn_; }

    void setVtx2ElmConn(IntElm2VtxView vtx2ElmConn) { vtx2ElmConn_ = vtx2ElmConn; }
    void setElm2VtxConn(IntVtx2ElmView elm2VtxConn) { elm2VtxConn_ = elm2VtxConn; }
};

KOKKOS_INLINE_FUNCTION
Vector2 calcElmCenter(const int iElm, IntVtx2ElmView elm2VtxConn, Vector2View vtxCoords){
    int numVtx = elm2VtxConn(iElm,0);
    double sum_x = 0.0, sum_y = 0.0;
    for(int i=1; i<= numVtx; i++){
        sum_x += vtxCoords(elm2VtxConn(iElm,i)-1)[0];
        sum_y += vtxCoords(elm2VtxConn(iElm,i)-1)[1];
    }
    return Vector2(sum_x/numVtx, sum_y/numVtx);
}

}//end namespace polyMpmTest

#endif
