#ifndef POLYMPMTEST_MESH_H
#define POLYMPMTEST_MESH_H

#include "pmt_utils.hpp"

namespace polyMpmTest{

#define maxVtxsPerElm 8
#define maxElmsPerVtx 5

using IntVtx2ElmView = Kokkos::View<int*[maxVtxsPerElm+1]>;
using IntElm2VtxView = Kokkos::View<int*[maxElmsPerVtx+1]>;

using DoubleVec2DView = Kokkos::View<double*[2]>;
using DoubleVec3DView = Kokkos::View<double*[3]>;
using DoubleSymMat3DView = Kokkos::View<double*[6]>;

enum MeshFieldIndex{
    MeshF_Invalid = -2,
    MeshF_Unsupported,
    MeshF_Vel,
    MeshF_Cur_Pos_XYZ
};

class Mesh {
  private:
    int numVtxs_;
    int numElms_;
    Vector2View vtxCoords_;
    IntVtx2ElmView elm2VtxConn_;
    IntElm2VtxView vtx2ElmConn_;
    //start of meshFields
    DoubleVec2DView vtxVel_; 
    DoubleVec3DView vtxCurXYZ_;
    //DoubleMat2DView vtxStress_;

    static Mesh readMPASMesh(int ncid);

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
          vtx2ElmConn_(vtx2ElmConn){
            vtxVel_ = DoubleVec2DView("vtxVelocity",numVtxs);    
            vtxCurXYZ_ = DoubleVec3DView("vtxCurrentPositionXYZ",numVtxs);    
        }
    static Mesh readMPASMesh(std::string filename);

    int getNumVertices() { return numVtxs_; }
    int getNumElements() { return numElms_; }
    Vector2View getVtxCoords() { return vtxCoords_; }
    IntVtx2ElmView getElm2VtxConn() { return elm2VtxConn_; }
    IntElm2VtxView getVtx2ElmConn() { return vtx2ElmConn_; }
    template<MeshFieldIndex index> auto getMeshField();

    //void setVtx2ElmConn(IntElm2VtxView vtx2ElmConn) { vtx2ElmConn_ = vtx2ElmConn; }
    //void setElm2VtxConn(IntVtx2ElmView elm2VtxConn) { elm2VtxConn_ = elm2VtxConn; }
    //void setAssemblyReturn(DoubleView asmReturn) { Kokkos::resize(assemblyReturn_, numVtxs_);  Kokkos::deep_copy(assemblyReturn_,asmReturn); }
    
};

template<MeshFieldIndex index>
auto Mesh::getMeshField(){
    if constexpr (index==MeshF_Invalid){
        fprintf(stderr,"Mesh Field Invalid!\n");
        exit(1);
    }
    else if constexpr (index==MeshF_Unsupported){
        fprintf(stderr,"Mesh Field Unsupported!\n");
        exit(1);
    }
    else if constexpr (index==MeshF_Vel){
        return vtxVel_;                  
    }
    else if constexpr (index==MeshF_Cur_Pos_XYZ){
        return vtxCurXYZ_;
    }
    fprintf(stderr,"Mesh Field Index error!\n");
    exit(1);
}


}

#endif
