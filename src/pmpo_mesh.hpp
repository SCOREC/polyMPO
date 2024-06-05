#ifndef POLYMPO_MESH_H
#define POLYMPO_MESH_H

#include "pmpo_utils.hpp"

namespace polyMPO{

#define maxVtxsPerElm 8
#define maxElmsPerVtx 5

using IntVtx2ElmView = Kokkos::View<int*[maxVtxsPerElm+1]>;
using IntElm2VtxView = Kokkos::View<int*[maxElmsPerVtx+1]>;
using IntElm2ElmView = Kokkos::View<int*[maxVtxsPerElm+1]>;

enum MeshFieldIndex{
    MeshF_Invalid = -2,
    MeshF_Unsupported,
    MeshF_VtxCoords,
    MeshF_VtxRotLat,
    MeshF_Vel,
    MeshF_VtxMass,
    MeshF_ElmMass,
    MeshF_OnSurfVeloIncr,
    MeshF_OnSurfDispIncr,
    MeshF_RotLatLonIncr
};
enum MeshFieldType{
    MeshFType_Invalid = -2,
    MeshFType_Unsupported,
    MeshFType_VtxBased,
    MeshFType_ElmBased
};

template <MeshFieldIndex> struct meshFieldToType;
template <> struct meshFieldToType < MeshF_VtxCoords         > { using type = Kokkos::View<vec3d_t*>; };
template <> struct meshFieldToType < MeshF_VtxRotLat         > { using type = DoubleView; };
template <> struct meshFieldToType < MeshF_Vel               > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_VtxMass           > { using type = Kokkos::View<doubleSclr_t*>; };
template <> struct meshFieldToType < MeshF_ElmMass           > { using type = Kokkos::View<doubleSclr_t*>; };
template <> struct meshFieldToType < MeshF_OnSurfVeloIncr    > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_OnSurfDispIncr    > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_RotLatLonIncr     > { using type = Kokkos::View<vec2d_t*>; };

template <MeshFieldIndex index>
using MeshFView = typename meshFieldToType<index>::type;

const std::map<MeshFieldIndex, std::pair<MeshFieldType,
                                         std::string>> meshFields2TypeAndString = 
              {{MeshF_Invalid,          {MeshFType_Invalid,"MeshField_InValid!"}},
               {MeshF_Unsupported,      {MeshFType_Unsupported,"MeshField_Unsupported"}},
               {MeshF_VtxCoords,        {MeshFType_VtxBased,"MeshField_VerticesCoords"}},
               {MeshF_VtxRotLat,        {MeshFType_VtxBased,"MeshField_VerticesLatitude"}},
               {MeshF_Vel,              {MeshFType_VtxBased,"MeshField_Velocity"}},
               {MeshF_VtxMass,          {MeshFType_VtxBased,"MeshField_VerticesMass"}},
               {MeshF_ElmMass,          {MeshFType_ElmBased,"MeshField_ElementsMass"}},
               {MeshF_OnSurfVeloIncr,   {MeshFType_VtxBased,"MeshField_OnSurfaceVelocityIncrement"}},
               {MeshF_OnSurfDispIncr,   {MeshFType_VtxBased,"MeshField_OnSurfaceDisplacementIncrement"}},
               {MeshF_RotLatLonIncr,    {MeshFType_VtxBased,"MeshField_RotationalLatitudeLongitudeIncreasement"}}};

enum mesh_type {mesh_unrecognized_lower = -1,
                mesh_general_polygonal, //other meshes
                mesh_CVT_polygonal,     //MPAS meshes
                mesh_unrecognized_upper};

enum geom_type {geom_unrecognized_lower = -1,
                geom_planar_surf,
                geom_spherical_surf,
                geom_unrecognized_upper};

class Mesh {
  private:
    bool meshEdit_ = false;
    mesh_type meshType_ = mesh_unrecognized_lower;
    geom_type geomType_ = geom_unrecognized_lower;

    double sphereRadius_;
    int numVtxs_;
    int numElms_;
    //IntView nEdgesPerElm_;
    IntVtx2ElmView elm2VtxConn_;
    IntElm2ElmView elm2ElmConn_;
  
    //start of meshFields
    MeshFView<MeshF_VtxCoords> vtxCoords_;
    MeshFView<MeshF_VtxRotLat> vtxRotLat_;
    MeshFView<MeshF_Vel> vtxVel_;
    MeshFView<MeshF_VtxMass> vtxMass_;
    MeshFView<MeshF_ElmMass> elmMass_;
    MeshFView<MeshF_OnSurfVeloIncr> vtxOnSurfVeloIncr_;
    MeshFView<MeshF_OnSurfDispIncr> vtxOnSurfDispIncr_;
    MeshFView<MeshF_RotLatLonIncr> vtxRotLatLonIncr_;
    //DoubleMat2DView vtxStress_;

  public:
    Mesh(){};
    Mesh( mesh_type meshType,
          geom_type geomType,
          double sphereRadius,
          int numVtxs,
          int numElms,
          MeshFView<MeshF_VtxCoords> vtxCoords,
          IntVtx2ElmView elm2VtxConn,
          IntElm2ElmView elm2ElmConn ):
          meshType_(meshType),
          geomType_(geomType),
          sphereRadius_(sphereRadius),
          numVtxs_(numVtxs),
          numElms_(numElms),
          elm2VtxConn_(elm2VtxConn),
          elm2ElmConn_(elm2ElmConn){
            meshEdit_ = true;
            setMeshVtxBasedFieldSize();
            setMeshElmBasedFieldSize();
            meshEdit_ = false;
            vtxCoords_ = vtxCoords;
        }

    bool meshEditable(){ return meshEdit_; }
    bool checkMeshType(int meshType);
    bool checkGeomType(int geomType);

    mesh_type getMeshType() { return meshType_; }
    geom_type getGeomType() { return geomType_; }
    double getSphereRadius() { return sphereRadius_; }
    int getNumVertices() { return numVtxs_; }
    int getNumElements() { return numElms_; }
    IntVtx2ElmView getElm2VtxConn() { return elm2VtxConn_; }
    IntElm2ElmView getElm2ElmConn() { return elm2ElmConn_; }
    template<MeshFieldIndex index> auto getMeshField();
    void setMeshVtxBasedFieldSize();
    void setMeshElmBasedFieldSize();

    void setMeshEdit(bool meshEdit) { meshEdit_ = meshEdit; }
    //onec MeshType/GeomType is set to valid types, we can't change them anymore
    void setMeshType(mesh_type meshType) {PMT_ALWAYS_ASSERT(!checkMeshType(meshType_) && meshEdit_);
                                          meshType_ = meshType;}
    void setGeomType(geom_type geomType) {PMT_ALWAYS_ASSERT(!checkGeomType(geomType_) && meshEdit_);
                                          geomType_ = geomType;}
    void setSphereRadius(double sphereRadius) {PMT_ALWAYS_ASSERT(meshEdit_);
                                               sphereRadius_ = sphereRadius;}
    void setNumVtxs(int numVtxs) {PMT_ALWAYS_ASSERT(meshEdit_);
                                  numVtxs_ = numVtxs;}
    void setNumElms(int numElms) {PMT_ALWAYS_ASSERT(meshEdit_);
                                  numElms_ = numElms;}
    void setElm2VtxConn(IntVtx2ElmView elm2VtxConn) {PMT_ALWAYS_ASSERT(meshEdit_);
                                                     elm2VtxConn_ = elm2VtxConn; }
    void setElm2ElmConn(IntElm2ElmView elm2ElmConn) {PMT_ALWAYS_ASSERT(meshEdit_);
                                                     elm2ElmConn_ = elm2ElmConn; }
    
    void computeRotLatLonIncr();
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
    else if constexpr (index==MeshF_VtxCoords){
        return vtxCoords_;
    }
    else if constexpr (index==MeshF_VtxRotLat){
        return vtxRotLat_;
    }
    else if constexpr (index==MeshF_Vel){
        return vtxVel_;
    }
    else if constexpr (index==MeshF_VtxMass){
        return vtxMass_;
    }
    else if constexpr (index==MeshF_ElmMass){
        return elmMass_;
    }
    else if constexpr (index==MeshF_OnSurfVeloIncr){
        return vtxOnSurfVeloIncr_;
    }
    else if constexpr (index==MeshF_OnSurfDispIncr){
        return vtxOnSurfDispIncr_;
    }
    else if constexpr (index==MeshF_RotLatLonIncr){
        return vtxRotLatLonIncr_;
    }
    fprintf(stderr,"Mesh Field Index error!\n");
    exit(1);
}

}

#endif
