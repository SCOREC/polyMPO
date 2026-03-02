#ifndef POLYMPO_MESH_H
#define POLYMPO_MESH_H

#include "pmpo_utils.hpp"
#include <pumipic_kktypes.hpp>
#include <particle_structs.hpp>

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
    MeshF_ElmCenterXYZ,
    MeshF_DualTriangleArea,
    MeshF_Vel,
    MeshF_VtxMass,
    MeshF_ElmMass,
    MeshF_OnSurfVeloIncr,
    MeshF_OnSurfDispIncr,
    MeshF_RotLatLonIncr,
    MeshF_VtxGnomProj,
    MeshF_ElmCenterGnomProj,
    MeshF_TanLatVertexRotatedOverRadius,
    MeshF_SolveStress,
    MeshF_SolveVelocity,
    MeshF_InteriorVertex,
    MeshF_StressDivergence,
    MeshF_TotalMassVtx,
    MeshF_AirStress,
    MeshF_SurfaceTilt,
    MeshF_TotalMassFVtx,
    MeshF_OceanStress,
    MeshF_OceanStressCoeff
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
template <> struct meshFieldToType < MeshF_ElmCenterXYZ      > { using type = Kokkos::View<vec3d_t*>; };
template <> struct meshFieldToType < MeshF_DualTriangleArea  > { using type = Kokkos::View<doubleSclr_t*>; };
template <> struct meshFieldToType < MeshF_Vel               > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_VtxMass           > { using type = Kokkos::View<doubleSclr_t*>; };
template <> struct meshFieldToType < MeshF_ElmMass           > { using type = Kokkos::View<doubleSclr_t*>; };
template <> struct meshFieldToType < MeshF_OnSurfVeloIncr    > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_OnSurfDispIncr    > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_RotLatLonIncr     > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_VtxGnomProj       > { using type = Kokkos::View<double*[maxVtxsPerElm][2]>; };
template <> struct meshFieldToType < MeshF_ElmCenterGnomProj > { using type = Kokkos::View<double*[4]>; };
template <> struct meshFieldToType < MeshF_TanLatVertexRotatedOverRadius > { using type = Kokkos::View<doubleSclr_t*>; };
template <> struct meshFieldToType < MeshF_SolveStress       > { using type = IntView; };
template <> struct meshFieldToType < MeshF_SolveVelocity     > { using type = IntView; };
template <> struct meshFieldToType < MeshF_InteriorVertex    > { using type = IntView; };
template <> struct meshFieldToType < MeshF_StressDivergence  > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_TotalMassVtx      > { using type = Kokkos::View<doubleSclr_t*>; };
template <> struct meshFieldToType < MeshF_AirStress         > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_SurfaceTilt       > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_TotalMassFVtx     > { using type = Kokkos::View<doubleSclr_t*>; };
template <> struct meshFieldToType < MeshF_OceanStress       > { using type = Kokkos::View<vec2d_t*>; };
template <> struct meshFieldToType < MeshF_OceanStressCoeff  > { using type = Kokkos::View<doubleSclr_t*>; };

template <MeshFieldIndex index>
using MeshFView = typename meshFieldToType<index>::type;

const std::map<MeshFieldIndex, std::pair<MeshFieldType, std::string>> meshFields2TypeAndString = {
        {MeshF_Invalid,          {MeshFType_Invalid,"MeshField_InValid!"}},
        {MeshF_Unsupported,      {MeshFType_Unsupported,"MeshField_Unsupported"}},
        {MeshF_VtxCoords,        {MeshFType_VtxBased,"MeshField_VerticesCoords"}},
        {MeshF_VtxRotLat,        {MeshFType_VtxBased,"MeshField_VerticesLatitude"}},
        {MeshF_ElmCenterXYZ,     {MeshFType_ElmBased,"MeshField_ElementCenterXYZ"}},
        {MeshF_DualTriangleArea, {MeshFType_VtxBased,"MeshField_DualTriangleArea"}},
        {MeshF_Vel,              {MeshFType_VtxBased,"MeshField_Velocity"}},
        {MeshF_VtxMass,          {MeshFType_VtxBased,"MeshField_VerticesMass"}},
        {MeshF_ElmMass,          {MeshFType_ElmBased,"MeshField_ElementsMass"}},
        {MeshF_OnSurfVeloIncr,   {MeshFType_VtxBased,"MeshField_OnSurfaceVelocityIncrement"}},
        {MeshF_OnSurfDispIncr,   {MeshFType_VtxBased,"MeshField_OnSurfaceDisplacementIncrement"}},
        {MeshF_RotLatLonIncr,    {MeshFType_VtxBased,"MeshField_RotationalLatitudeLongitudeIncreasement"}},
        {MeshF_VtxGnomProj,      {MeshFType_ElmBased,"MeshField_VertexGnomonicProjection"}},
        {MeshF_ElmCenterGnomProj,{MeshFType_ElmBased,"MeshField_ElementCenterGnomonicprojection"}},
        {MeshF_TanLatVertexRotatedOverRadius, {MeshFType_VtxBased,"MeshField_TanLatVertexRotatedOverRadius"}},
        {MeshF_SolveStress,      {MeshFType_ElmBased,"MeshField_SolveStress"}},
        {MeshF_SolveVelocity,    {MeshFType_VtxBased,"MeshField_SolveVelocity"}},
        {MeshF_InteriorVertex,   {MeshFType_VtxBased,"MeshField_InteriorVertex"}},
        {MeshF_StressDivergence, {MeshFType_VtxBased,"MeshField_StressDivergence"}},
        {MeshF_TotalMassVtx,     {MeshFType_VtxBased,"MeshField_TotalMassVtx"}},
        {MeshF_AirStress,        {MeshFType_VtxBased,"MeshField_AirStress"}},
        {MeshF_SurfaceTilt,      {MeshFType_VtxBased,"MeshField_SurfaceTilt"}},
        {MeshF_TotalMassFVtx,    {MeshFType_VtxBased,"MeshField_TotalMassFVtx"}},
        {MeshF_OceanStress,      {MeshFType_VtxBased,"MeshField_OceanStress"}},
        {MeshF_OceanStressCoeff, {MeshFType_VtxBased,"MeshField_OceanStressCoeff"}}
};

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
    IntView owningProc_;
    IntView owningProcVertex_;
    IntView globalElm_;
    IntView globalVtx_;
    //start of meshFields
    MeshFView<MeshF_VtxCoords> vtxCoords_;
    MeshFView<MeshF_VtxRotLat> vtxRotLat_;
    MeshFView<MeshF_ElmCenterXYZ> elmCenterXYZ_;
    MeshFView<MeshF_DualTriangleArea> dualTriangleArea_;
    MeshFView<MeshF_StressDivergence> stressDivergence_;
    MeshFView<MeshF_Vel> vtxVel_;
    MeshFView<MeshF_VtxMass> vtxMass_;
    MeshFView<MeshF_ElmMass> elmMass_;
    MeshFView<MeshF_OnSurfVeloIncr> vtxOnSurfVeloIncr_;
    MeshFView<MeshF_OnSurfDispIncr> vtxOnSurfDispIncr_;
    MeshFView<MeshF_RotLatLonIncr> vtxRotLatLonIncr_;
    //GnomonicProjection
    MeshFView<MeshF_VtxGnomProj> vtxGnomProj_;
    MeshFView<MeshF_ElmCenterGnomProj> elmCenterGnomProj_;
    
    MeshFView<MeshF_TanLatVertexRotatedOverRadius> tanLatVertexRotatedOverRadius_;
    MeshFView<MeshF_SolveStress> solveStress_;
    MeshFView<MeshF_SolveVelocity> solveVelocity_;
    MeshFView<MeshF_InteriorVertex> interiorVertex_;
    MeshFView<MeshF_TotalMassVtx> totalMassVtx_;
    MeshFView<MeshF_AirStress> airStress_;
    MeshFView<MeshF_SurfaceTilt> surfaceTilt_;
    MeshFView<MeshF_TotalMassFVtx> totalMassFVtx_;
    MeshFView<MeshF_OceanStress> oceanStress_;
    MeshFView<MeshF_OceanStressCoeff> oceanStressCoeff_;

    bool isRotatedFlag = false;
    double elasticTimeStep_;
    double dynamicTimeStep_;
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

    void setOwningProc(IntView owningProc){
      PMT_ALWAYS_ASSERT(meshEdit_);
      owningProc_ = owningProc;
    }
    void setOwningProcVertex(IntView owningProcVertex){
      owningProcVertex_ = owningProcVertex; 
    } 
    IntView getElm2Process() {return owningProc_;}
    IntView getVtx2Process() {return owningProcVertex_;}

    mesh_type getMeshType() { return meshType_; }
    geom_type getGeomType() { return geomType_; }
    double getSphereRadius() { return sphereRadius_; }
    int getNumVertices() { return numVtxs_; }
    int getNumElements() { return numElms_; }
    IntVtx2ElmView getElm2VtxConn() { return elm2VtxConn_; }
    IntElm2ElmView getElm2ElmConn() { return elm2ElmConn_; }
    template<MeshFieldIndex index> auto getMeshField();
    template<MeshFieldIndex index> void fillMeshField(int size, int numEntries, double val);
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


    void setElmGlobal(IntView globalElm) {globalElm_ = globalElm;}
    void setVtxGlobal(IntView globalVtx) {globalVtx_ = globalVtx;}
    IntView getElmGlobal() {return globalElm_;}
    IntView getVtxGlobal() {return globalVtx_;}

    void setGnomonicProjection(bool isRotated);

    void computeRotLatLonIncr();

    bool getRotatedFlag() {
      return isRotatedFlag;
    }
    void setRotatedFlag(bool flagSet) {
      isRotatedFlag = flagSet;
    }
    
    void setElasticTimeStep(double elasticTimeStep){
      elasticTimeStep_ = elasticTimeStep;
    }
    double getElasticTimeStep(){
      return elasticTimeStep_;
    }
   
    void setDynamicTimeStep(double dynamicTimeStep){
      dynamicTimeStep_ = dynamicTimeStep;
    }
    double getDynamicTimeStep(){
      return dynamicTimeStep_;
    }

    void gridSolveGPU();
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
    else if constexpr (index==MeshF_ElmCenterXYZ){
        return elmCenterXYZ_;
    }
    else if constexpr (index==MeshF_DualTriangleArea){
        return dualTriangleArea_;
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
    else if constexpr (index==MeshF_VtxGnomProj){
        return vtxGnomProj_;
    }
    else if constexpr (index==MeshF_ElmCenterGnomProj){
        return elmCenterGnomProj_;
    }
    else if constexpr (index==MeshF_TanLatVertexRotatedOverRadius){
        return tanLatVertexRotatedOverRadius_;
    }
    else if constexpr (index==MeshF_SolveStress){
        return solveStress_;
    }
    else if constexpr (index==MeshF_SolveVelocity){
        return solveVelocity_;
    }
    else if constexpr (index==MeshF_InteriorVertex){
        return interiorVertex_;
    }
    else if constexpr (index==MeshF_StressDivergence){
        return stressDivergence_;
    }
    else if constexpr (index==MeshF_TotalMassVtx){
        return totalMassVtx_;
    }
    else if constexpr (index==MeshF_AirStress){
        return airStress_;
    }
    else if constexpr (index==MeshF_SurfaceTilt){
        return surfaceTilt_;
    }
    else if constexpr (index==MeshF_TotalMassFVtx){
        return totalMassFVtx_;
    }
    else if constexpr (index==MeshF_OceanStress){
        return oceanStress_;
    }
    else if constexpr (index==MeshF_OceanStressCoeff){
        return oceanStressCoeff_;
    }
    fprintf(stderr,"Mesh Field Index error!\n");
    exit(1);
}

template<MeshFieldIndex index> 
void Mesh::fillMeshField(int size, int numEntries, double val){
    auto meshField = getMeshField<index>();
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},{size, numEntries});
    Kokkos::parallel_for("fill mesh field", policy, KOKKOS_LAMBDA(const int i, const int j){
        meshField(i, j) = val;
    });
}

KOKKOS_INLINE_FUNCTION
void computeGnomonicProjectionAtPoint(const Vec3d& Coord, 
     const Kokkos::View<double[4], Kokkos::LayoutStride, Kokkos::MemoryTraits<Kokkos::Unmanaged>>& gnomProjElmCenter_sub, 
     double& outX, double& outY){   
  const double iDen = 1.0 / (gnomProjElmCenter_sub(1) * gnomProjElmCenter_sub(3) * Coord[0] +
                             gnomProjElmCenter_sub(0) * gnomProjElmCenter_sub(3) * Coord[1] +
                             gnomProjElmCenter_sub(2) * Coord[2]);
  outX = iDen * (Coord[1] * gnomProjElmCenter_sub(1) -
                 Coord[0] * gnomProjElmCenter_sub(0));
  outY = iDen * (Coord[2] * gnomProjElmCenter_sub(3) - Coord[1] * gnomProjElmCenter_sub(2) * gnomProjElmCenter_sub(0) -
                 Coord[0] * gnomProjElmCenter_sub(1) * gnomProjElmCenter_sub(2));
}
}

#endif
