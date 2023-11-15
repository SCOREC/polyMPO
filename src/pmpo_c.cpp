#include "pmpo_createTestMPMesh.hpp"
#include "pmpo_c.h"
#include <stdio.h>

#define MP_DETACHED -1 //TODO: consider other ways later, like enum
#define MP_ACTIVE 1

namespace{
  std::vector<MPMesh_ptr> p_mpmeshes;////store the p_mpmeshes that is legal
    
  void checkMPMeshValid(MPMesh_ptr p_mpmesh){
    auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
    PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  }
}

void polympo_initialize_f() {
  int isMPIInit;
  MPI_Initialized(&isMPIInit);
  PMT_ALWAYS_ASSERT(isMPIInit);
  Kokkos::initialize();
}

void polympo_finalize_f() {
  Kokkos::finalize();
}

MPMesh_ptr polympo_createMPMesh_f(int testMeshOption, int testMPOption) {
  polyMPO::Mesh* p_mesh;
  if(testMeshOption){
    int replicateFactor = 1;
    p_mesh = polyMPO::initTestMesh(testMeshOption, replicateFactor);
  }else{
    p_mesh = new polyMPO::Mesh();
  }
  polyMPO::MaterialPoints* p_mps;
  if(testMPOption){
    PMT_ALWAYS_ASSERT(testMeshOption >= 1);
    p_mps = polyMPO::initTestMPs(p_mesh, testMPOption);
  }else{
    p_mps = new polyMPO::MaterialPoints();  
  }
  MPMesh_ptr p_mpMeshReturn = (MPMesh_ptr) new polyMPO::MPMesh(p_mesh, p_mps);
  p_mpmeshes.push_back(p_mpMeshReturn);
  return p_mpMeshReturn;
}

void polympo_deleteMPMesh_f(MPMesh_ptr p_mpmesh) {
  //check mpMesh is valid
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  p_mpmeshes.erase(p_mpmeshIter);
  delete (polyMPO::MPMesh*)p_mpmesh;
}

void polympo_setMPICommunicator_f(MPI_Fint fcomm){
    MPI_Comm comm = MPI_Comm_f2c(fcomm);
    int commSize;
    MPI_Comm_size(comm,&commSize);
}

/**
 * Attention: this typedef is LayoutLeft, meaning that the first
 * index is the contiguous one. This matches the Fortran and GPU conventions for
 * allocations.
 */
//TODO: order of these typedefs to be done later
typedef Kokkos::View<
          double*,
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>
        > kkDblViewHostU;//TODO:put it to mesh.hpp
                         
typedef Kokkos::View<
          double*[vec2d_nEntries],
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>
        > kkVec2dViewHostU;//TODO:put it to mesh.hpp

typedef Kokkos::View<
          double**,
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>
        > kkDbl2dViewHostU;//TODO:put it somewhere else (maybe)

typedef Kokkos::View<
          int**,
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>
        > kkInt2dViewHostU;//TODO:put it somewhere else (maybe)

typedef Kokkos::View<
          int*,
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>
        > kkIntViewHostU;//TODO:put it somewhere else (maybe)

void polympo_createMPs_f(MPMesh_ptr p_mpmesh,
                       int numElms,
                       int numMPs, // >= number of active MPs
                       int* mpsPerElm,
                       int* mp2Elm,
                       int* isMPActive) {
  checkMPMeshValid(p_mpmesh);

  //the mesh must be fixed/set before adding MPs
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  PMT_ALWAYS_ASSERT(!p_mesh->meshEditable());
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements() == numElms);

  int firstElmWithMPs=-1;
  for (int i=0; i<numElms; i++) {
    if(mpsPerElm[i]) {
      firstElmWithMPs = i;
      break;
    }
  }

  int minElmID = numElms+1;
  for(int i = 0; i < numMPs; i++) {
    if(isMPActive[i] == MP_ACTIVE) {
      if(mp2Elm[i] < minElmID) {
        minElmID = mp2Elm[i];
      }
    }
  }
  int offset = -1;
  if(minElmID-firstElmWithMPs==1) {
    offset = 1;
  }else if (minElmID-firstElmWithMPs==0){
    offset = 0;
  }else {
    fprintf(stderr,"The minElmID is incorrect! Offset is wrong!\n");
    exit(1);
  }

  std::vector<int> active_mpIDs(numMPs);
  std::vector<int> active_mp2Elm(numMPs);
  int numActiveMPs = 0;
  for(int i=0; i<numMPs; i++) {
    if(isMPActive[i] == MP_ACTIVE) {
      active_mpIDs[numActiveMPs] = i;
      active_mp2Elm[numActiveMPs] = mp2Elm[i]-offset; //adjust for 1 based indexing if needed
      numActiveMPs++;
    }
  }

  //TODO do we care about empty ranks? check just in case...
  PMT_ALWAYS_ASSERT(numActiveMPs>0);

  using space_t = Kokkos::DefaultExecutionSpace::memory_space;
  kkIntViewHostU mpsPerElm_h(mpsPerElm,numElms);
  auto mpsPerElm_d = Kokkos::create_mirror_view_and_copy(space_t(), mpsPerElm_h);

  kkIntViewHostU active_mp2Elm_h(active_mp2Elm.data(),numActiveMPs);
  auto active_mp2Elm_d = Kokkos::create_mirror_view_and_copy(space_t(), active_mp2Elm_h);

  kkIntViewHostU active_mpIDs_h(active_mpIDs.data(),numActiveMPs);
  auto active_mpIDs_d = Kokkos::create_mirror_view_and_copy(space_t(), active_mpIDs_h);

  delete ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  ((polyMPO::MPMesh*)p_mpmesh)->p_MPs =
     new polyMPO::MaterialPoints(numElms, numActiveMPs, mpsPerElm_d, active_mp2Elm_d, active_mpIDs_d);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  p_MPs->setElmIDoffset(offset);
}

void polympo_getMPCurElmID_f(MPMesh_ptr p_mpmesh,
                           int numMPs,
                           int* elmIDs){
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  auto mpCurElmID = p_MPs->getData<polyMPO::MPF_Cur_Elm_ID>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  auto elmIDoffset = p_MPs->getElmIDoffset();

  kkIntViewHostU arrayHost(elmIDs,numMPs);
  polyMPO::IntView mpCurElmIDCopy("mpCurElmIDNewValue",numMPs);

  auto getElmId = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
        mpCurElmIDCopy(mpAppID(mp)) = mpCurElmID(mp)+elmIDoffset;
    }
  };
  p_MPs->parallel_for(getElmId, "get mpCurElmID");
  Kokkos::deep_copy( arrayHost, mpCurElmIDCopy);
}

void polympo_setMPLatLonRotatedFlag_f(MPMesh_ptr p_mpmesh, int isRotateFlag){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->p_MPs->setRotatedFlag(isRotateFlag>0);

}

void polympo_setMPPositions_f(MPMesh_ptr p_mpmesh,
                            int numComps,
                            int numMPs,
                            double* mpPositionsIn){
  static int callCount = 0;
  PMT_ALWAYS_ASSERT(callCount == 0);
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numComps == vec3d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpPositions = p_MPs->getData<polyMPO::MPF_Cur_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  kkDbl2dViewHostU mpPositionsIn_h(mpPositionsIn,numComps,numMPs);
  Kokkos::View<double**> mpPositionsIn_d("mpPositionsDevice",vec3d_nEntries,numMPs);
  Kokkos::deep_copy(mpPositionsIn_d, mpPositionsIn_h);
  auto setPos = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpPositions(mp,0) = mpPositionsIn_d(0, mpAppID(mp));
      mpPositions(mp,1) = mpPositionsIn_d(1, mpAppID(mp));
      mpPositions(mp,2) = mpPositionsIn_d(2, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setPos, "setMPPositions");
  callCount++;
}

void polympo_getMPPositions_f(MPMesh_ptr p_mpmesh,
                            int numComps,
                            int numMPs,
                            double* mpPositionsHost){
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numComps == vec3d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpPositions = p_MPs->getData<polyMPO::MPF_Cur_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpPositionsCopy("mpPositionsCopy",vec3d_nEntries,numMPs);
  auto getPos = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpPositionsCopy(0,mpAppID(mp)) = mpPositions(mp,0);
      mpPositionsCopy(1,mpAppID(mp)) = mpPositions(mp,1);
      mpPositionsCopy(2,mpAppID(mp)) = mpPositions(mp,2);
    }
  };
  p_MPs->parallel_for(getPos, "getMPPositions");
  kkDbl2dViewHostU arrayHost(mpPositionsHost,numComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpPositionsCopy);
}

void polympo_setMPRotLatLon_f(MPMesh_ptr p_mpmesh,
                         int numComps,
                         int numMPs,
                         double* mpRotLatLonIn){
  static int callCount = 0;
  PMT_ALWAYS_ASSERT(callCount == 0);
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpRotLatLon = p_MPs->getData<polyMPO::MPF_Cur_Pos_Rot_Lat_Lon>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  kkDbl2dViewHostU mpRotLatLonIn_h(mpRotLatLonIn,numComps,numMPs);
  Kokkos::View<double**> mpRotLatLonIn_d("mpRotLatLonDevice",vec2d_nEntries,numMPs);
  Kokkos::deep_copy(mpRotLatLonIn_d, mpRotLatLonIn_h);
  auto setPos = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpRotLatLon(mp,0) = mpRotLatLonIn_d(0, mpAppID(mp));
      mpRotLatLon(mp,1) = mpRotLatLonIn_d(1, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setPos, "setMPRotLatLon");
  callCount++;
}

void polympo_getMPRotLatLon_f(MPMesh_ptr p_mpmesh,
                         int numComps,
                         int numMPs,
                         double* mpRotLatLonHost){
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpRotLatLon = p_MPs->getData<polyMPO::MPF_Cur_Pos_Rot_Lat_Lon>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpRotLatLonCopy("mpRotLatLonCopy",vec2d_nEntries,numMPs);
  auto getPos = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpRotLatLonCopy(0,mpAppID(mp)) = mpRotLatLon(mp,0);
      mpRotLatLonCopy(1,mpAppID(mp)) = mpRotLatLon(mp,1);
    }
  };
  p_MPs->parallel_for(getPos, "getMPRotLatLon");
  kkDbl2dViewHostU arrayHost(mpRotLatLonHost,numComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpRotLatLonCopy);
}

void polympo_setMPVel_f(MPMesh_ptr p_mpmesh, int size, double* array) {
  fprintf(stderr,"%s is no longer supported\n", __func__);
  PMT_ALWAYS_ASSERT(false);
  (void)p_mpmesh;// to silence the unused param warning
  (void)size;
  (void)array;
}

void polympo_getMPVel_f(MPMesh_ptr p_mpmesh, int size, double* array) {
  fprintf(stderr,"%s is no longer supported\n", __func__);
  PMT_ALWAYS_ASSERT(false);
  (void)p_mpmesh;// to silence the unused param warning
  (void)size;
  (void)array;
}

void polympo_startMeshFill_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->p_mesh->setMeshEdit(true);  
}

void polympo_endMeshFill_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh; 
  PMT_ALWAYS_ASSERT(p_mesh->meshEditable());
  p_mesh->setMeshEdit(false);  
}

void polympo_checkMeshMaxSettings_f(MPMesh_ptr p_mpmesh, int maxEdges, int vertexDegree){
  checkMPMeshValid(p_mpmesh);
  PMT_ALWAYS_ASSERT(maxEdges <= maxVtxsPerElm);
  PMT_ALWAYS_ASSERT(vertexDegree <=  maxElmsPerVtx);
}

void polympo_setMeshNumVtxs_f(MPMesh_ptr p_mpmesh, int numVtxs){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  p_mesh->setNumVtxs(numVtxs);
  p_mesh->setMeshVtxBasedFieldSize(); 
}

int polympo_getMeshNumVtxs_f(MPMesh_ptr p_mpmesh) {
  checkMPMeshValid(p_mpmesh); //chech vailidity
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  int nVtxs = p_mesh->getNumVertices();
  return nVtxs;
}

void polympo_setMeshNumElms_f(MPMesh_ptr p_mpmesh, int numElms){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  auto elm2Vtx = polyMPO::IntVtx2ElmView("MeshElementsToVertices",numElms); 
  auto elm2Elm = polyMPO::IntElm2ElmView("MeshElementsToElements",numElms); 

  p_mesh->setNumElms(numElms);
  p_mesh->setElm2VtxConn(elm2Vtx);
  p_mesh->setElm2ElmConn(elm2Elm);
}

int polympo_getMeshNumElms_f(MPMesh_ptr p_mpmesh) {
  checkMPMeshValid(p_mpmesh); //chech vailidity
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  int nElms = p_mesh->getNumElements();
  return nElms;
}

void polympo_setMeshTypeGeneralPoly_f(MPMesh_ptr p_mpmesh){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->p_mesh->setMeshType(polyMPO::mesh_general_polygonal);
}

void polympo_setMeshTypeCVTPoly_f(MPMesh_ptr p_mpmesh){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->p_mesh->setMeshType(polyMPO::mesh_CVT_polygonal);
}

void polympo_setMeshGeomTypePlanar_f(MPMesh_ptr p_mpmesh){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->p_mesh->setGeomType(polyMPO::geom_planar_surf);
}

void polympo_setMeshGeomTypeSpherical_f(MPMesh_ptr p_mpmesh){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->p_mesh->setGeomType(polyMPO::geom_spherical_surf);
}

void polympo_setMeshSphereRadius_f(MPMesh_ptr p_mpmesh, double sphereRadius){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  PMT_ALWAYS_ASSERT(sphereRadius >= 0);
  p_mesh->setSphereRadius(sphereRadius);
}

void polympo_setMeshNumEdgesPerElm_f(MPMesh_ptr p_mpmesh, int nCells, int* array){
  //chech vailidity
  checkMPMeshValid(p_mpmesh);
  kkIntViewHostU arrayHost(array,nCells);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  PMT_ALWAYS_ASSERT(p_mesh->meshEditable());

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==nCells);
  
  polyMPO::IntView nEdgesPerElm("MeshNumEdgesPerElm",nCells);
  Kokkos::deep_copy(nEdgesPerElm, arrayHost);
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto elm2ElmConn = p_mesh->getElm2ElmConn();
  Kokkos::parallel_for("set nEdgesPerElm", nCells, KOKKOS_LAMBDA(const int elm){
    elm2VtxConn(elm,0) = nEdgesPerElm(elm);
    elm2ElmConn(elm,0) = nEdgesPerElm(elm);
  });
}

void polympo_setMeshElm2VtxConn_f(MPMesh_ptr p_mpmesh, int maxEdges, int nCells, int* array){
  //chech vailidity
  checkMPMeshValid(p_mpmesh);
  kkInt2dViewHostU arrayHost(array,maxEdges,nCells); 
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh; 
  PMT_ALWAYS_ASSERT(p_mesh->meshEditable());

  //check the size
  PMT_ALWAYS_ASSERT(maxEdges <= maxVtxsPerElm);
  PMT_ALWAYS_ASSERT(nCells == p_mesh->getNumElements());
  
  Kokkos::View<int**> elm2VtxArray("MeshElementsToVertices",maxEdges,nCells);
  Kokkos::deep_copy(elm2VtxArray, arrayHost);
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  Kokkos::parallel_for("set elm2VtxConn", nCells, KOKKOS_LAMBDA(const int elm){
    for(int i=0; i<maxEdges; i++){
        elm2VtxConn(elm,i+1) = elm2VtxArray(i,elm);
    }
  });
}

void polympo_setMeshElm2ElmConn_f(MPMesh_ptr p_mpmesh, int maxEdges, int nCells, int* array){
  //chech vailidity
  checkMPMeshValid(p_mpmesh);
  kkInt2dViewHostU arrayHost(array,maxEdges,nCells); //Fortran is column-major
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh; 
  PMT_ALWAYS_ASSERT(p_mesh->meshEditable());

  //check the size
  PMT_ALWAYS_ASSERT(maxEdges <= maxVtxsPerElm);
  PMT_ALWAYS_ASSERT(nCells == p_mesh->getNumElements());
  
  Kokkos::View<int**> elm2ElmArray("MeshElementsToVertices",maxEdges,nCells);
  Kokkos::deep_copy(elm2ElmArray, arrayHost);
  auto elm2ElmConn = p_mesh->getElm2ElmConn();
  Kokkos::parallel_for("set elm2ElmConn", nCells, KOKKOS_LAMBDA(const int elm){
    for(int i=0; i<maxEdges; i++){
        elm2ElmConn(elm,i+1) = elm2ElmArray(i,elm);
    }  
  });
}

void polympo_setMeshVtxCoords_f(MPMesh_ptr p_mpmesh, int nVertices, double* xArray, double* yArray, double* zArray){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices); 

  //copy the host array to the device
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  auto h_coordsArray = Kokkos::create_mirror_view(coordsArray);
  for(int i=0; i<nVertices; i++){
    h_coordsArray(i,0) = xArray[i];
    h_coordsArray(i,1) = yArray[i];
    h_coordsArray(i,2) = zArray[i];
  }
  Kokkos::deep_copy(coordsArray, h_coordsArray);
}

void polympo_getMeshVtxCoords_f(MPMesh_ptr p_mpmesh, int nVertices, double* xArray, double* yArray, double* zArray){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices); 
  
  //copy the device to host 
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  auto h_coordsArray = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                           coordsArray);
  for(int i=0; i<nVertices; i++){
    xArray[i] = h_coordsArray(i,0);
    yArray[i] = h_coordsArray(i,1);
    zArray[i] = h_coordsArray(i,2);
  }
}

void polympo_setMeshVtxRotLatLon_f(MPMesh_ptr p_mpmesh, int nVertices, double* latitude, double* longitude){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices); 

  //copy the host array to the device
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_VtxRotLatLon>();
  auto h_coordsArray = Kokkos::create_mirror_view(coordsArray);
  for(int i=0; i<nVertices; i++){
    h_coordsArray(i,0) = latitude[i];
    h_coordsArray(i,1) = longitude[i];
  }
  Kokkos::deep_copy(coordsArray, h_coordsArray);
}

void polympo_getMeshVtxRotLatLon_f(MPMesh_ptr p_mpmesh, int nVertices, double* latitude, double* longitude){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices); 
  
  //copy the device to host 
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_VtxRotLatLon>();
  auto h_coordsArray = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                           coordsArray);
  for(int i=0; i<nVertices; i++){
    latitude[i] = h_coordsArray(i,0);
    longitude[i] = h_coordsArray(i,1);
  }
}

void polympo_setMeshOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, int nComps, int nVertices, double* array) {
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkVec2dViewHostU arrayHost(array,nVertices);

  auto vtxField = p_mesh->getMeshField<polyMPO::MeshF_OnSurfVeloIncr>();

  //check the size
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the host array to the device
  Kokkos::deep_copy(vtxField,arrayHost);
}

void polympo_getMeshOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, int nComps, int nVertices, double* array) {
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkVec2dViewHostU arrayHost(array,nVertices);

  auto vtxField = p_mesh->getMeshField<polyMPO::MeshF_OnSurfVeloIncr>();

  //check the size
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices() == nVertices); 
  PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the device array to the host
  Kokkos::deep_copy(arrayHost, vtxField);
}

void polympo_setMeshOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, int nComps, int nVertices, double* array) {
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkDbl2dViewHostU arrayHost(array,nComps,nVertices);
  Kokkos::View<double**> array_d("meshDispIncrDevice",nComps,nVertices);
  Kokkos::deep_copy(array_d, arrayHost);

  auto vtxField = p_mesh->getMeshField<polyMPO::MeshF_OnSurfDispIncr>();

  //check the size
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the host array to the device
  Kokkos::parallel_for("set mesh dispIncr", nVertices, KOKKOS_LAMBDA(const int iVtx){
    vtxField(iVtx,0) = array_d(0,iVtx);
    vtxField(iVtx,1) = array_d(1,iVtx);
  });
}

void polympo_getMeshOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, int nComps, int nVertices, double* array) {
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkDbl2dViewHostU arrayHost(array,nComps,nVertices);
  Kokkos::View<double**> array_d("meshDispIncrDevice",nComps,nVertices);

  auto vtxField = p_mesh->getMeshField<polyMPO::MeshF_OnSurfDispIncr>();

  //check the size
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices() == nVertices); 
  PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the device array to the host
  Kokkos::parallel_for("get mesh dispIncr", nVertices, KOKKOS_LAMBDA(const int iVtx){
    array_d(0,iVtx) = vtxField(iVtx,0);
    array_d(1,iVtx) = vtxField(iVtx,1);
  });
  Kokkos::deep_copy(arrayHost, array_d);
}

void polympo_push_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh) ->push();
}
