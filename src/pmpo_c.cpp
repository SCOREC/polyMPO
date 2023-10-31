#include "pmpo_createTestMPMesh.hpp"
#include "pmpo_defines.h"
#include "pmpo_c.h"
#include <stdio.h>

using space_t = Kokkos::DefaultExecutionSpace::memory_space;

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

MPMesh_ptr polympo_createMPMesh_f(const int testMeshOption, const int testMPOption) {
  polyMPO::Mesh* p_mesh;
  if(testMeshOption){
    int replicateFactor = 1;
    p_mesh = polyMPO::initTestMesh(testMeshOption, replicateFactor);
  }else{
    p_mesh = new polyMPO::Mesh();
  }
  polyMPO::MaterialPoints* p_mps;
  if(testMPOption){
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
template<typename DataT>
using kkViewHostU = Kokkos::View<
          DataT,
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

typedef kkViewHostU<double*> kkDblViewHostU;//TODO:put it to mesh.hpp             
typedef kkViewHostU<double*[vec2d_nEntries]> kkVec2dViewHostU;//TODO:put it to mesh.hpp
typedef kkViewHostU<double**> kkDbl2dViewHostU;//TODO:put it somewhere else (maybe)
typedef kkViewHostU<int**> kkInt2dViewHostU;//TODO:put it somewhere else (maybe)
typedef kkViewHostU<int*> kkIntViewHostU;//TODO:put it somewhere else (maybe)

template <typename DataT>
auto create_mirror_view_and_copy(DataT array, const int size){
  kkViewHostU<DataT> temp_host(array, size);
  return Kokkos::create_mirror_view_and_copy(space_t(), temp_host);
}

void polympo_createMPs_f(MPMesh_ptr p_mpmesh,
                       const int numElms,
                       const int numMPs, // total number of MPs which is GREATER than or equal to number of active MPs
                       const int* mpsPerElm,
                       const int* mp2Elm,
                       const int* isMPActive) {
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
      active_mpIDs[numActiveMPs] = i; //creates unique IDs
      active_mp2Elm[numActiveMPs] = mp2Elm[i]-offset; //adjust for 1 based indexing if needed
      numActiveMPs++;
    }
  }

  //TODO do we care about empty ranks? check just in case...
  PMT_ALWAYS_ASSERT(numActiveMPs>0);

  auto mpsPerElm_d = create_mirror_view_and_copy(mpsPerElm, numElms);
  auto active_mp2Elm_d = create_mirror_view_and_copy(active_mp2Elm.data(), numActiveMPs);
  auto active_mpIDs_d = create_mirror_view_and_copy(active_mpIDs.data(), numActiveMPs);

  delete ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  ((polyMPO::MPMesh*)p_mpmesh)->p_MPs =
     new polyMPO::MaterialPoints(numElms, numActiveMPs, mpsPerElm_d, active_mp2Elm_d, active_mpIDs_d);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  p_MPs->setElmIDoffset(offset);
}

void polympo_rebuildMPs_f(MPMesh_ptr p_mpmesh,
                         const int numMPs, // total number of MPs which is GREATER than or equal to number of active MPs
                         const int* allMP2Elm,
                         const int* addedMPMask) {
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  int offset = p_MPs->getElmIDoffset();
  std::vector<int> added_mpIDs(numMPs);
  std::vector<int> added_mp2Elm(numMPs);
  int numAddedMPs = 0;
  for(int i=0; i<numMPs; i++) {
    if(addedMPMask[i] == MP_ACTIVE) {
      added_mpIDs[numAddedMPs] = i;
      added_mp2Elm[numAddedMPs] = allMP2Elm[i]-offset; //adjust for 1 based indexing if needed
      numAddedMPs++;
    }
  }

  int internalMPCapacity = p_MPs->getCapacity(); // pumipic expects full capacity to rebuild
  Kokkos::View<int*> mp2Elm("mp2Elm", internalMPCapacity);
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();

  auto added_mp2Elm_d = create_mirror_view_and_copy(added_mp2Elm.data(), numAddedMPs);
  auto added_mpIDs_d = create_mirror_view_and_copy(added_mpIDs.data(), numAddedMPs);
  auto addedMPMask_d = create_mirror_view_and_copy(addedMPMask, numMPs);
  auto mpMP2ElmIn_d = create_mirror_view_and_copy(allMP2Elm, numMPs);

  auto setMP2Elm = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) {
      if (addedMPMask_d[mpAppID(mp)] == MP_ACTIVE) //two MPs can not occupy the same slot
        mp2Elm(mp) = MP_DELETE;
      else
        mp2Elm(mp) = mpMP2ElmIn_d(mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setMP2Elm, "setMP2Elm");
  p_MPs->rebuild(mp2Elm, numAddedMPs, added_mp2Elm_d, added_mpIDs_d);

  // check mpAppID is unique (on GPUs)
  if (p_MPs->getOpMode() == polyMPO::MP_DEBUG){
    mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
    Kokkos::View<int*> mpAppIDCount("mpAppIDCount", p_MPs->getCount());
    auto checkAppIDs = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) {
        int prev = Kokkos::atomic_fetch_add(&mpAppIDCount(mpAppID(mp)), 1);
        assert(prev == 0);
      }
    };
    p_MPs->parallel_for(checkAppIDs, "checkAppIDs");
  }
}

void polympo_getMPCurElmID_f(MPMesh_ptr p_mpmesh,
                           const int numMPs,
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

void polympo_setMPPositions_f(MPMesh_ptr p_mpmesh,
                            const int numComps,
                            const int numMPs,
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
                            const int numComps,
                            const int numMPs,
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

void polympo_setMPVel_f(MPMesh_ptr p_mpmesh, const int size, double* array) {
  fprintf(stderr,"%s is no longer supported\n", __func__);
  PMT_ALWAYS_ASSERT(false);
  (void)p_mpmesh;// to silence the unused param warning
  (void)size;
  (void)array;
}

void polympo_getMPVel_f(MPMesh_ptr p_mpmesh, const int size, double* array) {
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

void polympo_checkMeshMaxSettings_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int vertexDegree){
  checkMPMeshValid(p_mpmesh);
  PMT_ALWAYS_ASSERT(maxEdges <= maxVtxsPerElm);
  PMT_ALWAYS_ASSERT(vertexDegree <=  maxElmsPerVtx);
}

void polympo_setMeshNumVtxs_f(MPMesh_ptr p_mpmesh, const int numVtxs){
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

void polympo_setMeshNumElms_f(MPMesh_ptr p_mpmesh, const int numElms){
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

void polympo_setMeshSphereRadius_f(MPMesh_ptr p_mpmesh, const double sphereRadius){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  PMT_ALWAYS_ASSERT(sphereRadius >= 0);
  p_mesh->setSphereRadius(sphereRadius);
}

void polympo_setMeshNumEdgesPerElm_f(MPMesh_ptr p_mpmesh, const int nCells, int* array){
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

void polympo_setMeshElm2VtxConn_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int nCells, int* array){
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

void polympo_setMeshElm2ElmConn_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int nCells, int* array){
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

void polympo_setMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, double* xArray, double* yArray, double* zArray){
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

void polympo_getMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, double* xArray, double* yArray, double* zArray){
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

void polympo_setMeshOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array) {
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

void polympo_getMeshOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array) {
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

void polympo_setMeshOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array) {
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkVec2dViewHostU arrayHost(array,nVertices);

  auto vtxField = p_mesh->getMeshField<polyMPO::MeshF_OnSurfDispIncr>();

  //check the size
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the host array to the device
  Kokkos::deep_copy(vtxField,arrayHost);
}

void polympo_getMeshOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array) {
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkVec2dViewHostU arrayHost(array,nVertices);

  auto vtxField = p_mesh->getMeshField<polyMPO::MeshF_OnSurfDispIncr>();

  //check the size
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices() == nVertices); 
  PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the device array to the host

  Kokkos::deep_copy(arrayHost, vtxField);
}

