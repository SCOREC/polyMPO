#include "pmpo_createTestMPMesh.hpp"
#include "pmpo_defines.h"
#include "pmpo_c.h"
#include "pmpo_MPMesh_assembly.hpp"
#include <stdio.h>

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

void polympo_setMPICommunicator_f(MPMesh_ptr p_mpmesh, MPI_Fint fcomm){
  checkMPMeshValid(p_mpmesh);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  p_MPs->setMPIComm(comm);
}

void polympo_createMPs_f(MPMesh_ptr p_mpmesh,
                       const int numElms,
                       const int numMPs, // total number of MPs which is GREATER than or equal to number of active MPs
                       int* mpsPerElm,
                       const int* mp2Elm,
                       const int* isMPActive) {
  checkMPMeshValid(p_mpmesh);

  //the mesh must be fixed/set before adding MPs
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  PMT_ALWAYS_ASSERT(!p_mesh->meshEditable());
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements() == numElms);

  int numActiveMPs = 0;
  int minElmID = numElms+1;
  for(int i = 0; i < numMPs; i++) {
    if(isMPActive[i] == MP_ACTIVE) {
      if(mp2Elm[i] < minElmID) {
        minElmID = mp2Elm[i];
        numActiveMPs++;
      }
    }
  }
  //TODO do we care about empty ranks? check just in case...
  PMT_ALWAYS_ASSERT(numActiveMPs>0);

  int firstElmWithMPs=-1;
  for (int i=0; i<numElms; i++) {
    if(mpsPerElm[i]) {
      firstElmWithMPs = i;
      break;
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
  numActiveMPs = 0;
  for(int i=0; i<numMPs; i++) {
    if(isMPActive[i] == MP_ACTIVE) {
      active_mpIDs[numActiveMPs] = i; //creates unique IDs
      active_mp2Elm[numActiveMPs] = mp2Elm[i]-offset; //adjust for 1 based indexing if needed
      numActiveMPs++;
    }
  }

  auto mpsPerElm_d = create_mirror_view_and_copy(mpsPerElm, numElms);
  auto active_mp2Elm_d = create_mirror_view_and_copy(active_mp2Elm.data(), numActiveMPs);
  auto active_mpIDs_d = create_mirror_view_and_copy(active_mpIDs.data(), numActiveMPs);

  delete ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  ((polyMPO::MPMesh*)p_mpmesh)->p_MPs =
     new polyMPO::MaterialPoints(numElms, numActiveMPs, mpsPerElm_d, active_mp2Elm_d, active_mpIDs_d);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  p_MPs->setElmIDoffset(offset);
}

void polympo_startRebuildMPs_f(MPMesh_ptr p_mpmesh,
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

  Kokkos::View<int*> numDeletedMPs_d("numDeletedMPs", 1);
  auto setMP2Elm = PS_LAMBDA(const int&, const int& mp, const int& mask) {
    if(mask) {
      if (addedMPMask_d[mpAppID(mp)] == MP_ACTIVE) //two MPs can not occupy the same slot
        mp2Elm(mp) = MP_DELETE;
      else
        mp2Elm(mp) = mpMP2ElmIn_d(mpAppID(mp));
      if (mp2Elm(mp) == MP_DELETE)
        Kokkos::atomic_increment(&numDeletedMPs_d(0));
    }
  };
  p_MPs->parallel_for(setMP2Elm, "setMP2Elm");

  int numDeletedMPs = pumipic::getLastValue(numDeletedMPs_d);
  PMT_ALWAYS_ASSERT(numAddedMPs > 0 || numDeletedMPs > 0);

  p_MPs->startRebuild(mp2Elm, numAddedMPs, added_mp2Elm_d, added_mpIDs_d, addedMPMask_d);

  // check mpAppID is unique (on GPUs)
  if (p_MPs->getOpMode() == polyMPO::MP_DEBUG){
    mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
    Kokkos::View<int*> mpAppIDCount("mpAppIDCount", p_MPs->getCount());
    auto checkAppIDs = PS_LAMBDA(const int&, const int& mp, const int& mask){
      if(mask) {
        int prev = Kokkos::atomic_fetch_add(&mpAppIDCount(mpAppID(mp)), 1);
        assert(prev == 0);
      }
    };
    p_MPs->parallel_for(checkAppIDs, "checkAppIDs");
  }
}

void polympo_finishRebuildMPs_f(MPMesh_ptr p_mpmesh)
{
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  p_MPs->finishRebuild();
}

void polympo_setAppIDFunc_f(MPMesh_ptr p_mpmesh, IntVoidFunc getNext, void* appIDs) {
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  polyMPO::IntFunc getNextAppID = [getNext, appIDs]() { return getNext(appIDs); };
  p_MPs->setAppIDFunc(getNextAppID);
}

void polympo_getMPCurElmID_f(MPMesh_ptr p_mpmesh,
                           const int numMPs,
                           int* elmIDs){
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());
  auto mpCurElmID = p_MPs->getData<polyMPO::MPF_Cur_Elm_ID>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  auto elmIDoffset = p_MPs->getElmIDoffset();

  kkIntViewHostU arrayHost(elmIDs,numMPs);
  polyMPO::IntView mpCurElmIDCopy("mpCurElmIDNewValue",numMPs);

  auto getElmId = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask){
        mpCurElmIDCopy(mpAppID(mp)) = mpCurElmID(mp)+elmIDoffset;
    }
  };
  p_MPs->parallel_for(getElmId, "get mpCurElmID");
  Kokkos::deep_copy( arrayHost, mpCurElmIDCopy);
}

void polympo_setMPLatLonRotatedFlag_f(MPMesh_ptr p_mpmesh, const int isRotateFlag){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->p_MPs->setRotatedFlag(isRotateFlag>0);

}

template <polyMPO::MaterialPointSlice mpSlice>
void setMPData(MPMesh_ptr p_mpmesh,
              const int nComps,
              const int numMPs,
              const double* mpDataIn){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == polyMPO::mpSliceToNumEntries<mpSlice>());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());
  
  kkViewHostU<const double**> mpDataIn_h(mpDataIn,nComps,numMPs);

  if (mpSlice == polyMPO::MPF_Cur_Pos_XYZ && p_MPs->rebuildOngoing()) {
    p_MPs->setRebuildMPSlice<polyMPO::MPF_Cur_Pos_XYZ>(mpDataIn_h);
    return;
  }

  auto mpData = p_MPs->getData<mpSlice>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpData_d("mpData_d",nComps,numMPs);
  Kokkos::deep_copy(mpData_d, mpDataIn_h);

  auto setData = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask){
      for (int i=0; i<nComps; i++)
        mpData(mp, i) = mpData_d(i, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setData, "setMPData");
  pumipic::RecordTime("PolyMPO_setMPData", timer.seconds());
}

template <polyMPO::MaterialPointSlice mpSlice>
void getMPData(MPMesh_ptr p_mpmesh,
                      const int nComps,
                      const int numMPs,
                      double* mpDataOut){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == polyMPO::mpSliceToNumEntries<mpSlice>());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpData = p_MPs->getData<mpSlice>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpDataCopy("mpDataCopy",nComps,numMPs);
  auto getData = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask){
      for (int i=0; i<nComps; i++)
        mpDataCopy(i,mpAppID(mp)) = mpData(mp,i);
    }
  };
  p_MPs->parallel_for(getData, "getMPData");
  kkViewHostU<double**> arrayHost(mpDataOut,nComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpDataCopy);
  pumipic::RecordTime("PolyMPO_getMPData", timer.seconds());
}

using setMPFunc = void (*)(MPMesh_ptr, const int, const int, const double*);
std::map<polyMPO::MaterialPointSlice, setMPFunc> setMPMap = {
  {polyMPO::MPF_Cur_Pos_Rot_Lat_Lon, setMPData<polyMPO::MPF_Cur_Pos_Rot_Lat_Lon>},
  {polyMPO::MPF_Tgt_Pos_Rot_Lat_Lon, setMPData<polyMPO::MPF_Tgt_Pos_Rot_Lat_Lon>},
  {polyMPO::MPF_Cur_Pos_XYZ, setMPData<polyMPO::MPF_Cur_Pos_XYZ>},
  {polyMPO::MPF_Tgt_Pos_XYZ, setMPData<polyMPO::MPF_Tgt_Pos_XYZ>},
  {polyMPO::MPF_Mass, setMPData<polyMPO::MPF_Mass>},
  {polyMPO::MPF_Vel, setMPData<polyMPO::MPF_Vel>},
  {polyMPO::MPF_Rot_Lat_Lon_Incr, setMPData<polyMPO::MPF_Rot_Lat_Lon_Incr>},
  {polyMPO::MPF_Strain_Rate, setMPData<polyMPO::MPF_Strain_Rate>},
  {polyMPO::MPF_Stress, setMPData<polyMPO::MPF_Stress>}
};

void polympo_setMPData_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpDataIn, const int mpDataType){
  polyMPO::MaterialPointSlice type = static_cast<polyMPO::MaterialPointSlice>(mpDataType);
  (*setMPMap[type])(p_mpmesh, nComps, numMPs, mpDataIn);
}

void polympo_setMPPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpPositionsIn){
  setMPData<polyMPO::MPF_Cur_Pos_XYZ>(p_mpmesh, nComps, numMPs, mpPositionsIn);
}
void polympo_getMPPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpPositionsHost){
  getMPData<polyMPO::MPF_Cur_Pos_XYZ>(p_mpmesh, nComps, numMPs, mpPositionsHost);
}
void polympo_setMPRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpRotLatLonIn){
  static int callCount = 0;
  PMT_ALWAYS_ASSERT(callCount == 0);
  setMPData<polyMPO::MPF_Cur_Pos_Rot_Lat_Lon>(p_mpmesh, nComps, numMPs, mpRotLatLonIn);
  callCount++;
}
void polympo_getMPRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpRotLatLonHost){
  getMPData<polyMPO::MPF_Cur_Pos_Rot_Lat_Lon>(p_mpmesh, nComps, numMPs, mpRotLatLonHost);
}
void polympo_setMPMass_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpMassIn) {
  setMPData<polyMPO::MPF_Mass>(p_mpmesh, nComps, numMPs, mpMassIn);
}
void polympo_getMPMass_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpMassHost) {
  getMPData<polyMPO::MPF_Mass>(p_mpmesh, nComps, numMPs, mpMassHost);
}
void polympo_setMPVel_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpVelIn) {
  setMPData<polyMPO::MPF_Vel>(p_mpmesh, nComps, numMPs, mpVelIn);
}
void polympo_getMPVel_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpVelHost) {
  getMPData<polyMPO::MPF_Vel>(p_mpmesh, nComps, numMPs, mpVelHost);
}
void polympo_setMPStrainRate_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpStrainRateIn){
  setMPData<polyMPO::MPF_Strain_Rate>(p_mpmesh, nComps, numMPs, mpStrainRateIn);
}
void polympo_getMPStrainRate_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpStrainRateHost){
  getMPData<polyMPO::MPF_Strain_Rate>(p_mpmesh, nComps, numMPs, mpStrainRateHost);
}
void polympo_setMPStress_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpStressIn){
  setMPData<polyMPO::MPF_Stress>(p_mpmesh, nComps, numMPs, mpStressIn);
}
void polympo_getMPStress_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpStressHost){
  getMPData<polyMPO::MPF_Stress>(p_mpmesh, nComps, numMPs, mpStressHost);
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

void polympo_setMeshNumVtxs_f(MPMesh_ptr p_mpmesh, const int numVtxs){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  p_mesh->setNumVtxs(numVtxs);
  p_mesh->setMeshVtxBasedFieldSize(); 
}

void polympo_getMeshNumVtxs_f(MPMesh_ptr p_mpmesh, int & numVtxs) {
  checkMPMeshValid(p_mpmesh); //chech vailidity
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  numVtxs = p_mesh->getNumVertices();
}

void polympo_setMeshNumElms_f(MPMesh_ptr p_mpmesh, const int numElms){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  auto elm2Vtx = polyMPO::IntVtx2ElmView("MeshElementsToVertices",numElms); 
  auto elm2Elm = polyMPO::IntElm2ElmView("MeshElementsToElements",numElms); 

  p_mesh->setNumElms(numElms);
  p_mesh->setElm2VtxConn(elm2Vtx);
  p_mesh->setElm2ElmConn(elm2Elm);
  p_mesh->setMeshElmBasedFieldSize();
}

void polympo_getMeshNumElms_f(MPMesh_ptr p_mpmesh, int & numElms) {
  checkMPMeshValid(p_mpmesh); //chech vailidity
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  numElms = p_mesh->getNumElements();
}

void polympo_setMeshNumEdgesPerElm_f(MPMesh_ptr p_mpmesh, const int nCells, const int* array){
  //chech vailidity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  PMT_ALWAYS_ASSERT(p_mesh->meshEditable());

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==nCells);
  auto nEdgesPerElm = create_mirror_view_and_copy(array, nCells);
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto elm2ElmConn = p_mesh->getElm2ElmConn();
  Kokkos::parallel_for("set nEdgesPerElm", nCells, KOKKOS_LAMBDA(const int elm){
    elm2VtxConn(elm,0) = nEdgesPerElm(elm);
    elm2ElmConn(elm,0) = nEdgesPerElm(elm);
  });
}

void polympo_setMeshElm2VtxConn_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int nCells, const int* array){
  //chech vailidity
  checkMPMeshValid(p_mpmesh);
  kkViewHostU<const int**> arrayHost(array,maxEdges,nCells); 
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

void polympo_setMeshElm2ElmConn_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int nCells, const int* array){
  //chech vailidity
  checkMPMeshValid(p_mpmesh);
  kkViewHostU<const int**> arrayHost(array,maxEdges,nCells); //Fortran is column-major
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

int polympo_getMeshFVtxType_f() {
  return polyMPO::MeshFType_VtxBased;
}

int polympo_getMeshFElmType_f() {
  return polyMPO::MeshFType_ElmBased;
}

template<polyMPO::MeshFieldIndex fieldType, typename... Args>
void setMeshData(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, Args... arrayArgs){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  const double* arrayIn[] = {arrayArgs...};

  //check the size
  // PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices);
  // PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==nCells);

  //copy the host array to the device
  auto meshField = p_mesh->getMeshField<fieldType>();
  auto meshField_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), meshField);
  for(int i=0; i<nVertices; i++)
  for(int j=0; j<nComps; j++)
    meshField_h(i, j) = arrayIn[j][i];
  Kokkos::deep_copy(meshField, meshField_h);
  pumipic::RecordTime("PolyMPO_setMeshData", timer.seconds());
}

template<polyMPO::MeshFieldIndex fieldType, typename... Args>
void getMeshData(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, Args... arrayArgs){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  double* arrayOut[] = {arrayArgs...};

  //check the size
  // PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices);
  // PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==nCells);
  
  //copy the device to host 
  auto meshField = p_mesh->getMeshField<fieldType>();
  auto meshField_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), meshField);
  for(int i=0; i<nVertices; i++)
  for(int j=0; j<nComps; j++)
    arrayOut[j][i] = meshField_h(i,j);
  pumipic::RecordTime("PolyMPO_getMeshData", timer.seconds());
}

template<polyMPO::MeshFieldIndex fieldType>
void setMeshDataContiguous(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* arrayIn){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkViewHostU<const double**> hostView(arrayIn,nComps,nVertices);
  Kokkos::View<double**> deviceView("meshDeviceView",nComps,nVertices);
  Kokkos::deep_copy(deviceView, hostView);

  auto vtxField = p_mesh->getMeshField<fieldType>();

  // //check the size
  // PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  // PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the host array to the device
  Kokkos::parallel_for("set mesh field", nVertices, KOKKOS_LAMBDA(const int iVtx){
    for (int j=0; j<nComps; j++)
      vtxField(iVtx,j) = deviceView(j,iVtx);
  });
  pumipic::RecordTime("PolyMPO_setMeshDataContiguous", timer.seconds());
}

template<polyMPO::MeshFieldIndex fieldType>
void getMeshDataContiguous(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* arrayOut){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkDbl2dViewHostU hostViewOut(arrayOut,nComps,nVertices);
  Kokkos::View<double**> deviceView("meshDeviceView",nComps,nVertices);

  auto vtxField = p_mesh->getMeshField<fieldType>();

  //check the size
  // PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  // PMT_ALWAYS_ASSERT(p_mesh->getNumVertices() == nVertices); 
  // PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the device array to the host
  Kokkos::parallel_for("get mesh field", nVertices, KOKKOS_LAMBDA(const int iVtx){
    for (int j=0; j<nComps; j++)
      deviceView(j,iVtx) = vtxField(iVtx,j);
  });
  Kokkos::deep_copy(hostViewOut, deviceView);
  pumipic::RecordTime("PolyMPO_getMeshDataContiguous", timer.seconds());
}

void polympo_setMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* xArray, const double* yArray, const double* zArray){
  setMeshData<polyMPO::MeshF_VtxCoords>(p_mpmesh, 3, nVertices, xArray, yArray, zArray);
}
void polympo_getMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, double* xArray, double* yArray, double* zArray){
  getMeshData<polyMPO::MeshF_VtxCoords>(p_mpmesh, 3, nVertices, xArray, yArray, zArray);
}
void polympo_setMeshVtxRotLat_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* latitude){
  setMeshData<polyMPO::MeshF_VtxRotLat>(p_mpmesh, 1, nVertices, latitude);
}
void polympo_getMeshVtxRotLat_f(MPMesh_ptr p_mpmesh, const int nVertices, double* latitude){
  getMeshData<polyMPO::MeshF_VtxRotLat>(p_mpmesh, 1, nVertices, latitude);
}
void polympo_setMeshElmCenter_f(MPMesh_ptr p_mpmesh, const int nCells, const double* xArray, const double* yArray, const double* zArray){
  setMeshData<polyMPO::MeshF_ElmCenterXYZ>(p_mpmesh, 3, nCells, xArray, yArray, zArray);
}
void polympo_getMeshElmCenter_f(MPMesh_ptr p_mpmesh, const int nCells, double* xArray, double* yArray, double* zArray){
  getMeshData<polyMPO::MeshF_ElmCenterXYZ>(p_mpmesh, 3, nCells, xArray, yArray, zArray);
}
void polympo_setMeshVtxVel_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* uVelIn, const double* vVelIn){
  setMeshData<polyMPO::MeshF_Vel>(p_mpmesh, 2, nVertices, uVelIn, vVelIn);
}
void polympo_getMeshVtxVel_f(MPMesh_ptr p_mpmesh, const int nVertices, double* uVelOut, double* vVelOut){
  getMeshData<polyMPO::MeshF_Vel>(p_mpmesh, 2, nVertices, uVelOut, vVelOut);
}
void polympo_setMeshVtxMass_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* vtxMass){
  setMeshData<polyMPO::MeshF_VtxMass>(p_mpmesh, 1, nVertices, vtxMass);
}
void polympo_getMeshVtxMass_f(MPMesh_ptr p_mpmesh, const int nVertices, double* vtxMass){
  getMeshData<polyMPO::MeshF_VtxMass>(p_mpmesh, 1, nVertices, vtxMass);
}
void polympo_setMeshElmMass_f(MPMesh_ptr p_mpmesh, const int nCells, const double* elmMass){
  setMeshData<polyMPO::MeshF_ElmMass>(p_mpmesh, 1, nCells, elmMass);
}
void polympo_getMeshElmMass_f(MPMesh_ptr p_mpmesh, const int nCells, double* elmMass){
  getMeshData<polyMPO::MeshF_ElmMass>(p_mpmesh, 1, nCells, elmMass);
}
void polympo_setMeshVtxOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* array) {
  setMeshDataContiguous<polyMPO::MeshF_OnSurfVeloIncr>(p_mpmesh, nComps, nVertices, array);
}
void polympo_getMeshVtxOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array) {
  getMeshDataContiguous<polyMPO::MeshF_OnSurfVeloIncr>(p_mpmesh, nComps, nVertices, array);
}
void polympo_setMeshVtxOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* array) {
  setMeshDataContiguous<polyMPO::MeshF_OnSurfDispIncr>(p_mpmesh, nComps, nVertices, array);
}
void polympo_getMeshVtxOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array) {
  getMeshDataContiguous<polyMPO::MeshF_OnSurfDispIncr>(p_mpmesh, nComps, nVertices, array);
}

void polympo_push_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh) ->push();
}

//TODO skeleton of reconstruction functions
void polympo_setReconstructionOfMass_f(MPMesh_ptr p_mpmesh, const int order, const int meshEntType){
  checkMPMeshValid(p_mpmesh);
  auto mpmesh = ((polyMPO::MPMesh*)p_mpmesh);
  polyMPO::MeshFieldType type = static_cast<polyMPO::MeshFieldType>(meshEntType);
  if (type == polyMPO::MeshFType_VtxBased)
    mpmesh->setReconstructSlice<polyMPO::MeshF_VtxMass>(order, type);
  if (type == polyMPO::MeshFType_ElmBased)
    mpmesh->setReconstructSlice<polyMPO::MeshF_ElmMass>(order, type);
}

void polympo_setReconstructionOfVel_f(MPMesh_ptr p_mpmesh, const int order, const int meshEntType){
  checkMPMeshValid(p_mpmesh);
  auto mpmesh = ((polyMPO::MPMesh*)p_mpmesh);
  polyMPO::MeshFieldType type = static_cast<polyMPO::MeshFieldType>(meshEntType);
  if (type == polyMPO::MeshFType_VtxBased)
    mpmesh->setReconstructSlice<polyMPO::MeshF_Vel>(order, type);
  else {
    std::cerr << "Error: This reconstruction is not supported\n";
    exit(1);
  }
}

void polympo_setReconstructionOfStrainRate_f(MPMesh_ptr p_mpmesh, const int order, const int meshEntType){
  checkMPMeshValid(p_mpmesh);
  std::cerr << "Error: This routine is not implemented yet\n";
  exit(1);
  (void)order;
  (void)meshEntType;
}

void polympo_setReconstructionOfStress_f(MPMesh_ptr p_mpmesh, const int order, const int meshEntType){
  checkMPMeshValid(p_mpmesh);
  std::cerr << "Error: This routine is not implemented yet\n";
  exit(1);
  (void)order;
  (void)meshEntType;
}

void polympo_applyReconstruction_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  auto mpmesh = ((polyMPO::MPMesh*)p_mpmesh);
  mpmesh->reconstructSlices();
}

void polympo_setOwningProc_f(MPMesh_ptr p_mpmesh, const int nCells, const int* array){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh; 
  PMT_ALWAYS_ASSERT(p_mesh->meshEditable());
  kkViewHostU<const int*> arrayHost(array,nCells); 

  //check the size
  PMT_ALWAYS_ASSERT(nCells == p_mesh->getNumElements());

  Kokkos::View<int*> owningProc("owningProc",nCells);
  Kokkos::deep_copy(owningProc, arrayHost);
  p_mesh->setOwningProc(owningProc);
}

void polympo_enableTiming_f(){
  pumipic::EnableTiming();
}

void polympo_summarizeTime_f(){
  pumipic::SummarizeTime();
}

void polympo_setTimingVerbosity_f(int v){
  pumipic::SetTimingVerbosity(v);
}

