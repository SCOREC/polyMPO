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

//initialize and finalize
void polympo_initialize_f() {
  int isMPIInit;
  MPI_Initialized(&isMPIInit);
  PMT_ALWAYS_ASSERT(isMPIInit);
  Kokkos::initialize();
}

void polympo_finalize_f() {
  Kokkos::finalize();
}

//create/delete MpMesh object
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

//set MPI communicator
void polympo_setMPICommunicator_f(MPMesh_ptr p_mpmesh, MPI_Fint fcomm){
  checkMPMeshValid(p_mpmesh);
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  p_MPs->setMPIComm(comm);
}

void polympo_startCommunication_f(MPMesh_ptr p_mpmesh){
  auto mpmesh = ((polyMPO::MPMesh*)p_mpmesh);
  mpmesh->startCommunication();  //Temporary last API in set mesh needs to be separate API
}

//MP info
void polympo_createMPs_f(MPMesh_ptr p_mpmesh,
                         const int numElms,
                         const int numMPs, // total number of MPs which is >= number of active MPs
                         int* mpsPerElm,
                         const int* mp2Elm,
                         const int* isMPActive){
  checkMPMeshValid(p_mpmesh);
  //the mesh must be fixed/set before adding MPs
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  PMT_ALWAYS_ASSERT(!p_mesh->meshEditable());
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements() == numElms);

  //Find the total no of MPs across all ranks
  //And loop over all MPs and find the smallest element id associated across a MP
  int numActiveMPs = 0;
  int minElmID = INT_MAX;
  for(int i = 0; i < numMPs; i++) {
    if(isMPActive[i] == MP_ACTIVE) {
      numActiveMPs++;
      if(mp2Elm[i] < minElmID) 
        minElmID = mp2Elm[i];
    }
  }
  int globalNumActiveMPs = 0;
  int globalMinElmID;
  MPI_Allreduce(&numActiveMPs, &globalNumActiveMPs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  MPI_Allreduce(&minElmID, &globalMinElmID, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); 
  PMT_ALWAYS_ASSERT(globalNumActiveMPs>0);

  //Loop over all mesh elements 0,1,... and find the first element that has an associated MP
  int firstElmWithMPs=INT_MAX;
  for (int i=0; i<numElms; i++) {
    if(mpsPerElm[i]) {
      firstElmWithMPs = i;
      break;
    }
  }
  int globalFirstElmWithMPs;
  MPI_Allreduce(&firstElmWithMPs, &globalFirstElmWithMPs, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  int offset = -1;
  if(globalMinElmID-globalFirstElmWithMPs==1) {
    offset = 1;
  }else if (globalMinElmID-globalFirstElmWithMPs==0){
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
  auto elm2global = p_mesh->getElmGlobal();
  auto mpsPerElm_d = create_mirror_view_and_copy(mpsPerElm, numElms);
  auto active_mp2Elm_d = create_mirror_view_and_copy(active_mp2Elm.data(), numActiveMPs);
  auto active_mpIDs_d = create_mirror_view_and_copy(active_mpIDs.data(), numActiveMPs);

  delete ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  ((polyMPO::MPMesh*)p_mpmesh)->p_MPs =
     new polyMPO::MaterialPoints(numElms, numActiveMPs, mpsPerElm_d, active_mp2Elm_d, active_mpIDs_d, elm2global);

  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  p_MPs->setElmIDoffset(offset);  
}

void polympo_startRebuildMPs_f(MPMesh_ptr p_mpmesh,
                               const int numMPs, // Total MPs which is GREATER than or equal to number of active MPs
                               const int* allMP2Elm,
                               const int* addedMPMask){
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

void polympo_startRebuildMPs2_f(MPMesh_ptr p_mpmesh,
                                const int sizeMP2elm,
                                const int* elem_ids,
                                const int nMPs_delete,
                                const int nMPs_add,
                                int* recvMPs_elm,
                                int* recvMPs_ids) {

  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  int offset = p_MPs->getElmIDoffset();

  for (int k=0; k<nMPs_add; k++){
    recvMPs_elm[k]=recvMPs_elm[k]-offset;
    recvMPs_ids[k]=recvMPs_ids[k]-offset;
  }
  
  auto elem_ids_d = create_mirror_view_and_copy(elem_ids, sizeMP2elm);
  auto recvMPs_elm_d = create_mirror_view_and_copy(recvMPs_elm, nMPs_add);
  auto recvMPs_ids_d = create_mirror_view_and_copy(recvMPs_ids, nMPs_add);

  Kokkos::View<int*> mp2Elm("mp2Elm", p_MPs->getCapacity());
  Kokkos::View<int*> numDeletedMPs_d("numDeletedMPs", 1);
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  auto setMP2Elm = PS_LAMBDA(const int& e, const int& mp, const int& mask) {
    if(mask) {
      mp2Elm(mp) = elem_ids_d(mpAppID(mp))-offset;
      if (mp2Elm(mp) == MP_DELETE){
          Kokkos::atomic_increment(&numDeletedMPs_d(0));
      }
    }
  };
  p_MPs->parallel_for(setMP2Elm, "setMP2Elm");
  int numDeletedMPs = pumipic::getLastValue(numDeletedMPs_d);
  assert(nMPs_delete==numDeletedMPs);
  p_MPs->startRebuild(mp2Elm, nMPs_add, recvMPs_elm_d, recvMPs_ids_d);
  pumipic::RecordTime("polympo_startRebuildMPs2_f", timer.seconds());
}

int polympo_getMPCount_f(MPMesh_ptr p_mpmesh) {
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  return p_MPs->getCount();
}

void polympo_finishRebuildMPs_f(MPMesh_ptr p_mpmesh){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  p_MPs->finishRebuild();
  pumipic::RecordTime("polympo_finishRebuildMPs_f", timer.seconds());
}

void polympo_setAppIDFunc_f(MPMesh_ptr p_mpmesh, IntVoidFunc getNext, void* appIDs) {
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  polyMPO::IntFunc getNextAppID = [getNext, appIDs]() { return getNext(appIDs); };
  p_MPs->setAppIDFunc(getNextAppID);
}

void polympo_getMPTgtElmID_f(MPMesh_ptr p_mpmesh, const int numMPs, int* elmIDs){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());
  auto mpTgtElmID = p_MPs->getData<polyMPO::MPF_Tgt_Elm_ID>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  auto elmIDoffset = p_MPs->getElmIDoffset();

  kkIntViewHostU arrayHost(elmIDs,numMPs);
  polyMPO::IntView mpTgtElmIDCopy("mpTgtElmIDNewValue",numMPs);

  auto setTgtElmId = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask){
        mpTgtElmIDCopy(mpAppID(mp)) = mpTgtElmID(mp)+elmIDoffset;
    }
  };
  p_MPs->parallel_for(setTgtElmId, "set mpTgtElmID");
  Kokkos::deep_copy( arrayHost, mpTgtElmIDCopy);
  pumipic::RecordTime("PolyMPO_getMPTgtElmID", timer.seconds());
}

void polympo_getMPCurElmID_f(MPMesh_ptr p_mpmesh, const int numMPs, int* elmIDs){
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());
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
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->p_mesh->setRotatedFlag(isRotateFlag>0);
}

//MP Slices
//Positions
void polympo_setMPPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpPositionsIn){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec3d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());
  kkViewHostU<const double**> mpPositionsIn_h(mpPositionsIn,nComps,numMPs);

  if (p_MPs->rebuildOngoing()) {
    p_MPs->setRebuildMPSlice<polyMPO::MPF_Cur_Pos_XYZ>(mpPositionsIn_h);
    return;
  }

  auto mpPositions = p_MPs->getData<polyMPO::MPF_Cur_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  
  Kokkos::View<double**> mpPositionsIn_d("mpPositionsDevice",vec3d_nEntries,numMPs);
  Kokkos::deep_copy(mpPositionsIn_d, mpPositionsIn_h);
  auto setPos = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask){
      mpPositions(mp,0) = mpPositionsIn_d(0, mpAppID(mp));
      mpPositions(mp,1) = mpPositionsIn_d(1, mpAppID(mp));
      mpPositions(mp,2) = mpPositionsIn_d(2, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setPos, "setMPPositions");
  pumipic::RecordTime("PolyMPO_setMPPositions", timer.seconds());
}

void polympo_getMPPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpPositionsHost){
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec3d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpPositions = p_MPs->getData<polyMPO::MPF_Cur_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpPositionsCopy("mpPositionsCopy",vec3d_nEntries,numMPs);
  auto getPos = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask){
      mpPositionsCopy(0,mpAppID(mp)) = mpPositions(mp,0);
      mpPositionsCopy(1,mpAppID(mp)) = mpPositions(mp,1);
      mpPositionsCopy(2,mpAppID(mp)) = mpPositions(mp,2);
    }
  };
  p_MPs->parallel_for(getPos, "getMPPositions");
  kkDbl2dViewHostU arrayHost(mpPositionsHost,nComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpPositionsCopy);
}

void polympo_setMPTgtPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpPositionsIn){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec3d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());
  kkViewHostU<const double**> mpPositionsIn_h(mpPositionsIn,nComps,numMPs);

  if (p_MPs->rebuildOngoing()) {
    p_MPs->setRebuildMPSlice<polyMPO::MPF_Cur_Pos_XYZ>(mpPositionsIn_h);
    return;
  }

  auto  mpPositions = p_MPs->getData<polyMPO::MPF_Tgt_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpPositionsIn_d("mpPositionsDevice",vec3d_nEntries,numMPs);
  Kokkos::deep_copy(mpPositionsIn_d, mpPositionsIn_h);
  auto setPos = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask){
      mpPositions(mp,0) = mpPositionsIn_d(0, mpAppID(mp));
      mpPositions(mp,1) = mpPositionsIn_d(1, mpAppID(mp));
      mpPositions(mp,2) = mpPositionsIn_d(2, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setPos, "setMPPositions");
  pumipic::RecordTime("PolyMPO_setMPTgtPositions", timer.seconds());
}

void polympo_getMPTgtPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpPositionsHost){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec3d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpPositions = p_MPs->getData<polyMPO::MPF_Tgt_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpPositionsCopy("mpPositionsCopy",vec3d_nEntries,numMPs);
  auto getPos = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask){
      mpPositionsCopy(0,mpAppID(mp)) = mpPositions(mp,0);
      mpPositionsCopy(1,mpAppID(mp)) = mpPositions(mp,1);
      mpPositionsCopy(2,mpAppID(mp)) = mpPositions(mp,2);
    }
  };
  p_MPs->parallel_for(getPos, "getMPPositions");
  kkDbl2dViewHostU arrayHost(mpPositionsHost,nComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpPositionsCopy);
  pumipic::RecordTime("PolyMPO_getMPTgtPositions", timer.seconds());
}

//LatLon
void polympo_setMPRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpRotLatLonIn){
  Kokkos::Timer timer;
  static int callCount = 0;
  PMT_ALWAYS_ASSERT(callCount == 0);
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpRotLatLon = p_MPs->getData<polyMPO::MPF_Cur_Pos_Rot_Lat_Lon>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  kkViewHostU<const double**> mpRotLatLonIn_h(mpRotLatLonIn,nComps,numMPs);
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
  pumipic::RecordTime("PolyMPO_setMPRotLatLon", timer.seconds());
}

void polympo_getMPRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpRotLatLonHost){
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
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
  kkDbl2dViewHostU arrayHost(mpRotLatLonHost,nComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpRotLatLonCopy);
}

void polympo_setMPTgtRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpRotLatLonIn){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpRotLatLon = p_MPs->getData<polyMPO::MPF_Tgt_Pos_Rot_Lat_Lon>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  kkViewHostU<const double**> mpRotLatLonIn_h(mpRotLatLonIn,nComps,numMPs);
  Kokkos::View<double**> mpRotLatLonIn_d("mpRotLatLonDevice",vec2d_nEntries,numMPs);
  Kokkos::deep_copy(mpRotLatLonIn_d, mpRotLatLonIn_h);
  auto setPos = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpRotLatLon(mp,0) = mpRotLatLonIn_d(0, mpAppID(mp));
      mpRotLatLon(mp,1) = mpRotLatLonIn_d(1, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setPos, "setMPTgtRotLatLon");
  pumipic::RecordTime("PolyMPO_setMPTgtRotLatLon", timer.seconds());
}

void polympo_getMPTgtRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpRotLatLonHost){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpRotLatLon = p_MPs->getData<polyMPO::MPF_Tgt_Pos_Rot_Lat_Lon>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpRotLatLonCopy("mpRotLatLonCopy",vec2d_nEntries,numMPs);
  auto getPos = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpRotLatLonCopy(0,mpAppID(mp)) = mpRotLatLon(mp,0);
      mpRotLatLonCopy(1,mpAppID(mp)) = mpRotLatLon(mp,1);
    }
  };
  p_MPs->parallel_for(getPos, "getMPRotLatLon");
  kkDbl2dViewHostU arrayHost(mpRotLatLonHost,nComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpRotLatLonCopy);
  pumipic::RecordTime("PolyMPO_getMPTgtRotLatLon", timer.seconds());
}

//MP Fields
void polympo_setMPMass_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpMassIn) {
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  int self;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  
  PMT_ALWAYS_ASSERT(nComps == 1); //TODO mp_sclr_t
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpMass = p_MPs->getData<polyMPO::MPF_Mass>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  kkViewHostU<const double**> mpMassIn_h(mpMassIn,nComps,numMPs);
  Kokkos::View<double**> mpMassIn_d("mpMassDevice",nComps,numMPs);
  Kokkos::deep_copy(mpMassIn_d, mpMassIn_h);
  auto setMPMass = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpMass(mp,0) = mpMassIn_d(0, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setMPMass, "setMPMass");
  pumipic::RecordTime("PolyMPO_setMPMass" + std::to_string(self), timer.seconds());
}

void polympo_getMPMass_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpMassHost) {
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == 1); //TODO mp_sclr_t
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpMass = p_MPs->getData<polyMPO::MPF_Mass>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpMassCopy("mpMassCopy",nComps,numMPs);
  auto getMPMass = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpMassCopy(0,mpAppID(mp)) = mpMass(mp,0);
    }
  };
  p_MPs->parallel_for(getMPMass, "getMPMass");
  kkDbl2dViewHostU arrayHost(mpMassHost,nComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpMassCopy);
  pumipic::RecordTime("PolyMPO_getMPMass", timer.seconds());
}

void polympo_setMPVel_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpVelIn) {
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpVel = p_MPs->getData<polyMPO::MPF_Vel>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  kkViewHostU<const double**> mpVelIn_h(mpVelIn,nComps,numMPs);
  Kokkos::View<double**> mpVelIn_d("mpVelDevice",vec2d_nEntries,numMPs);
  Kokkos::deep_copy(mpVelIn_d, mpVelIn_h);
  auto setMPVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpVel(mp,0) = mpVelIn_d(0, mpAppID(mp));
      mpVel(mp,1) = mpVelIn_d(1, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setMPVel, "setMPVel");
  pumipic::RecordTime("PolyMPO_setMPVel", timer.seconds());
}

void polympo_getMPVel_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpVelHost) {
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  auto mpVel = p_MPs->getData<polyMPO::MPF_Vel>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<double**> mpVelCopy("mpVelCopy",vec2d_nEntries,numMPs);
  auto getMPVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpVelCopy(0,mpAppID(mp)) = mpVel(mp,0);
      mpVelCopy(1,mpAppID(mp)) = mpVel(mp,1);
    }
  };
  p_MPs->parallel_for(getMPVel, "getMPVel");
  kkDbl2dViewHostU arrayHost(mpVelHost,nComps,numMPs);
  Kokkos::deep_copy(arrayHost, mpVelCopy);
  pumipic::RecordTime("PolyMPO_getMPVel", timer.seconds());
}

void polympo_setMPStrainRate_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  auto mpMesh = ((polyMPO::MPMesh*)p_mpmesh);
  mpMesh->calculateStrain(); 
}

void polympo_getMPStrainRate_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpStrainRateHost){
  checkMPMeshValid(p_mpmesh);
  std::cerr << "Error: This routine is not implemented yet\n";
  exit(1);
  (void)p_mpmesh;
  (void)nComps;
  (void)numMPs;
  (void)mpStrainRateHost;
}

void polympo_setMPStress_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  auto mpMesh = ((polyMPO::MPMesh*)p_mpmesh);
  mpMesh->calculateStress(); 
}

void polympo_getMPStress_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpStressHost){
  checkMPMeshValid(p_mpmesh);
  std::cerr << "Error: This routine is not implemented yet\n";
  exit(1);
  (void)p_mpmesh;
  (void)nComps;
  (void)numMPs;
  (void)mpStressHost;
}

void polympo_setAreaMP_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* areaMPHost){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  //Rank information
  int self;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  //Asserts
  PMT_ALWAYS_ASSERT(nComps == 1);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //MP Data
  auto mpArea = p_MPs->getData<polyMPO::MPF_Area>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  //Copy to device
  kkViewHostU<const double**> mpAreaIn_h(areaMPHost, nComps, numMPs);
  Kokkos::View<double**> mpAreaIn_d("mpAreaDevice", nComps, numMPs);
  Kokkos::deep_copy(mpAreaIn_d, mpAreaIn_h);
  //Set in PS
  auto setMPArea = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpArea(mp,0) = mpAreaIn_d(0, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setMPArea, "setMPArea");
  pumipic::RecordTime("PolyMPO_setMPArea" + std::to_string(self), timer.seconds());
}

void polympo_setIcePressureMP_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* icePressureMPHost){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  //Rank information
  int self;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  //Asserts
  PMT_ALWAYS_ASSERT(nComps == 1);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //MP Data
  auto mpIcePressure = p_MPs->getData<polyMPO::MPF_IcePressure>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  //Copy to device
  kkViewHostU<const double**> mpIcePressure_h(icePressureMPHost, nComps, numMPs);
  Kokkos::View<double**> mpIcePressure_d("mpIcePressureDevice", nComps, numMPs);
  Kokkos::deep_copy(mpIcePressure_d, mpIcePressure_h);
  //Set in PS
  auto setMPIcePressure = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpIcePressure(mp,0) = mpIcePressure_d(0, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setMPIcePressure, "setIcePressure");
  pumipic::RecordTime("PolyMPO_setIcePressure" + std::to_string(self), timer.seconds());
}

void polympo_setReplacementPressureMP_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* replacementPressureMPHost){
  Kokkos::Timer timer;
  checkMPMeshValid(p_mpmesh);
  
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  //Rank information
  int self;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  //Asserts
  PMT_ALWAYS_ASSERT(nComps == 1);
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  //MP Data
  auto mpReplacementPressure = p_MPs->getData<polyMPO::MPF_ReplacementPressure>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  //Copy to device
  kkViewHostU<const double**> mpReplacementPressure_h(replacementPressureMPHost, nComps, numMPs);
  Kokkos::View<double**> mpReplacementPressure_d("mpIcePressureDevice", nComps, numMPs);
  Kokkos::deep_copy(mpReplacementPressure_d, mpReplacementPressure_h);
  //Set in PS
  auto setMPReplacementPressure = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      mpReplacementPressure(mp,0) = mpReplacementPressure_d(0, mpAppID(mp));
    }
  };
  p_MPs->parallel_for(setMPReplacementPressure, "setReplacementPressure");
  pumipic::RecordTime("PolyMPO_setReplacementPressure" + std::to_string(self), timer.seconds());
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

//Owning Process and Global IDs
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

void polympo_setOwningProcVertex_f(MPMesh_ptr p_mpmesh, const int nVertices, const int* array){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh; 
  kkViewHostU<const int*> arrayHost(array,nVertices); 

  //check the size
  PMT_ALWAYS_ASSERT(nVertices == p_mesh->getNumVertices());

  Kokkos::View<int*> owningProcVertex("owningProcVerterx", nVertices);
  Kokkos::deep_copy(owningProcVertex, arrayHost);
  p_mesh->setOwningProcVertex(owningProcVertex);
}

void polympo_setElmGlobal_f(MPMesh_ptr p_mpmesh, const int nCells, const int* array){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh; 
  Kokkos::View<int*, Kokkos::HostSpace> arrayHost("arrayHost", nCells);
  for (int i = 0; i < nCells; i++) {
    arrayHost(i) = array[i] - 1;  // TODO right now elmID offset is set after MPs initialized
  }
  //check the size
  PMT_ALWAYS_ASSERT(nCells == p_mesh->getNumElements());

  Kokkos::View<int*> elmGlobal("elmGlobal",nCells);
  Kokkos::deep_copy(elmGlobal, arrayHost);
  p_mesh->setElmGlobal(elmGlobal);
}

void polympo_setVtxGlobal_f(MPMesh_ptr p_mpmesh, const int nVertices, const int* array){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh; 
  PMT_ALWAYS_ASSERT(nVertices==p_mesh->getNumVertices());
  Kokkos::View<int*, Kokkos::HostSpace> arrayHost("arrayHost", nVertices);
  for (int i = 0; i < nVertices; i++) {
    arrayHost(i) = array[i] - 1;  // TODO right now elmID offset is set after MPs initialized
  }

  Kokkos::View<int*> vtxGlobal("vtxGlobal",nVertices);
  Kokkos::deep_copy(vtxGlobal, arrayHost);
  p_mesh->setVtxGlobal(vtxGlobal);
}

//Mesh Fields
int polympo_getMeshFVtxType_f() {
  return polyMPO::MeshFType_VtxBased;
}

int polympo_getMeshFElmType_f() {
  return polyMPO::MeshFType_ElmBased;
}

void polympo_setMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* xArray, const double* yArray, const double* zArray){
  Kokkos::Timer timer;
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
  pumipic::RecordTime("PolyMPO_setMeshVtxCoords", timer.seconds());
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

void polympo_setMeshVtxRotLat_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* latitude){
  Kokkos::Timer timer;
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices); 

  //copy the host array to the device
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_VtxRotLat>();
  auto h_coordsArray = Kokkos::create_mirror_view(coordsArray);
  for(int i=0; i<nVertices; i++){
    h_coordsArray(i) = latitude[i];
  }
  Kokkos::deep_copy(coordsArray, h_coordsArray);
  pumipic::RecordTime("PolyMPO_setMeshVtxRotLat", timer.seconds());
}

void polympo_getMeshVtxRotLat_f(MPMesh_ptr p_mpmesh, const int nVertices, double* latitude){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices); 
  
  //copy the device to host 
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_VtxRotLat>();
  auto h_coordsArray = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                           coordsArray);
  for(int i=0; i<nVertices; i++){
    latitude[i] = h_coordsArray(i);
  }
}

void polympo_setMeshVtxVel_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* uVelIn, const double* vVelIn){
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices); 

  //copy the host array to the device
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_Vel>();
  auto h_coordsArray = Kokkos::create_mirror_view(coordsArray);
  for(int i=0; i<nVertices; i++){
    h_coordsArray(i,0) = uVelIn[i];
    h_coordsArray(i,1) = vVelIn[i];
  }
  Kokkos::deep_copy(coordsArray, h_coordsArray);
}

void polympo_getMeshVtxVel_f(MPMesh_ptr p_mpmesh, const int nVertices, double* uVelOut, double* vVelOut){
  Kokkos::Timer timer;
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices() == nVertices); 

  //copy the device array to the host
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_Vel>();
  auto h_coordsArray = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),coordsArray);
  for(int i=0; i<nVertices; i++){
    uVelOut[i] = h_coordsArray(i,0);
    vVelOut[i] = h_coordsArray(i,1);
  }
  pumipic::RecordTime("PolyMPO_getMeshVtxVel", timer.seconds());
}

void polympo_setMeshVtxMass_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* vtxMass){
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices); 

  //copy the host array to the device
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_VtxMass>();
  auto h_coordsArray = Kokkos::create_mirror_view(coordsArray);
  for(int i=0; i<nVertices; i++){
    h_coordsArray(i,0) = vtxMass[i];
  }
  Kokkos::deep_copy(coordsArray, h_coordsArray);
}

void polympo_getMeshVtxMass_f(MPMesh_ptr p_mpmesh, const int nVertices, double* vtxMass){
  Kokkos::Timer timer;
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  int self;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  
  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices() == nVertices); 

  //copy the device array to the host
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_VtxMass>();
  auto h_coordsArray = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),coordsArray);
  for(int i=0; i<nVertices; i++){
    vtxMass[i] = h_coordsArray(i,0);
  }
  pumipic::RecordTime("PolyMPO_getMeshVtxMass" + std::to_string(self), timer.seconds());
}

void polympo_setMeshElmMass_f(MPMesh_ptr p_mpmesh, const int nCells, const double* elmMass){
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==nCells); 

  //copy the host array to the device
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_ElmMass>();
  auto h_coordsArray = Kokkos::create_mirror_view(coordsArray);
  kkViewHostU<const double*> arrayHost(elmMass, nCells);
  for(int i=0; i<nCells; i++){
    h_coordsArray(i,0) = arrayHost[i];
  }
  Kokkos::deep_copy(coordsArray, h_coordsArray);
}

void polympo_getMeshElmMass_f(MPMesh_ptr p_mpmesh, const int nCells, double* elmMass){
  Kokkos::Timer timer;
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements() == nCells); 

  //copy the device array to the host
  auto coordsArray = p_mesh->getMeshField<polyMPO::MeshF_ElmMass>();
  auto h_coordsArray = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),coordsArray);
  for(int i=0; i<nCells; i++){
    elmMass[i] = h_coordsArray(i,0);
  }
  pumipic::RecordTime("PolyMPO_getMeshElmMass", timer.seconds());
}

//Increments in vertex velcoity and displacement
void polympo_setMeshVtxOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* array) {
  
  Kokkos::Timer timer;
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkViewHostU<const double**> arrayHost(array,nComps,nVertices);
  Kokkos::View<double**> array_d("meshVelIncrDevice",nComps,nVertices);
  Kokkos::deep_copy(array_d, arrayHost);

  auto vtxField = p_mesh->getMeshField<polyMPO::MeshF_OnSurfVeloIncr>();

  //check the size
  PMT_ALWAYS_ASSERT(nComps == vec2d_nEntries);
  PMT_ALWAYS_ASSERT(static_cast<size_t>(nVertices*vec2d_nEntries)==vtxField.size());

  //copy the host array to the device
  Kokkos::parallel_for("set mesh dispIncr", nVertices, KOKKOS_LAMBDA(const int iVtx){
    vtxField(iVtx,0) = array_d(0,iVtx);
    vtxField(iVtx,1) = array_d(1,iVtx);
  });
  pumipic::RecordTime("PolyMPO_setMeshVtxOnSurfVelIncr", timer.seconds());
}

void polympo_getMeshVtxOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array) {
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkDbl2dViewHostU arrayHost(array,nComps,nVertices);
  Kokkos::View<double**> array_d("meshVelIncrDevice",nComps,nVertices);

  auto vtxField = p_mesh->getMeshField<polyMPO::MeshF_OnSurfVeloIncr>();

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

void polympo_setMeshVtxOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* array) {
  Kokkos::Timer timer;
  //check mpMesh is valid
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  kkViewHostU<const double**> arrayHost(array,nComps,nVertices);
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
  pumipic::RecordTime("PolyMPO_setMeshVtxOnSurfDispIncr", timer.seconds());
}

void polympo_getMeshVtxOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array) {
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

void polympo_setMeshElmCenter_f(MPMesh_ptr p_mpmesh, const int nCells, const double* xArray, const double* yArray, const double* zArray){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==nCells); 

  //copy the host array to the device
  auto elmCenter = p_mesh->getMeshField<polyMPO::MeshF_ElmCenterXYZ>();
  auto h_elmCenter = Kokkos::create_mirror_view(elmCenter);
  for(int i=0; i<nCells; i++){
    h_elmCenter(i,0) = xArray[i];
    h_elmCenter(i,1) = yArray[i];
    h_elmCenter(i,2) = zArray[i];
  }
  Kokkos::deep_copy(elmCenter, h_elmCenter);
}

void polympo_getMeshElmCenter_f(MPMesh_ptr p_mpmesh, const int nCells, double* xArray, double* yArray, double* zArray){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==nCells); 

  //copy the device to host 
  auto elmCenter = p_mesh->getMeshField<polyMPO::MeshF_ElmCenterXYZ>();
  auto h_elmCenter = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), elmCenter);
  for(int i=0; i<nCells; i++){
    xArray[i] = h_elmCenter(i,0);
    yArray[i] = h_elmCenter(i,1);
    zArray[i] = h_elmCenter(i,2);
  }
}

void polympo_setMeshDualTriangleArea_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* areaTriangle){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices);
  //copy the host array to the device
  auto dualArea = p_mesh->getMeshField<polyMPO::MeshF_DualTriangleArea>();
  auto h_dualArea = Kokkos::create_mirror_view(dualArea);
  for(int i=0; i<nVertices; i++)
    h_dualArea(i,0) = areaTriangle[i];
  Kokkos::deep_copy(dualArea, h_dualArea);
}

void polympo_getMeshDualTriangleArea_f(MPMesh_ptr p_mpmesh, const int nVertices, double* areaTriangle){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices);
  //copy the device to host 
  auto dualArea = p_mesh->getMeshField<polyMPO::MeshF_DualTriangleArea>();
  auto h_dualArea = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dualArea);
  for(int i=0; i<nVertices; i++)
    areaTriangle[i] = h_dualArea(i,0);
}

void polympo_setGnomonicProjection_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  bool isRotated = p_mesh->getRotatedFlag();
  p_mesh->setGnomonicProjection(isRotated);
}

void polyMPO_setTanLatVertexRotatedOverRadius_f(MPMesh_ptr p_mpmesh, const int nVertices, double* array){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==nVertices);
  //copy the host array to the device
  auto tanLatVertexRotatedOverRadius = p_mesh->getMeshField<polyMPO::MeshF_TanLatVertexRotatedOverRadius>();
  auto h_tanLatVertexRotatedOverRadius = Kokkos::create_mirror_view(tanLatVertexRotatedOverRadius);
  for(int i=0; i<nVertices; i++)
    h_tanLatVertexRotatedOverRadius(i, 0) = array[i];
  Kokkos::deep_copy(tanLatVertexRotatedOverRadius, h_tanLatVertexRotatedOverRadius);
}

void polympo_setElasticTimeStep_f(MPMesh_ptr p_mpmesh, const double elasticTimeStep){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  p_mesh->setElasticTimeStep(elasticTimeStep);
}

void polympo_setDynamicTimeStep_f(MPMesh_ptr p_mpmesh, const double dynamicTimeStep){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;
  p_mesh->setDynamicTimeStep(dynamicTimeStep);
}

void polympo_setSolveStressMesh_f(MPMesh_ptr p_mpmesh, const int nCells, int* array){
  //chech validity
  checkMPMeshValid(p_mpmesh);
  auto p_mesh = ((polyMPO::MPMesh*)p_mpmesh)->p_mesh;

  PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==nCells);
  //copy the host array to the device
  auto solveStress = p_mesh->getMeshField<polyMPO::MeshF_SolveStress>();
  auto h_solveStress = Kokkos::create_mirror_view(solveStress);
  for(int i=0; i<nCells; i++)
    h_solveStress(i) = array[i];
  Kokkos::deep_copy(solveStress, h_solveStress);
}

//Advection Calcualtions
void polympo_push_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->push();
}

void polympo_push_ahead_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->push_ahead();
}

bool polympo_push1P_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  bool is_migrating=((polyMPO::MPMesh*)p_mpmesh)->push1P();
  return is_migrating;
}

void polympo_push_swap_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->push_swap();
}

void polympo_push_swap_pos_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  ((polyMPO::MPMesh*)p_mpmesh)->push_swap_pos();
}

// Reconstruction of variables from MPs to mesh vertices
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

//Simpler/cleaner way of calling reconstruction
void polympo_reconstruct_coeff_with_MPI_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  auto mpmesh = ((polyMPO::MPMesh*)p_mpmesh);
  mpmesh->reconstruct_coeff_full();
}

void polympo_reconstruct_iceArea_with_MPI_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  auto mpmesh = ((polyMPO::MPMesh*)p_mpmesh);
  mpmesh->assemblyVtx1<polyMPO::MeshF_VtxMass>();
}

void polympo_reconstruct_velocity_with_MPI_f(MPMesh_ptr p_mpmesh){
  checkMPMeshValid(p_mpmesh);
  auto mpmesh = ((polyMPO::MPMesh*)p_mpmesh);
  mpmesh->assemblyVtx1<polyMPO::MeshF_Vel>();
}

//Timing
void polympo_enableTiming_f(){
  pumipic::EnableTiming();
}

void polympo_summarizeTime_f(){
  pumipic::SummarizeTime();
}

void polympo_setTimingVerbosity_f(int v){
  pumipic::SetTimingVerbosity(v);
}

