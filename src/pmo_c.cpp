#include "pmo_createTestMPMesh.hpp"
#include "pmo_c.h"
#include <stdio.h>

namespace{
    std::vector<mpmesh> mpMeshes;////store the mpMeshes that is legal
}

void polympo_initialize() {
  printf("polympo_initialize c++\n");
  int isMpiInit;
  MPI_Initialized(&isMpiInit);
  PMT_ALWAYS_ASSERT(isMpiInit);
  Kokkos::initialize();
}

void polympo_finalize() {
  printf("polympo_finalize c++\n");
  Kokkos::finalize();
}

mpmesh polympo_createMpMesh(int testMeshOption, int testMPOption) {
  polyMpmTest::Mesh* mesh;
  if(testMeshOption){
    int scaleFactor = 1;
    mesh = polyMpmTest::initTestMesh(scaleFactor, testMeshOption);
  }else{
    mesh = new polyMpmTest::Mesh();
  }
  polyMpmTest::MaterialPoints* mps;
  if(testMPOption){
    mps = polyMpmTest::initTestMPs(mesh, testMPOption);
  }else{
    mps = new polyMpmTest::MaterialPoints();  
  }
  mpmesh mpMeshReturn = (mpmesh) new polyMpmTest::MPMesh(mesh, mps);
  mpMeshes.push_back(mpMeshReturn);
  return mpMeshReturn;
}

void polympo_deleteMpMesh(mpmesh mpMeshIn) {
  //check mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());

  mpMeshes.erase(mpMeshInIter);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  delete mpMesh;
}

/**
 * Attention: this typedef is LayoutLeft, meaning that the first
 * index is the contiguous one. This matches the Fortran and GPU conventions for
 * allocations.
 */
typedef Kokkos::View<
          double*[vec2d_nEntries],
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>
        > kkDblViewHostU;//TODO:put it to mesh.hpp
                         
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

void polympo_setMP2dVelArray(mpmesh mpMeshIn, int rank1size, int rank2size, double* array) {
  //check mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_setMP2dVelArray c++ size %d %d\n", rank1size, rank2size);
  kkDbl2dViewHostU arrayHost(array,rank1size,rank2size);
  for(int i=0; i<rank1size; i++) {
    for(int j=0; j<rank2size; j++) {
      double kkval = arrayHost(i,j);
      double expected = i*rank2size+j;
      fprintf(stderr, "%d %d kkVal %.3f expected %.3f\n", i, j, kkval, expected);
    }
  }
}

void polympo_setMPVelArray(mpmesh mpMeshIn, int size, double* array) {
  //check mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_setMPVelArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);
  auto MPs = mpMesh->MPs;
  auto mpVel = MPs->getData<polyMpmTest::MPF_Vel>();
 
  auto mpVelCopy = polyMpmTest::DoubleVec2DView("mpVelNewValue",size);
  //copy the host array to the device
  Kokkos::deep_copy(mpVelCopy,arrayHost);
  
  //modify the MP array with the mpVelCopy copied from the host array
  auto setVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask) { 
      for(int i=0; i<vec2d_nEntries; i++){
        mpVel(mp,i) = mpVelCopy(mp,i);
      }
    }
  };
  MPs->parallel_for(setVel, "setVel to array");
}

void polympo_setMeshVelArray(mpmesh mpMeshIn, int size, double* array) {
  //check mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_setMeshVelArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);

  auto mesh = *(mpMesh->mesh);
  auto vtxField = mesh.getMeshField<polyMpmTest::MeshF_Vel>();

  //check the size
  PMT_ALWAYS_ASSERT(static_cast<size_t>(size*vec2d_nEntries)==vtxField.size());

  //copy the host array to the device
  Kokkos::deep_copy(vtxField,arrayHost);
}

void polympo_getMPVelArray(mpmesh mpMeshIn, int size, double* array) {
  //check mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_getMPVelArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);

  auto mesh = *(mpMesh->mesh);
  auto MPs = mpMesh->MPs;
  auto mpVel = MPs->getData<polyMpmTest::MPF_Vel>();
  auto mpVelCopy = polyMpmTest::DoubleVec2DView("copyOfMPVel",size);
    
  PMT_ALWAYS_ASSERT(static_cast<size_t>(MPs->getCount()*vec2d_nEntries)==mpVelCopy.size());

  //the pumipic 'slices' are not guaranteed to be contiguous, so we run a
  //pumipic parallel_for loop to fill a contigious array which will be copied to
  //a contiguous host array
  auto copyVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask) { 
      for(int i=0; i<vec2d_nEntries; i++){
        mpVelCopy(mp,i) = mpVel(mp,i);
      }
    }
  };
  MPs->parallel_for(copyVel, "copy mpVel to mpVelCopy");

  //copy the device array to the host
  Kokkos::deep_copy(arrayHost,mpVelCopy); 
}

void polympo_getMeshVelArray(mpmesh mpMeshIn, int size, double* array) {
  //check mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_getMeshVelArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);

  auto mesh = *mpMesh->mesh;
  auto vtxField = mesh.getMeshField<polyMpmTest::MeshF_Vel>();

  //check the size
  PMT_ALWAYS_ASSERT(static_cast<size_t>(size*vec2d_nEntries)==vtxField.size());

  //copy the device array to the host
  Kokkos::deep_copy(arrayHost, vtxField);
}

void polympo_setMPICommunicator(MPI_Fint fcomm){
    MPI_Comm comm = MPI_Comm_f2c(fcomm);
    int commSize;
    MPI_Comm_size(comm,&commSize);
    printf("polympo_setMPICommunicator with a communicator size: %d\n",commSize);
}

void polympo_setMeshVtxCoords(mpmesh mpMeshIn, int size, double* xArray, double* yArray, double* zArray){
  //chech validity
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  PMT_ALWAYS_ASSERT(size == mpMesh->mesh->getNumVertices());

  auto mesh = mpMesh->mesh;
  kkDblViewHostU xArrayHost(xArray,size); 
  kkDblViewHostU yArrayHost(yArray,size); 
  kkDblViewHostU zArrayHost(zArray,size); 

  //copy the host array to the device
  auto coordsArray = polyMpmTest::Vector2View("MeshVtxCoords",size);
  polyMpmTest::Vector2View::HostMirror h_coordsArray = Kokkos::create_mirror_view(coordsArray);
  for(int i=0; i<size; i++){
    //XXX: we only have Vector2 now,so zArray is not used
    h_coordsArray(i)[0] = xArrayHost(i,0);
    h_coordsArray(i)[1] = xArrayHost(i,1);
  }
  Kokkos::deep_copy(coordsArray, h_coordsArray);
  mesh->setVtxCoords(coordsArray);
}

void polympo_setMeshElm2VtxConn(mpmesh mpMeshIn, int size1, int size2, int* array){
  //chech vailidity
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  kkInt2dViewHostU arrayHost(array,size1,size2);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  auto mesh = *(mpMesh->mesh); 
  PMT_ALWAYS_ASSERT(size1 == mesh.getNumElements());
  PMT_ALWAYS_ASSERT(size2 <=  maxVtxsPerElm);
  
  Kokkos::View<int**> elm2VtxArray("MeshElementsToVertices",size1,size2);
  Kokkos::deep_copy(elm2VtxArray, arrayHost);
  auto elm2VtxConn = mesh.getElm2VtxConn();
  Kokkos::parallel_for("set elm2VtxConn", size1, KOKKOS_LAMBDA(const int elm){
    for(int i=0; i<size2; i++){
        elm2VtxConn(elm,i+1) = elm2VtxArray(elm,i);//TODO XXX will change to i,j =i,j
    }  
  });
}

void polympo_setMeshNumEdgesPerElm(mpmesh mpMeshIn, int size, int* array){
  //chech vailidity
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  kkIntViewHostU arrayHost(array,size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  auto mesh = *(mpMesh->mesh); 
  PMT_ALWAYS_ASSERT(size == mesh.getNumElements());
  
  polyMpmTest::IntView nEdgesPerElm("MeshNumEdgesPerElm",size);
  Kokkos::deep_copy(nEdgesPerElm, arrayHost);
  auto elm2VtxConn = mesh.getElm2VtxConn();
  Kokkos::parallel_for("set nEdgesPerElm", size, KOKKOS_LAMBDA(const int elm){
    elm2VtxConn(elm,0) = nEdgesPerElm(elm);//TODO XXX will change to i,j =i,j
  });
}

void polympo_setMeshVtx2ElmConn(mpmesh mpMeshIn, int size1, int size2, int* array){
  //chech validity
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  kkInt2dViewHostU arrayHost(array,size1,size2); 
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  auto mesh = *(mpMesh->mesh); 
  PMT_ALWAYS_ASSERT(size1 == mesh.getNumVertices());
  PMT_ALWAYS_ASSERT(size2 <=  maxElmsPerVtx);
  
  auto vtx2Elm = polyMpmTest::IntElm2VtxView("MeshVerticesToElement",size1); 
  polyMpmTest::IntElm2VtxView::HostMirror h_Vtx2Elm = Kokkos::create_mirror_view(vtx2Elm);
  for(int i=0; i<size1; i++){
    for(int j=0; j<size2; j++){
      h_Vtx2Elm(i,j+1) =  arrayHost(i,j);//TODO XXX will change to i,j =i,j
    } 
  }
  Kokkos::deep_copy(vtx2Elm, h_Vtx2Elm);
  mesh.setVtx2ElmConn(vtx2Elm);
}

void polympo_setMeshElm2ElmConn(mpmesh mpMeshIn, int size1, int size2, int* array){
  //chech validity
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  kkInt2dViewHostU arrayHost(array,size1,size2);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  auto mesh = *(mpMesh->mesh);
  PMT_ALWAYS_ASSERT(size1 == mesh.getNumElements());
  PMT_ALWAYS_ASSERT(size2 <=  maxVtxsPerElm);

  auto elm2Elm = polyMpmTest::IntElm2ElmView("MeshElementsToElements",size1); 
  polyMpmTest::IntElm2ElmView::HostMirror h_Elm2Elm = Kokkos::create_mirror_view(elm2Elm);
  for(int i=0; i<size1; i++){
    for(int j=0; j<size2; j++){
        h_Elm2Elm(i,j+1) = arrayHost(i,j);
    }  
  }
  Kokkos::deep_copy(elm2Elm, h_Elm2Elm);
  mesh.setElm2ElmConn(elm2Elm);
}

void polympo_checkMeshSetting(mpmesh mpMeshIn, int maxEdges, int vertexDegree){
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  PMT_ALWAYS_ASSERT(maxEdges <= maxVtxsPerElm);
  PMT_ALWAYS_ASSERT(vertexDegree <=  maxElmsPerVtx);
}

void polympo_setMeshNumVtxs(mpmesh mpMeshIn, int numVtxs){
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  mpMesh->mesh->setNumVtxs(numVtxs);
}

void polympo_setMeshNumElms(mpmesh mpMeshIn, int numElms){
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  auto elm2Vtx = polyMpmTest::IntVtx2ElmView("MeshElementsToVertices",numElms); 
  auto elm2Elm = polyMpmTest::IntElm2ElmView("MeshElementsToElements",numElms); 
  auto mesh = mpMesh->mesh;
  mesh->setNumElms(numElms);
  mesh->setElm2VtxConn(elm2Vtx);
  mesh->setElm2ElmConn(elm2Elm);
}

