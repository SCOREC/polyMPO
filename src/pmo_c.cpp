#include "pmo_createTestMPMesh.hpp"
#include "pmo_c.h"
#include <stdio.h>

namespace{
    std::vector<MPMesh_ptr> p_mpmeshes;////store the p_mpmeshes that is legal
}

void polympo_initialize() {
  int isMPIInit;
  MPI_Initialized(&isMPIInit);
  PMT_ALWAYS_ASSERT(isMPIInit);
  Kokkos::initialize();
}

void polympo_finalize() {
  Kokkos::finalize();
}

MPMesh_ptr polympo_createMPMesh(int testMeshOption, int testMPOption) {
  polyMpmTest::Mesh* p_mesh;
  if(testMeshOption){
    int scaleFactor = 1;
    p_mesh = polyMpmTest::initTestMesh(testMeshOption, scaleFactor);
  }else{
    p_mesh = new polyMpmTest::Mesh();
  }
  polyMpmTest::MaterialPoints* p_mps;
  if(testMPOption){
    p_mps = polyMpmTest::initTestMPs(p_mesh, testMPOption);
  }else{
    p_mps = new polyMpmTest::MaterialPoints();  
  }
  MPMesh_ptr p_mpMeshReturn = (MPMesh_ptr) new polyMpmTest::MPMesh(p_mesh, p_mps);
  p_mpmeshes.push_back(p_mpMeshReturn);
  return p_mpMeshReturn;
}

void polympo_deleteMPMesh(MPMesh_ptr p_mpmesh) {
  //check mpMesh is valid
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());

  p_mpmeshes.erase(p_mpmeshIter);
  delete (polyMpmTest::MPMesh*)p_mpmesh;
}

void polympo_setMPICommunicator(MPI_Fint fcomm){
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

void polympo_setMPVelArray(MPMesh_ptr p_mpmesh, int size, double* array) {
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  auto p_MPs = ((polyMpmTest::MPMesh*)p_mpmesh)->p_MPs;
  kkVec2dViewHostU arrayHost(array,size);
  auto mpVel = p_MPs->getData<polyMpmTest::MPF_Vel>();
 
  auto mpVelCopy = polyMpmTest::DoubleVec2DView("mpVelNewValue",size);
  //copy the host array to the device
  Kokkos::deep_copy(mpVelCopy,arrayHost);//TODO: double check the cpu memory peak usage 
  
  //modify the MP array with the mpVelCopy copied from the host array
  auto setVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask) { 
      for(int i=0; i<vec2d_nEntries; i++){
        mpVel(mp,i) = mpVelCopy(mp,i);
      }
    }
  };
  p_MPs->parallel_for(setVel, "setVel to array");
}

void polympo_getMPVelArray(MPMesh_ptr p_mpmesh, int size, double* array) {
  //check mpMesh is valid
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  kkVec2dViewHostU arrayHost(array,size);

  auto p_MPs = ((polyMpmTest::MPMesh*)p_mpmesh)->p_MPs;
  auto mpVel = p_MPs->getData<polyMpmTest::MPF_Vel>();
  auto mpVelCopy = polyMpmTest::DoubleVec2DView("copyOfMPVel",size);
   
  PMT_ALWAYS_ASSERT(p_MPs->getCount()==size); 
  PMT_ALWAYS_ASSERT(static_cast<size_t>(p_MPs->getCount()*vec2d_nEntries)==mpVelCopy.size());

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
  p_MPs->parallel_for(copyVel, "copy mpVel to mpVelCopy");

  //copy the device array to the host
  Kokkos::deep_copy(arrayHost,mpVelCopy); 
}

void polympo_checkMeshSetting(MPMesh_ptr p_mpmesh, int maxEdges, int vertexDegree){
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  PMT_ALWAYS_ASSERT(maxEdges <= maxVtxsPerElm);
  PMT_ALWAYS_ASSERT(vertexDegree <=  maxElmsPerVtx);
}

void polympo_setMeshNumVtxs(MPMesh_ptr p_mpmesh, int numVtxs){
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  ((polyMpmTest::MPMesh*)p_mpmesh)->p_mesh->setNumVtxs(numVtxs);
}

void polympo_setMeshNumElms(MPMesh_ptr p_mpmesh, int numElms){
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  auto p_mesh = ((polyMpmTest::MPMesh*)p_mpmesh)->p_mesh;

  auto elm2Vtx = polyMpmTest::IntVtx2ElmView("MeshElementsToVertices",numElms); 
  auto elm2Elm = polyMpmTest::IntElm2ElmView("MeshElementsToElements",numElms); 

  p_mesh->setNumElms(numElms);
  p_mesh->setElm2VtxConn(elm2Vtx);
  p_mesh->setElm2ElmConn(elm2Elm);
}

void polympo_setMeshVtxCoords(MPMesh_ptr p_mpmesh, int size, double* xArray, double* yArray, double* zArray){
  //chech validity
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  auto p_mesh = ((polyMpmTest::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==size); 

  kkDblViewHostU xArrayHost(xArray,size); 
  kkDblViewHostU yArrayHost(yArray,size); 
  kkDblViewHostU zArrayHost(zArray,size); 
  
  //copy the host array to the device
  auto coordsArray = polyMpmTest::Vector2View("MeshVtxCoords",size);
  polyMpmTest::Vector2View::HostMirror h_coordsArray = Kokkos::create_mirror_view(coordsArray);
  for(int i=0; i<size; i++){
    //we only have Vector2 now,so zArray is not used
    h_coordsArray(i)[0] = xArrayHost(i);
    h_coordsArray(i)[1] = yArrayHost(i);
  }
  Kokkos::deep_copy(coordsArray, h_coordsArray);
  p_mesh->setVtxCoords(coordsArray);
}

void polympo_setMeshNumEdgesPerElm(MPMesh_ptr p_mpmesh, int size, int* array){
  //chech vailidity
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  kkIntViewHostU arrayHost(array,size);
  auto p_mesh = ((polyMpmTest::MPMesh*)p_mpmesh)->p_mesh;

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumElements()==size);
  
  polyMpmTest::IntView nEdgesPerElm("MeshNumEdgesPerElm",size);
  Kokkos::deep_copy(nEdgesPerElm, arrayHost);
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto elm2ElmConn = p_mesh->getElm2ElmConn();
  Kokkos::parallel_for("set nEdgesPerElm", size, KOKKOS_LAMBDA(const int elm){
    elm2VtxConn(elm,0) = nEdgesPerElm(elm);
    elm2ElmConn(elm,0) = nEdgesPerElm(elm);
  });
}

void polympo_setMeshElm2VtxConn(MPMesh_ptr p_mpmesh, int size1, int size2, int* array){
  //chech vailidity
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  kkInt2dViewHostU arrayHost(array,size1,size2);
  auto p_mesh = ((polyMpmTest::MPMesh*)p_mpmesh)->p_mesh; 

  //check the size
  PMT_ALWAYS_ASSERT(size1 == p_mesh->getNumElements());
  PMT_ALWAYS_ASSERT(size2 <= maxVtxsPerElm);
  
  Kokkos::View<int**> elm2VtxArray("MeshElementsToVertices",size1,size2);
  Kokkos::deep_copy(elm2VtxArray, arrayHost);
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  Kokkos::parallel_for("set elm2VtxConn", size1, KOKKOS_LAMBDA(const int elm){
    for(int i=0; i<size2; i++){
        elm2VtxConn(elm,i+1) = elm2VtxArray(elm,i);
    }  
  });
}

void polympo_setMeshElm2ElmConn(MPMesh_ptr p_mpmesh, int size1, int size2, int* array){
  //chech vailidity
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  kkInt2dViewHostU arrayHost(array,size1,size2);
  auto p_mesh = ((polyMpmTest::MPMesh*)p_mpmesh)->p_mesh; 

  //check the size
  PMT_ALWAYS_ASSERT(size1 == p_mesh->getNumElements());
  PMT_ALWAYS_ASSERT(size2 <= maxVtxsPerElm);
  
  Kokkos::View<int**> elm2ElmArray("MeshElementsToVertices",size1,size2);
  Kokkos::deep_copy(elm2ElmArray, arrayHost);
  auto elm2ElmConn = p_mesh->getElm2ElmConn();
  Kokkos::parallel_for("set elm2ElmConn", size1, KOKKOS_LAMBDA(const int elm){
    for(int i=0; i<size2; i++){
        elm2ElmConn(elm,i+1) = elm2ElmArray(elm,i);
    }  
  });
}

void polympo_setMeshVelArray(MPMesh_ptr p_mpmesh, int size, double* array) {
  //check mpMesh is valid
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  auto p_mesh = ((polyMpmTest::MPMesh*)p_mpmesh)->p_mesh;
  kkVec2dViewHostU arrayHost(array,size);

  auto vtxField = p_mesh->getMeshField<polyMpmTest::MeshF_Vel>();

  //check the size
  PMT_ALWAYS_ASSERT(static_cast<size_t>(size*vec2d_nEntries)==vtxField.size());

  //copy the host array to the device
  Kokkos::deep_copy(vtxField,arrayHost);
}

void polympo_getMeshVelArray(MPMesh_ptr p_mpmesh, int size, double* array) {
  //check mpMesh is valid
  auto p_mpmeshIter = std::find(p_mpmeshes.begin(),p_mpmeshes.end(),p_mpmesh);
  PMT_ALWAYS_ASSERT(p_mpmeshIter != p_mpmeshes.end());
  auto p_mesh = ((polyMpmTest::MPMesh*)p_mpmesh)->p_mesh;
  kkVec2dViewHostU arrayHost(array,size);

  auto vtxField = p_mesh->getMeshField<polyMpmTest::MeshF_Vel>();

  //check the size
  PMT_ALWAYS_ASSERT(p_mesh->getNumVertices()==size); 
  PMT_ALWAYS_ASSERT(static_cast<size_t>(size*vec2d_nEntries)==vtxField.size());

  //copy the device array to the host
  Kokkos::deep_copy(arrayHost, vtxField);
}

