#include "pmpo_MPMesh.hpp"
#include "pmpo_assembly.hpp"
#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

using namespace polyMPO;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    const int testMeshOption = 1;
    const int replicateFactor = 1;
    const int testMPOption = 1;
    {
      auto testMesh = initTestMesh(testMeshOption,replicateFactor); 
      auto mpMesh = initTestMPMesh(testMesh, testMPOption); //creates test MPs based on testMesh and mpPerElement
      
      auto p_mesh = mpMesh.p_mesh;
      PMT_ALWAYS_ASSERT(p_mesh->getNumVertices() == 19);
      PMT_ALWAYS_ASSERT(p_mesh->getNumElements() == 10);

      /* run the weighted assembly for Vec2d and scalar vertex fields */
      auto vtxVec2Field = polyMPO::wtVec2Assembly<MPF_Cur_Pos_XYZ>(mpMesh);
      
      auto nVtxs = p_mesh->getNumVertices();

      //check the result
      auto vtxVec2Field_h = Kokkos::create_mirror_view(vtxVec2Field);
      Kokkos::deep_copy(vtxVec2Field_h, vtxVec2Field);

      printf("Number of vertices %d\n", nVtxs);
      for (int i = 0; i < nVtxs; i++) {
        std::cout << i << ", " << vtxVec2Field_h(i)[0] << ", " << vtxVec2Field_h(i)[1] << "\n";
      }
    }
    Kokkos::finalize();
    MPI_Finalize();

    return 0;
}
