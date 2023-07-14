#include "pmt_MPMesh.hpp"
#include "pmt_assembly.hpp"
#include "pmo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

using namespace polyMpmTest;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    {
      auto testMesh = initTestMesh(1); //creates simple test mesh, '1' is a replication factor
      //auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
      auto mpMesh = initTestMPMesh(testMesh, 1); //creates test MPs based on testMesh and mpPerElement
      
      auto mesh = *mpMesh.mesh;
      PMT_ALWAYS_ASSERT(mesh.getNumVertices() == 19);
      PMT_ALWAYS_ASSERT(mesh.getNumElements() == 10);

      /* run the weighted assembly for Vector2 and scalar vertex fields */
      auto vtxVec2Field = polyMpmTest::wtVec2Assembly<MPF_Cur_Pos_XYZ>(mpMesh);
      //auto vtxScalarField = polyMpmTest::wtScaAssembly<MP_CUR_ELM_ID>(mpMesh); 
      
      auto nVtxs = mesh.getNumVertices();

      /* How can I do the loop over vertices of the mesh 
      and compute
      for (int i = 0; i < nVtxs, i++) {
        for (int j = 0; j < 2; j++) {
          vtxVec2Field(j) = vtxVec2Field(j) / vtxScalarField(i);
        }
      }
      */

      //check the result
      auto vtxVec2Field_h = Kokkos::create_mirror_view(vtxVec2Field);
      Kokkos::deep_copy(vtxVec2Field_h, vtxVec2Field);

      //auto nVtxs = testMesh.getNumVertices();
      
      printf("Number of vertices %d\n", nVtxs);
      for (int i = 0; i < nVtxs; i++) {
        std::cout << i << ", " << vtxVec2Field_h(i)[0] << ", " << vtxVec2Field_h(i)[1] << "\n";
      }
    }
    Kokkos::finalize();

    return 0;
}
