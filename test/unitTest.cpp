#include "pmt_MPMesh.hpp"
#include "pmt_assembly.hpp"
#include "testUtils.hpp"
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

//  test Vector2
    auto v = polyMpmTest::Vector2();
    v[0] = 1;
    v[1] = 2;
    PMT_ALWAYS_ASSERT(v[0] == 1);
    PMT_ALWAYS_ASSERT(v[1] == 2);
//  test Vector operator
    auto v1 = polyMpmTest::Vector2(1,2);
    auto v2 = polyMpmTest::Vector2(3,4);
    auto v3 = v1 + v2;
    PMT_ALWAYS_ASSERT(v3[0] == 4);
    PMT_ALWAYS_ASSERT(v3[1] == 6);
    auto v4 = v1 - v2;
    PMT_ALWAYS_ASSERT(v4[0] == -2);
    PMT_ALWAYS_ASSERT(v4[1] == -2);
    auto v5 = v1 * 2;
    PMT_ALWAYS_ASSERT(v5[0] == 2);
    PMT_ALWAYS_ASSERT(v5[1] == 4);
    auto v6 = -v1;
    PMT_ALWAYS_ASSERT(v6[0] == -1);
    PMT_ALWAYS_ASSERT(v6[1] == -2);
    auto v7 = v1.dot(v2);
    PMT_ALWAYS_ASSERT(v7 == 11);
    auto v8 = v1.cross(v2);
    PMT_ALWAYS_ASSERT(v8 == -2);
    auto v9 = v1.magnitude();
    PMT_ALWAYS_ASSERT(v9 - sqrt(5) < 1e-6);

    //run assembly and test Wachspress
    {
        //auto testMesh = Mesh::readMPASMesh("/path/to/mpas/mesh.nc"); //read from MPAS via netcdf
        auto testMesh = initTestMesh(1); //creates simple test mesh, '1' is a replication factor
        auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
        auto mpMesh = initTestMPMesh(testMesh, mpPerElement); //creates test MPs
        
        auto mesh = mpMesh.getMesh();
        PMT_ALWAYS_ASSERT(mesh.getNumVertices() == 19);
        PMT_ALWAYS_ASSERT(mesh.getNumElements() == 10);

        //run non-physical assembly (mp -to- mesh vertex) kernel
        auto vtxFieldOld = polyMpmTest::assembly(mpMesh);
        polyMpmTest::assemblyNew<MP_Cur_Pos_XYZ>(mpMesh);
        //interpolateWachspress(mpMesh);
        //auto vtxFieldBasis = polyMpmTest::assemblyNew<MP_Cur_Pos_XYZ>(mpMesh,true);
        auto vtxFieldNew = mpMesh.getMesh().getMeshField<meshFieldCurPosXYZ>();
        //check the result
        auto vtxField_h_Old = Kokkos::create_mirror_view(vtxFieldOld);
        auto vtxField_h_New = Kokkos::create_mirror_view(vtxFieldNew);
        Kokkos::deep_copy(vtxField_h_Old, vtxFieldOld);
        Kokkos::deep_copy(vtxField_h_New, vtxFieldNew);
        const std::vector<double> vtxFieldExpected = {
          1.768750, 4.528750, 17.660000, 8.228750,
          8.406250, 2.818750, 4.978750, 5.708750,
          6.486250, 10.551786, 13.014286, 15.114286,
          9.714286, 8.631786, 4.990000, 1.050000,
          2.830000, 5.694286, 3.914286
        };
        PMT_ALWAYS_ASSERT(vtxField_h_Old.size()*3 == vtxField_h_New.size());
        for(size_t i=0; i<vtxField_h_Old.size(); i++) {
          //auto res = polyMpmTest::isEqual(vtxField_h(i),vtxFieldExpected[i], 1e-6);
            auto res = polyMpmTest::isEqual(vtxField_h_Old(i),vtxField_h_New(i,0)+vtxField_h_New(i,1)+vtxField_h_New(i,2), 1e-6);
          //printf("%.f\n",vtxField_h_Basis[i]);
          if(!res) {
            fprintf(stderr, "%.6f = (%.6f, %.6f, %.6f)\n",
                    vtxField_h_Old(i),vtxField_h_New(i,0),vtxField_h_New(i,1),vtxField_h_New(i,2));
            //fprintf(stderr, "computed value for vtx %ld, %.6f, does not match expected value %.6f\n",
            //        i, vtxField_h(i), vtxFieldExpected[i]);
            //printf("computed value for vtx %ld, %.6f, does not match expected value %.6f\n",
            //        i, vtxField_h_Old(i), vtxField_h_New[i]);
          }
          PMT_ALWAYS_ASSERT(res);
        }

        interpolateWachspress(mpMesh);              
    }
    Kokkos::finalize();

    return 0;
}
