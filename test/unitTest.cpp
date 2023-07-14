#include "pmt_MPMesh.hpp"
#include "pmt_assembly.hpp"
#include "pmo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

#define TEST_EPSILON 1e-6
using namespace polyMpmTest;

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
        //auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
        //                                      4 5 6 4 5 6 4 5 6 4 = 49 
        auto mpMesh = initTestMPMesh(testMesh, 1); //creates test MPs
        auto MPs = mpMesh.MPs;
        MPs->fillData<MPF_Mass>(1.0); //set MPF_Mass to 1.0
        MPs->fillData<MPF_Basis_Vals>(1.0);//TODO: change this to real basis value based on the basis computation routine
        //TODO: write PS_LAMBDA to assign 2 velocity component to be position[0] and [1]
        auto mpVel = MPs->getData<MPF_Vel>();
        auto mpCurPosXYZ = MPs->getData<MPF_Cur_Pos_XYZ>();
        auto setVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
            if(mask) { 
                for(int i=0; i<2; i++){
                    mpVel(mp,i) = mpCurPosXYZ(mp,i);
                }
            }
        };
        mpMesh.MPs->parallel_for(setVel, "setVel=CurPosXY");
        auto mesh = *mpMesh.mesh;
        PMT_ALWAYS_ASSERT(mesh.getNumVertices() == 19);
        PMT_ALWAYS_ASSERT(mesh.getNumElements() == 10);

        //run non-physical assembly (mp -to- mesh vertex) kernel
        polyMpmTest::assembly<MPF_Vel,MeshF_Vel>(mpMesh,false,false);//TODO: two flags not supported yet
        auto vtxField = mesh.getMeshField<MeshF_Vel>();
        //interpolateWachspress(mpMesh);
        //auto vtxFieldBasis = polyMpmTest::assemblyNew<MP_Cur_Pos_XYZ>(mpMesh,true);
        //check the result
        auto vtxField_h = Kokkos::create_mirror_view(vtxField);
        //auto vtxFieldFromMesh_h_ = Kokkos::create_mirror_view(vtxFieldFromMesh);
        Kokkos::deep_copy(vtxField_h, vtxField);
        const std::vector<std::vector<double>> vtxFieldExpected = {
            {1.415000, 1.450000}, {4.865000, 1.866667},
            {16.323333, 5.251333}, {9.305000, 3.774667},
            {9.380000, 6.418000}, {2.675000, 6.130000},
            {4.475000, 9.463333}, {4.995000, 7.816667},
            {6.720000, 7.543333}, {11.096429, 7.573714},
            {11.171429, 5.740381}, {11.564762, 5.532381},
            {7.964762, 4.305714}, {8.436429, 8.699048},
            {4.840000, 11.046667}, {1.260000, 4.680000},
            {3.040000, 7.713333}, {4.911429, 5.639048},
            {3.131429, 2.605714}}; 
        for(size_t i=0; i<vtxField_h.size(); i++) {
            int j = i/2;
            int k = i%2;
            auto res = polyMpmTest::isEqual(vtxField_h(j,k),vtxFieldExpected[j][k], TEST_EPSILON);
          if(!res) {
            fprintf(stderr, "expected != calc Value!\n\t[%d][%d]: %.6lf != %.6lf\n",
                                                j,k,vtxFieldExpected[j][k],vtxField_h(j,k));
          }
          PMT_ALWAYS_ASSERT(res);
        }

        interpolateWachspress(mpMesh);              
    }
    Kokkos::finalize();

    return 0;
}
