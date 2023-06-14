#include "pmt_MPMesh.hpp"
#include "pmt_assembly.hpp"
#include "testUtils.hpp"
#include <mpi.h>

#define TEST_EPSILON 1e-6
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
        auto mesh = mpMesh.getMesh();
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
        {1.768750 ,1.812500 }, {4.528750 ,2.145833},//TODO: remove 0.0
        {17.660000,5.803333 }, {8.228750 ,3.735833},
        {8.406250 ,5.952500 }, {2.818750 ,5.712500},
        {4.978750 ,9.712500 }, {5.708750 ,8.845833},
        {6.486250 ,7.395833 }, {10.551786,7.397143},
        {13.014286,6.687143 }, {15.114286,7.137143},
        {9.714286 ,5.297143 }, {8.631786 ,8.840476},
        {4.990000 ,10.933333}, {1.050000 ,3.900000},
        {2.830000 ,6.933333 }, {5.694286 ,6.290476},
        {3.914286 ,3.257143 }};
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
