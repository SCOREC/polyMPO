#include "pmpo_MPMesh.hpp"
#include "pmpo_mesh.hpp"
#include "pmpo_assembly.hpp"
#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

using namespace polyMPO;

void matrixMultiply(double matrix[3][3], Vec3d &v, Vec3d &result);

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    //test Vec2d
    auto v = polyMPO::Vec2d();
    v[0] = 1;
    v[1] = 2;
    PMT_ALWAYS_ASSERT(v[0] == 1);
    PMT_ALWAYS_ASSERT(v[1] == 2);
    //test Vector operator
    auto v1 = polyMPO::Vec2d(1,2);
    auto v2 = polyMPO::Vec2d(3,4);
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
    PMT_ALWAYS_ASSERT(v9 - sqrt(5) < TEST_EPSILON);

    //test Vec3d
    auto v10 = polyMPO::Vec3d();
    v10[0] = 1;
    v10[1] = 2;
    v10[2] = 3;
    PMT_ALWAYS_ASSERT(v10[0] == 1);
    PMT_ALWAYS_ASSERT(v10[1] == 2);
    PMT_ALWAYS_ASSERT(v10[2] == 3);
    //test Vec3d operators
    auto v11 = polyMPO::Vec3d(1,2,3);
    auto v12 = polyMPO::Vec3d(4,5,6);
    auto v13 = v11 + v12;
    PMT_ALWAYS_ASSERT(v13[0] == 5);
    PMT_ALWAYS_ASSERT(v13[1] == 7);
    PMT_ALWAYS_ASSERT(v13[2] == 9);
    auto v14 = v11 - v12;
    PMT_ALWAYS_ASSERT(v14[0] == -3);
    PMT_ALWAYS_ASSERT(v14[1] == -3);
    PMT_ALWAYS_ASSERT(v14[2] == -3);
    auto v15 = v11 * 2;
    PMT_ALWAYS_ASSERT(v15[0] == 2);
    PMT_ALWAYS_ASSERT(v15[1] == 4);
    PMT_ALWAYS_ASSERT(v15[2] == 6);
    auto v16 = -v11;
    PMT_ALWAYS_ASSERT(v16[0] == -1);
    PMT_ALWAYS_ASSERT(v16[1] == -2);
    PMT_ALWAYS_ASSERT(v16[2] == -3);
    auto v17 = v11.dot(v12);
    PMT_ALWAYS_ASSERT(v17 == 32);
    auto v18 = v11.cross(v12);
    PMT_ALWAYS_ASSERT(v18[0] == -3);
    PMT_ALWAYS_ASSERT(v18[1] == 6);
    PMT_ALWAYS_ASSERT(v18[2] == -3);
    auto v19 = v11.magnitude();
    PMT_ALWAYS_ASSERT(v19 - sqrt(14) < TEST_EPSILON); 

    //test calc sphericalTriangleArea sphericalTriangleArea2
    double radius = 6371229.0;
    auto a1 = polyMPO::Vec3d(radius,0,0);
    auto b1 = polyMPO::Vec3d(0,radius,0);
    auto c1 = polyMPO::Vec3d(0,0,radius);
    double areaRef = 4 * TEST_PI / 8;
    auto area1_1 = polyMPO::sphericalTriangleArea(a1,b1,c1,radius);
    auto area1_2 = polyMPO::sphericalTriangleArea2(a1,b1,c1,radius);
    printf("area1_1: %.16e\n",area1_1);
    printf("area1_2: %.16e\n",area1_2);
    PMT_ALWAYS_ASSERT(Kokkos::fabs(area1_1 - areaRef) < TEST_EPSILON); 
    PMT_ALWAYS_ASSERT(Kokkos::fabs(area1_2 - areaRef) < TEST_EPSILON); 
    //rotate Z by 30 degrees
    double rotateZ[3][3] = {{std::sqrt(3)/2, -1.0/2,         0.0}, 
                            {1.0/2,          std::sqrt(3)/2, 0.0},
                            {0.0,            0.0,            1.0}};
    //rotate Y by 45 degrees
    double rotateY[3][3] = {{std::sqrt(2)/2,  0.0, std::sqrt(2)/2},   
                            {0.0,             1.0, 0.0},
                            {-std::sqrt(2)/2, 0.0, std::sqrt(2)/2}};
    Vec3d a2, b2, c2;
    matrixMultiply(rotateZ, a1, a2);
    matrixMultiply(rotateZ, b1, b2);
    matrixMultiply(rotateZ, c1, c2);
    auto area2_1 = polyMPO::sphericalTriangleArea(a2,b2,c2,radius);
    auto area2_2 = polyMPO::sphericalTriangleArea2(a2,b2,c2,radius);
    PMT_ALWAYS_ASSERT(Kokkos::fabs(area2_1 - areaRef) < TEST_EPSILON); 
    PMT_ALWAYS_ASSERT(Kokkos::fabs(area2_2 - areaRef) < TEST_EPSILON); 
    Vec3d a3, b3, c3;
    matrixMultiply(rotateY, a2, a3);
    matrixMultiply(rotateY, b2, b3);
    matrixMultiply(rotateY, c2, c3);
    auto area3_1 = polyMPO::sphericalTriangleArea(a3,b3,c3,radius);
    auto area3_2 = polyMPO::sphericalTriangleArea2(a3,b3,c3,radius);
    PMT_ALWAYS_ASSERT(Kokkos::fabs(area3_1 - areaRef) < TEST_EPSILON); 
    PMT_ALWAYS_ASSERT(Kokkos::fabs(area3_2 - areaRef) < TEST_EPSILON); 
    
    //this test is only designed to work with the following option values:
    const int testMeshOption = 1;
    const int replicateFactor = 1;
    const int testMPOption = 1;
    //run assembly and test Wachspress
    {
        auto testMesh = initTestMesh(testMeshOption,replicateFactor); 
        PMT_ALWAYS_ASSERT(testMesh->getMeshType() == mesh_general_polygonal);
        PMT_ALWAYS_ASSERT(testMesh->getGeomType() == geom_planar_surf);
        auto mpMesh = initTestMPMesh(testMesh, testMPOption); //creates test MPs 
        auto p_MPs = mpMesh.p_MPs;
        p_MPs->fillData<MPF_Mass>(1.0); //set MPF_Mass to 1.0
        p_MPs->fillData<MPF_Basis_Vals>(1.0);//TODO: change this to real basis value based on the basis computation routine
        //TODO: write PS_LAMBDA to assign 2 velocity component to be position[0] and [1]
        auto mpVel = p_MPs->getData<MPF_Vel>();
        auto mpCurPosXYZ = p_MPs->getData<MPF_Cur_Pos_XYZ>();
        auto setVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
            if(mask) { 
                for(int i=0; i<2; i++){
                    mpVel(mp,i) = mpCurPosXYZ(mp,i);
                }
            }
        };
        mpMesh.p_MPs->parallel_for(setVel, "setVel=CurPosXY");
        auto p_mesh = mpMesh.p_mesh;
        PMT_ALWAYS_ASSERT(p_mesh->getNumVertices() == 19);
        PMT_ALWAYS_ASSERT(p_mesh->getNumElements() == 10);

        //run non-physical assembly (mp -to- mesh vertex) kernel
        polyMPO::assembly<MPF_Vel,MeshF_Vel>(mpMesh,false,false);//TODO: two flags not supported yet
        auto vtxField = p_mesh->getMeshField<MeshF_Vel>();
        //interpolateWachspress(mpMesh);
        //auto vtxFieldBasis = polyMPO::assemblyNew<MP_Cur_Pos_XYZ>(mpMesh,true);
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
            auto res = polyMPO::isEqual(vtxField_h(j,k),vtxFieldExpected[j][k], TEST_EPSILON);
          if(!res) {
            fprintf(stderr, "expected != calc Value!\n\t[%d][%d]: %.6lf != %.6lf\n",
                                                j,k,vtxFieldExpected[j][k],vtxField_h(j,k));
          }
          PMT_ALWAYS_ASSERT(res);
        }

        interpolateWachspress2DTest(mpMesh);              
    }
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

void matrixMultiply(double matrix[3][3], Vec3d &v, Vec3d &result) {
    for (int i = 0; i < 3; i++) {
        result[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            result[i] += matrix[i][j] * v[j];
        }
    }
}
