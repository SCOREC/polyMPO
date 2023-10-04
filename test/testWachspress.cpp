#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"
#include "pmpo_assembly.hpp"
#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    const int testMeshOption = 1;
    const int replicateFactor = 100;
    const int testMPOption = 1;
    //test init Test Mesh and run assembly and Wachspress
    {
        auto m = polyMPO::Mesh();
        auto mp = polyMPO::MaterialPoints();
        
        auto mesh = initTestMesh(testMeshOption, replicateFactor);
        auto mpMesh = initTestMPMesh(mesh,testMPOption);
        
        //test assembly in assembly.hpp
        polyMPO::assembly<MPF_Vel,MeshF_Vel>(mpMesh,false,false);
        interpolateWachspress2DTest(mpMesh);

        //test new sphericalInterpolation function
        if(argc == 2){
            void* p;
            setWithMPASMeshByFortran(&p, argv[1], (int)strlen(argv[1]));
            auto mpmesh = initTestMPMesh(((MPMesh*)p)->p_mesh,testMPOption);
            auto p_mesh = mpmesh.p_mesh;
            auto p_MPs = mpmesh.p_MPs;
            if(p_mesh->getGeomType() == geom_spherical_surf){
                auto meshVel = p_mesh->getMeshField<MeshF_Vel>();
                auto meshVtxCoords = p_mesh->getMeshField<MeshF_VtxCoords>(); 
                auto setVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
                    if(mask) { 
                        for(int i=0; i<2; i++){
                            meshVel(elm,i) = meshVtxCoords(elm,i);
                        }   
                    }   
                };
                p_MPs->parallel_for(setVel, "setVel=CurPosXY");
               
                sphericalInterpolation<MeshF_Vel, MPF_Vel>(mpmesh);
                
                auto mpVel = p_MPs->getData<MPF_Vel>();
                auto mpCurPosXYZ = p_MPs->getData<MPF_Cur_Pos_XYZ>();
                auto check = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
                    if(mask && elm<3) { 
                        for(int i=0; i<2; i++){
                            //TODO: test outcome is incorrect
                            //if(std::abs(mpVel(mp,i)-mpCurPosXYZ(mp,i))>TEST_EPSILON)
                            //    printf("Not equal! ");
                            printf("(%d,%d):%e, %e\n",elm,i,mpVel(mp,i),mpCurPosXYZ(mp,i));
                        }   
                    }   
                };
                p_MPs->parallel_for(check, "check MPF_Vel = MPF_");
            }else{ 
                printf("Warning: No further planar test!\n");
            }
        }else{
            printf("Warning: No further file test!\n");
        }
    }
    
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
