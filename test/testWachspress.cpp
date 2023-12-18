#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"
#include "pmpo_assembly.hpp"
#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

void test_planar(const int replicateFactor, const int testMPOption) {
    //test init Test Mesh and run assembly and Wachspress
    for (int testMeshOption = 1; testMeshOption <= 2; testMeshOption++){

        auto mesh = initTestMesh(testMeshOption,replicateFactor);
        auto mpMesh = initTestMPMesh(mesh,testMPOption);
        
        //test assembly in assembly.hpp
        polyMPO::assembly<MPF_Vel,MeshF_Vel>(mpMesh,false,false);
        interpolateWachspress2DTest(mpMesh);
        interpolateWachspress3DTest(mpMesh,testMeshOption);
    }
}

void test_spherical(const int testMPOption) {
    //test init Test Mesh and run assembly and Wachspress on spherical surface
    void* meshP;
    char* filename = (char *)malloc(sizeof(char) * 256); 
    
    // TODO: add relative path
    //sprintf(filename,"sample_mpas_meshes/spherical_cvt_642elms.nc"); 
    sprintf(filename,"/gpfs/u/home/MPMS/MPMSsngj/scratch/polyDev/polyMPO/test/sample_mpas_meshes/spherical_cvt_642elms.nc");
    setWithMPASMeshByFortran(&meshP, filename, (int)strlen(filename));
    auto mpmesh = (MPMesh*)meshP;
    auto mesh = mpmesh->p_mesh;
    auto mpMesh = initTestMPMesh(mesh, testMPOption);    
   
    interpolateWachspressSphericalTest(mpMesh);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    const int replicateFactor = 100;
    const int testMPOption = 1;
    test_planar(replicateFactor, testMPOption);
    test_spherical(testMPOption); 
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
