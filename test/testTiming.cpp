#include "pmpo_MPMesh.hpp"
#include "pmpo_assembly.hpp"
#include "pmpo_wachspressBasis.hpp"
#include "pmpo_c.h"

#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

int main(int argc, char* argv[] ) {
    PMT_ALWAYS_ASSERT(argc <= 3);
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc,argv);
    int testMeshOption = 0;
    if(argc == 2) // no mesh read
        testMeshOption = 1;
    int replicateFactor = 1;
    if(argc != 1)
        replicateFactor = atoi(argv[1]);
    const int testMPOption = 1;
    printf("Time assembly and wachspress with factor: %d\n",replicateFactor);
    {
        Mesh* mesh = NULL;
        MPMesh* mpmesh = NULL;
        if(testMeshOption == 0){
            void* p;
            polympo_setWithMPASMeshByFortran(&p, argv[2], (int)strlen(argv[2]));
            mpmesh = (MPMesh*)p;
            mesh = mpmesh->p_mesh;
            //mesh = replicateMesh(mesh, replicateFactor); 
        }else{
            mesh = initTestMesh(testMeshOption, replicateFactor);
        }
        auto mpMesh = initTestMPMesh(mesh,testMPOption);
     
        printf("Total MPs:%d\n",mpMesh.p_MPs->getCount());
        
        for(int i=0; i<1; i++){
            //polyMPO::assemblyV0(mpMesh);
            if(mesh->getGeomType() == geom_spherical_surf)
                interpolateWachspress3DTest(mpMesh);
            else
                interpolateWachspress2DTest(mpMesh);
        }
        //printVTP(mpMesh);
    } 
    Kokkos::finalize();
    return 0;
}
