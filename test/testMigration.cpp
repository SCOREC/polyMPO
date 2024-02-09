#include "pmpo_MPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

int main(int argc, char* argv[] ) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc,argv);

    {
        Mesh* mesh = NULL;
        MPMesh* mpmesh = NULL;
        void* p;
        setWithMPASMeshByFortran(&p, argv[1], (int)strlen(argv[1]));
        mpmesh = (MPMesh*)p;
        mesh = mpmesh->p_mesh;
    } 

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}