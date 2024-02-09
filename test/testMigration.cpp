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

        int numElms = mesh->getNumElements();
        auto elm2VtxConn = mesh->getElm2VtxConn();
        auto vtxCoords = mesh->getMeshField<polyMPO::MeshF_VtxCoords>(); 
        
        Vec3dView elmCenter("elementCenter",numElms);
        Kokkos::parallel_for("calcElementCenter", numElms, KOKKOS_LAMBDA(const int elm){  
            int numVtx = elm2VtxConn(elm,0);
            double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
            for(int i=1; i<= numVtx; i++){
                sum_x += vtxCoords(elm2VtxConn(elm,i)-1,0);
                sum_y += vtxCoords(elm2VtxConn(elm,i)-1,1);
                sum_z += vtxCoords(elm2VtxConn(elm,i)-1,2);
            }
            elmCenter(elm)[0] = sum_x/numVtx;
            elmCenter(elm)[1] = sum_y/numVtx;
            elmCenter(elm)[2] = sum_z/numVtx;
        });

        IntView owningProc("owningProc", numElms);
        Kokkos::parallel_for("setOwningProc", numElms, KOKKOS_LAMBDA(const int elm){
            owningProc(elm) = elmCenter(elm)[1] > 0 ? 1 : 0;
        });
        mesh->setMeshEdit(true);
        mesh->setOwningProc(owningProc);
    } 

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}