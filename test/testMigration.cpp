#include "pmpo_MPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

int main(int argc, char* argv[] ) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc,argv);

    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    {
        MPMesh* mpmesh = NULL;
        void* p;
        setWithMPASMeshByFortran(&p, argv[1], (int)strlen(argv[1]));
        mpmesh = (MPMesh*)p;
        Mesh* mesh = mpmesh->p_mesh;
        int numElms = mesh->getNumElements();
        int numMPs = numElms;

        IntView mpsPerElm("mpsPerElm", numElms);
        IntView activeMP2Elm("activeMP2Elm", numMPs);
        IntView activeMPIDs("activeMPIDs", numMPs);

        Kokkos::parallel_for("setMPsPerElm", numElms, KOKKOS_LAMBDA(const int elm){
            mpsPerElm(elm) = 1;
        });
        Kokkos::parallel_for("setMPs", numMPs, KOKKOS_LAMBDA(const int mp){
            activeMP2Elm(mp) = mp;
            activeMPIDs(mp) = mp;
        });

        MaterialPoints p_MPs(numElms, numMPs, mpsPerElm, activeMP2Elm, activeMPIDs);
        mpmesh->p_MPs = &p_MPs;

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
            owningProc(elm) = elmCenter(elm)[1] > 0 ? comm_size-1 : 0;
        });
        mesh->setMeshEdit(true);
        mesh->setOwningProc(owningProc);
        mpmesh->push();

        auto MPs2Proc = p_MPs.getData<MPF_Tgt_Proc_ID>();
        auto checkPostBackMigrate = PS_LAMBDA(const int& e, const int& mp, const bool& mask) {
            if (mask) {
                assert(MPs2Proc(mp) == comm_rank);
                assert(owningProc(e) == comm_rank);
            }
        };
        p_MPs.parallel_for(checkPostBackMigrate, "checkPostBackMigrate");
    } 

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}