#include "pmpo_MPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

IntView procWedges(Mesh* mesh, int numElms, int comm_size) {
    auto elm2VtxConn = mesh->getElm2VtxConn();
    auto vtxLat = mesh->getMeshField<polyMPO::MeshF_VtxRotLat>(); 
    
    DoubleView elmLat("elmLat",numElms);
    DoubleView max("max", 1);
    DoubleView min("min", 1);
    Kokkos::parallel_for("calcElementLat", numElms, KOKKOS_LAMBDA(const int elm){  
        int numVtx = elm2VtxConn(elm,0);
        double sum = 0.0;
        for(int i=1; i<= numVtx; i++){
            sum += vtxLat(elm2VtxConn(elm,i)-1);
        }
        elmLat(elm) = sum/numVtx;
        Kokkos::atomic_max(&max(0), elmLat(elm));
        Kokkos::atomic_min(&min(0), elmLat(elm));
    });

    IntView owningProc("owningProc", numElms);
    Kokkos::parallel_for("setOwningProc", numElms, KOKKOS_LAMBDA(const int elm){
        double normalizedLat = (elmLat(elm) - min(0)) / (max(0) - min(0)) * .9;
        owningProc(elm) = normalizedLat * comm_size;
    });
    return owningProc;
}

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

        IntView owningProc = procWedges(mesh, numElms, comm_size);
        mesh->setMeshEdit(true);
        mesh->setOwningProc(owningProc);
        for (int i=0; i<5; i++)
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