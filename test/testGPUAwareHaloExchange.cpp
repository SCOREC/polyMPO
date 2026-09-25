#include "pmpo_MPMesh.hpp"
#include "pmpo_createTestMPMesh.hpp"

#include <mpi.h>
#include <Kokkos_Core.hpp>

#include <iostream>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    int testResult = 0;

    {
        int rank = -1;
        int size = -1;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        if (rank == 0) {
            std::cout
                << "GPU-aware halo exchange test running with "
                << size << " MPI ranks."
                << std::endl;
        }

        // This test uses the following communication pattern:
        //
        // Rank 0 ----\
        // Rank 1 ----- > Rank 3
        // Rank 2 ----/
        //                |
        //                v
        //              Rank 0
        //
        // Therefore:
        // Rank 0 has one halo owned by Rank 3.
        // Rank 1 has no halos.
        // Rank 2 has no halos.
        // Rank 3 has three halos owned by Ranks 0, 1, and 2.

        if (size != 4) {
            if (rank == 0) {
                std::cerr
                    << "This test requires exactly 4 MPI ranks (got "
                    << size << "); skipping."
                    << std::endl;
            }
            Kokkos::finalize();
            MPI_Finalize();
            return 77;
        }

        // Create the existing polyMPO test mesh and MPMesh.

        polyMPO::Mesh* mesh = polyMPO::initTestMesh(1, 1);

        polyMPO::MPMesh mpMesh = polyMPO::initTestMPMesh(mesh, 1);

        // initTestMPMesh() creates MaterialPoints, but the test helper
        // does not initialize its MPI communicator.
        mpMesh.p_MPs->setMPIComm(MPI_COMM_WORLD);

        const int nVertices = mesh->getNumVertices();

        // startCommunication() expects local vertices to be ordered:

        polyMPO::IntView owningProcVertex("testOwningProcVertex", nVertices);

        polyMPO::IntView globalVtx("testGlobalVtx",nVertices);

        auto owningProcHost = Kokkos::create_mirror_view(owningProcVertex);

        auto globalVtxHost = Kokkos::create_mirror_view(globalVtx);

        // Number of locally-owned vertices differs by rank.
        // Rank 0: 18 owners + 1 halo
        // Rank 1: 19 owners + 0 halos
        // Rank 2: 19 owners + 0 halos
        // Rank 3: 16 owners + 3 halos

        int numLocalOwners = nVertices;

        if (rank == 0) {
            numLocalOwners = nVertices - 1;
        }
        else if (rank == 3) {
            numLocalOwners = nVertices - 3;
        }

        // Assign locally-owned vertices.
        // Rank 0 global IDs:   0,   1,   2, ...
        // Rank 1 global IDs: 100, 101, 102, ...
        // Rank 2 global IDs: 200, 201, 202, ...
        // Rank 3 global IDs: 300, 301, 302, ...

        for (int i = 0; i < numLocalOwners; ++i) {
            owningProcHost(i) = rank;
            globalVtxHost(i) = rank * 100 + i;
        }

        // Define halo vertices.

        if (rank == 0) {

            // Rank 0 receives Rank 3 owner vertex 0.
            // Rank 3 owner vertex 0 has global ID 300.

            owningProcHost(nVertices - 1) = 3;
            globalVtxHost(nVertices - 1) = 300;
        }
        else if (rank == 3) {

            // Rank 3 receives Rank 0 owner vertex 0.
            owningProcHost(nVertices - 3) = 0;
            globalVtxHost(nVertices - 3) = 0;

            // Rank 3 receives Rank 1 owner vertex 0.
            owningProcHost(nVertices - 2) = 1;
            globalVtxHost(nVertices - 2) = 100;

            // Rank 3 receives Rank 2 owner vertex 0.
            owningProcHost(nVertices - 1) = 2;
            globalVtxHost(nVertices - 1) = 200;
        }

        Kokkos::deep_copy(owningProcVertex, owningProcHost);

        Kokkos::deep_copy(globalVtx, globalVtxHost);

        mesh->setOwningProcVertex(owningProcVertex);
        mesh->setVtxGlobal(globalVtx);

        // Build owner/halo communication metadata
        mpMesh.startCommunication();

        std::cout
            << "Rank " << rank
            << ": owners = " << mpMesh.numOwnersTot
            << ", halos = " << mpMesh.numHalosTot
            << std::endl;

        // Create one scalar field value per vertex.

        const int nEntities = mpMesh.numOwnersTot + mpMesh.numHalosTot;

        const int numEntries = 1;

        Kokkos::View<double**> field("gpuAwareHaloTestField", nEntities, numEntries);

        const int numOwners = mpMesh.numOwnersTot;

        // Owner values are deterministic:
        // Rank 0 owner 0 =    0
        // Rank 1 owner 0 = 1000
        // Rank 2 owner 0 = 2000
        // Rank 3 owner 0 = 3000
        // All halos start at -1.

        Kokkos::parallel_for("initialize_gpu_aware_halo_test", nEntities,KOKKOS_LAMBDA(const int i)
            {
                if (i < numOwners) {
                    field(i, 0) =
                        1000.0 * static_cast<double>(rank)
                        + static_cast<double>(i);
                }
                else {
                    field(i, 0) = -1.0;
                }
            });

        Kokkos::fence();

        // Perform owner -> halo assignment.
        // mode = 1 : owner -> halo
        // op   = 1 : assignment

        const int mode = 1;
        const int op = 1;

        mpMesh.communicate_and_take_halo_contributions_gpu_aware(field, nEntities, numEntries, mode, op);

        Kokkos::fence();

        // Copy results to host and verify exact halo values.

        auto fieldHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), field);

        int localFailures = 0;

        if (rank == 0) {

            // Rank 0 receives Rank 3 owner vertex 0.
            const int haloIndex = nVertices - 1;
            const double expectedValue = 3000.0;

            if (fieldHost(haloIndex, 0) != expectedValue) {

                ++localFailures;

                std::cerr
                    << "Rank 0: halo value = "
                    << fieldHost(haloIndex, 0)
                    << ", expected = "
                    << expectedValue
                    << std::endl;
            }
        }
        else if (rank == 3) {

            // Rank 3 receives owner vertex 0 from Ranks 0, 1, and 2.

            const int haloFromRank0 = nVertices - 3;
            const int haloFromRank1 = nVertices - 2;
            const int haloFromRank2 = nVertices - 1;

            if (fieldHost(haloFromRank0, 0) != 0.0) {

                ++localFailures;

                std::cerr
                    << "Rank 3: halo from Rank 0 = "
                    << fieldHost(haloFromRank0, 0)
                    << ", expected = 0"
                    << std::endl;
            }

            if (fieldHost(haloFromRank1, 0) != 1000.0) {

                ++localFailures;

                std::cerr
                    << "Rank 3: halo from Rank 1 = "
                    << fieldHost(haloFromRank1, 0)
                    << ", expected = 1000"
                    << std::endl;
            }

            if (fieldHost(haloFromRank2, 0) != 2000.0) {

                ++localFailures;

                std::cerr
                    << "Rank 3: halo from Rank 2 = "
                    << fieldHost(haloFromRank2, 0)
                    << ", expected = 2000"
                    << std::endl;
            }
        }

        int globalFailures = 0;

        MPI_Allreduce(&localFailures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0) {

            if (globalFailures == 0) {
                std::cout
                    << "GPU-aware halo exchange test PASSED."
                    << std::endl;
            }
            else {
                std::cerr
                    << "GPU-aware halo exchange test FAILED with "
                    << globalFailures
                    << " halo errors."
                    << std::endl;
            }
        }

        if (globalFailures != 0) {
            testResult = 1;
        }

        // Do not delete mesh here.
        // mpMesh owns the Mesh and MaterialPoints objects.
    }

    Kokkos::finalize();
    MPI_Finalize();

    return testResult;
}
