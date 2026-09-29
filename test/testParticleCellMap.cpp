#include "pmpo_MPMesh.hpp"
#include "pmpo_createTestMPMesh.hpp"

#include <mpi.h>
#include <Kokkos_Core.hpp>

#include <cmath>
#include <iostream>
#include <vector>


// Unit test for particle-cell mapping of the category fields on a SPHERICAL mesh
//
// Mesh: the 10-cell polyMPO test mesh (same cells and connectivity), with its
// vertices placed on a sphere of Earth radius, away from the poles.
//
// MPs: built by this test so the number of MPs per cell is controlled.
// Cells are split into four groups by (cell mod 4):
//
//   Group A (0): 2 MPs per vertex, spread inside the cell, all carry base
//   Group B (1): 2 MPs per vertex, alternating base and base + 0.5
//   Group C (2): no MPs
//   Group D (3): exactly 1 MP, carries base
//
//   base(f, k, cell) = f * 100000 + (k + 1) * 1000 + cell
//
// Part 0: MP count per cell matches the layout above.
//
// Part 1: MPs -> cells (MPMesh::mapMPsToCells)
//   Group A, D -> cell value exactly base
//   Group B    -> cell value in [base, base + 0.5]   (clamping)
//   Group C    -> cell value exactly 0
//   All three fields are mapped before any is read back.
//
// Part 2: cells -> MPs (MPMesh::mapCellsToMPs)

namespace {

constexpr double SPHERE_RADIUS = 6371229.0;
constexpr double SPREAD = 0.5;
constexpr double REL_TOL = 1.0e-12;
constexpr int MPS_PER_VTX = 2;
constexpr int MAX_PRINT = 10;

int cellGroup(const int elm) { return elm % 4; }

KOKKOS_INLINE_FUNCTION
double baseValue(const int field, const int cat, const int elm)
{
    return field * 1.0e5 + (cat + 1) * 1.0e3 + elm;
}

KOKKOS_INLINE_FUNCTION
double cellToMPValue(const int cat)
{
    return 1.5 * (cat + 1);
}

// Spherical version of the 10-cell test mesh: same cells and connectivity,
// vertices moved onto a sphere of radius SPHERE_RADIUS.
polyMPO::Mesh* createSphericalTestMesh()
{
    polyMPO::Mesh* planar = polyMPO::initTestMesh(1, 1);

    const int nVertices = planar->getNumVertices();
    const int nCells = planar->getNumElements();

    auto planarCoords = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), planar->getMeshField<polyMPO::MeshF_VtxCoords>());

    polyMPO::MeshFView<polyMPO::MeshF_VtxCoords> vtxCoords("sphericalVtxCoords", nVertices);
    auto vtxCoordsHost = Kokkos::create_mirror_view(vtxCoords);
    for (int v = 0; v < nVertices; v++) {
        const double p[3] = {1.1, planarCoords(v, 0) - 0.5, planarCoords(v, 1) - 0.5};
        const double norm = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        for (int d = 0; d < 3; d++) {
            vtxCoordsHost(v, d) = SPHERE_RADIUS * p[d] / norm;
        }
    }
    Kokkos::deep_copy(vtxCoords, vtxCoordsHost);

    polyMPO::Mesh* mesh = new polyMPO::Mesh(planar->getMeshType(),
                                            polyMPO::geom_spherical_surf,
                                            SPHERE_RADIUS,
                                            nVertices,
                                            nCells,
                                            vtxCoords,
                                            planar->getElm2VtxConn(),
                                            planar->getElm2ElmConn());
    delete planar;
    return mesh;
}

// Cell centers (mean of vertices, projected onto the sphere) and cell areas
// (sum of triangles center-v_i-v_i+1). elm2VtxConn(elm, 0) is the vertex
// count, (elm, 1..n) are 1-based vertex IDs.
void setCellGeometry(polyMPO::Mesh* mesh)
{
    const int nElms = mesh->getNumElements();
    const double radius = mesh->getSphereRadius();
    auto elm2Vtx = mesh->getElm2VtxConn();
    auto vtxCoords = mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elmCenter = mesh->getMeshField<polyMPO::MeshF_ElmCenterXYZ>();
    auto cellArea = mesh->getMeshField<polyMPO::MeshF_CellArea>();

    Kokkos::parallel_for("setTestCellGeometry", nElms, KOKKOS_LAMBDA(const int elm) {
        const int nv = elm2Vtx(elm, 0);

        double c[3] = {0.0, 0.0, 0.0};
        for (int i = 1; i <= nv; i++) {
            const int v = elm2Vtx(elm, i) - 1;
            for (int d = 0; d < 3; d++) c[d] += vtxCoords(v, d);
        }
        const double norm = Kokkos::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
        for (int d = 0; d < 3; d++) c[d] *= radius / norm;

        double area = 0.0;
        for (int i = 1; i <= nv; i++) {
            const int a = elm2Vtx(elm, i) - 1;
            const int b = elm2Vtx(elm, (i % nv) + 1) - 1;
            double u[3], w[3];
            for (int d = 0; d < 3; d++) {
                u[d] = vtxCoords(a, d) - c[d];
                w[d] = vtxCoords(b, d) - c[d];
            }
            const double cx = u[1]*w[2] - u[2]*w[1];
            const double cy = u[2]*w[0] - u[0]*w[2];
            const double cz = u[0]*w[1] - u[1]*w[0];
            area += 0.5 * Kokkos::sqrt(cx*cx + cy*cy + cz*cz);
        }

        for (int d = 0; d < 3; d++) elmCenter(elm, d) = c[d];
        cellArea(elm, 0) = area;
    });
    Kokkos::fence();
}

int expectedMPs(const int elm, const int nVtx)
{
    switch (cellGroup(elm)) {
        case 2:  return 0;
        case 3:  return 1;
        default: return MPS_PER_VTX * nVtx;
    }
}

// MPs on the sphere, laid out per group (see top of file).
polyMPO::MaterialPoints* createTestMPs(polyMPO::Mesh* mesh, std::vector<int>& mpsPerCell)
{
    const int nElms = mesh->getNumElements();
    const double radius = mesh->getSphereRadius();

    auto elm2VtxHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mesh->getElm2VtxConn());
    auto vtxHost = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), mesh->getMeshField<polyMPO::MeshF_VtxCoords>());

    std::vector<double> px, py, pz;
    std::vector<int> owner;
    mpsPerCell.assign(nElms, 0);

    auto addMP = [&](const int elm, double x, double y, double z) {
        const double norm = std::sqrt(x*x + y*y + z*z);
        px.push_back(radius * x / norm);
        py.push_back(radius * y / norm);
        pz.push_back(radius * z / norm);
        owner.push_back(elm);
        mpsPerCell[elm]++;
    };

    for (int elm = 0; elm < nElms; elm++) {
        const int nv = elm2VtxHost(elm, 0);
        double c[3] = {0.0, 0.0, 0.0};
        for (int i = 1; i <= nv; i++) {
            const int v = elm2VtxHost(elm, i) - 1;
            for (int d = 0; d < 3; d++) c[d] += vtxHost(v, d) / nv;
        }

        const int group = cellGroup(elm);
        if (group == 2) {
            continue;                                   // C: no MPs
        }
        if (group == 3) {                               // D: one MP
            const int v = elm2VtxHost(elm, 1) - 1;
            addMP(elm, 0.5 * (c[0] + vtxHost(v, 0)),
                       0.5 * (c[1] + vtxHost(v, 1)),
                       0.5 * (c[2] + vtxHost(v, 2)));
            continue;
        }
        for (int i = 1; i <= nv; i++) {                 // A, B: spread inside
            const int v = elm2VtxHost(elm, i) - 1;
            for (int m = 1; m <= MPS_PER_VTX; m++) {
                const double t = static_cast<double>(m) / (MPS_PER_VTX + 1);
                addMP(elm, (1.0 - t) * c[0] + t * vtxHost(v, 0),
                           (1.0 - t) * c[1] + t * vtxHost(v, 1),
                           (1.0 - t) * c[2] + t * vtxHost(v, 2));
            }
        }
    }

    const int numMPs = static_cast<int>(owner.size());

    polyMPO::IntView numMPsPerElement("numMPsPerElement", nElms);
    polyMPO::IntView mpToElement("mpToElement", numMPs);
    polyMPO::MPSView<polyMPO::MPF_Cur_Pos_XYZ> positions("mpPositions", numMPs);
    polyMPO::MPSView<polyMPO::MPF_Cur_Pos_Rot_Lat_Lon> latLon("mpLatLon", numMPs);

    auto numMPsPerElementHost = Kokkos::create_mirror_view(numMPsPerElement);
    auto mpToElementHost = Kokkos::create_mirror_view(mpToElement);
    auto positionsHost = Kokkos::create_mirror_view(positions);
    auto latLonHost = Kokkos::create_mirror_view(latLon);

    for (int elm = 0; elm < nElms; elm++) numMPsPerElementHost(elm) = mpsPerCell[elm];
    for (int mp = 0; mp < numMPs; mp++) {
        mpToElementHost(mp) = owner[mp];
        positionsHost(mp, 0) = px[mp];
        positionsHost(mp, 1) = py[mp];
        positionsHost(mp, 2) = pz[mp];
        latLonHost(mp, 0) = std::asin(pz[mp] / radius);
        latLonHost(mp, 1) = std::atan2(py[mp], px[mp]);
    }
    Kokkos::deep_copy(numMPsPerElement, numMPsPerElementHost);
    Kokkos::deep_copy(mpToElement, mpToElementHost);
    Kokkos::deep_copy(positions, positionsHost);
    Kokkos::deep_copy(latLon, latLonHost);

    auto p_MPs = new polyMPO::MaterialPoints(nElms, numMPs, positions, numMPsPerElement, mpToElement);

    mesh->setRotatedFlag(false);
    auto mpLatLonField = p_MPs->getData<polyMPO::MPF_Cur_Pos_Rot_Lat_Lon>();
    auto setLatLon = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        mpLatLonField(mp, 0) = latLon(mp, 0);
        mpLatLonField(mp, 1) = latLon(mp, 1);
    };
    p_MPs->parallel_for(setLatLon, "setTestMPLatLon");

    return p_MPs;
}

// Set MP values for one category field (Part 1 pattern).
template <polyMPO::MaterialPointSlice mpSlice>
void setMPValues(polyMPO::MPMesh& mpMesh, const int field)
{
    auto p_MPs = mpMesh.p_MPs;
    auto mpField = p_MPs->getData<mpSlice>();

    auto setValues = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if (mask) {
            const bool addSpread = (elm % 4 == 1) && (mp % 2 == 1);
            for (int k = 0; k < nIceCategories; k++) {
                mpField(mp, k) = baseValue(field, k, elm) + (addSpread ? SPREAD : 0.0);
            }
        }
    };
    p_MPs->parallel_for(setValues, "setCategoryMPValues");
    Kokkos::fence();
}

// Check cell values of one category field after MPs -> cells.
template <polyMPO::MeshFieldIndex mfIndex>
int checkCellValues(polyMPO::MPMesh& mpMesh, const int field, const char* name, const int rank)
{
    auto meshField = mpMesh.p_mesh->getMeshField<mfIndex>();
    auto fieldHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), meshField);

    const int nElms = mpMesh.p_mesh->getNumElements();
    int failures = 0;

    for (int elm = 0; elm < nElms; elm++) {
        for (int k = 0; k < nIceCategories; k++) {

            const double value = fieldHost(elm, k);
            const double base = baseValue(field, k, elm);
            const double tol = REL_TOL * std::abs(base);

            bool ok = std::isfinite(value);
            const char* expected = "";

            switch (cellGroup(elm)) {
                case 0:
                    ok = ok && (std::abs(value - base) <= tol);
                    expected = "exactly base (group A)";
                    break;
                case 1:
                    ok = ok && (value >= base - tol) && (value <= base + SPREAD + tol);
                    expected = "within [base, base+0.5] (group B)";
                    break;
                case 2:
                    ok = ok && (value == 0.0);
                    expected = "0 (group C, no MPs)";
                    break;
                default:
                    ok = ok && (std::abs(value - base) <= tol);
                    expected = "exactly base (group D, 1 MP)";
                    break;
            }

            if (!ok) {
                ++failures;
                if (failures <= MAX_PRINT) {
                    std::cerr
                        << "Rank " << rank << ": " << name
                        << " cell " << elm << " cat " << k
                        << " = " << value
                        << ", base = " << base
                        << ", expected " << expected
                        << std::endl;
                }
            }
        }
    }
    return failures;
}

// Cells -> MPs for one category field, then check every MP.
template <polyMPO::MeshFieldIndex mfIndex, polyMPO::MaterialPointSlice mpSlice>
int checkCellsToMPs(polyMPO::MPMesh& mpMesh, const char* name, const int rank)
{
    auto p_mesh = mpMesh.p_mesh;
    auto p_MPs = mpMesh.p_MPs;
    const int nElms = p_mesh->getNumElements();

    auto meshField = p_mesh->getMeshField<mfIndex>();
    Kokkos::parallel_for("setConstantCellValues", nElms, KOKKOS_LAMBDA(const int elm) {
        for (int k = 0; k < nIceCategories; k++) {
            meshField(elm, k) = cellToMPValue(k);
        }
    });
    Kokkos::fence();

    mpMesh.mapCellsToMPs<mfIndex>();
    Kokkos::fence();

    auto mpField = p_MPs->getData<mpSlice>();
    Kokkos::View<int> badValues("badMPValues");

    auto checkValues = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if (mask) {
            for (int k = 0; k < nIceCategories; k++) {
                const double value = mpField(mp, k);
                const double expected = cellToMPValue(k);
                if (!(Kokkos::fabs(value - expected) <= REL_TOL * expected)) {
                    Kokkos::atomic_increment(&badValues());
                }
            }
        }
    };
    p_MPs->parallel_for(checkValues, "checkCellToMPValues");
    Kokkos::fence();

    auto badHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), badValues);
    if (badHost() != 0) {
        std::cerr
            << "Rank " << rank << ": " << name
            << " cells -> MPs: " << badHost()
            << " MP values differ from 1.5*(k+1)"
            << std::endl;
    }
    return badHost();
}

} // namespace

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
                << "Particle-cell mapping test running with "
                << size << " MPI ranks."
                << std::endl;
        }

        // Spherical 10-cell mesh with cell centers and areas.
        // The mapping is local to each rank, so every rank runs the same test.

        polyMPO::Mesh* mesh = createSphericalTestMesh();
        setCellGeometry(mesh);

        std::vector<int> mpsPerCell;
        polyMPO::MaterialPoints* p_MPs = createTestMPs(mesh, mpsPerCell);

        polyMPO::MPMesh mpMesh(mesh, p_MPs);
        mpMesh.p_MPs->setMPIComm(MPI_COMM_WORLD);

        // mapCellsToMPs uses the gnomonic projection of cells and MPs.
        mesh->setGnomonicProjection(mesh->getRotatedFlag());

        const int nElms = mesh->getNumElements();
        const bool spherical = (mesh->getGeomType() == polyMPO::geom_spherical_surf);

        std::cout
            << "Rank " << rank
            << ": geometry = " << (spherical ? "spherical" : "NOT spherical")
            << ", radius = " << mesh->getSphereRadius()
            << ", cells = " << nElms
            << ", MPs = " << mpMesh.p_MPs->getCount()
            << std::endl;

        // MPs -> cells weights each MP by its area: use uniform area.
        auto mpArea = mpMesh.p_MPs->getData<polyMPO::MPF_Area>();
        auto setArea = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
            if (mask) {
                mpArea(mp, 0) = 1.0;
            }
        };
        mpMesh.p_MPs->parallel_for(setArea, "setUniformMPArea");

        int localFailures = 0;

        if (!spherical) {
            std::cerr << "Rank " << rank << ": mesh is not spherical" << std::endl;
            ++localFailures;
        }

        // Part 0: MP count per cell matches the requested layout.

        Kokkos::View<int*> nMPs("nMPsPerCell", nElms);
        auto countMPs = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
            if (mask) {
                Kokkos::atomic_increment(&nMPs(elm));
            }
        };
        mpMesh.p_MPs->parallel_for(countMPs, "countMPsPerCell");
        Kokkos::fence();
        auto nMPsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nMPs);
        auto elm2VtxHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), mesh->getElm2VtxConn());

        for (int elm = 0; elm < nElms; elm++) {
            const int expected = expectedMPs(elm, elm2VtxHost(elm, 0));
            std::cout
                << "Rank " << rank << ": cell " << elm
                << " group " << "ABCD"[cellGroup(elm)]
                << ": MPs = " << nMPsHost(elm)
                << " (expected " << expected << ")"
                << std::endl;
            if (nMPsHost(elm) != expected) {
                ++localFailures;
            }
        }

        // Part 1: MPs -> cells.
        // Set and map all three fields first, then read all three back.

        setMPValues<polyMPO::MPF_IceAreaCategory>(mpMesh, 1);
        mpMesh.mapMPsToCells<polyMPO::MPF_IceAreaCategory>();

        setMPValues<polyMPO::MPF_IceVolumeCategory>(mpMesh, 2);
        mpMesh.mapMPsToCells<polyMPO::MPF_IceVolumeCategory>();

        setMPValues<polyMPO::MPF_SnowVolumeCategory>(mpMesh, 3);
        mpMesh.mapMPsToCells<polyMPO::MPF_SnowVolumeCategory>();

        Kokkos::fence();

        localFailures += checkCellValues<polyMPO::MeshF_IceAreaCategory>(
            mpMesh, 1, "iceAreaCategory", rank);
        localFailures += checkCellValues<polyMPO::MeshF_IceVolumeCategory>(
            mpMesh, 2, "iceVolumeCategory", rank);
        localFailures += checkCellValues<polyMPO::MeshF_SnowVolumeCategory>(
            mpMesh, 3, "snowVolumeCategory", rank);

        // Part 2: cells -> MPs.

        localFailures += checkCellsToMPs<polyMPO::MeshF_IceAreaCategory,
                                         polyMPO::MPF_IceAreaCategory>(mpMesh, "iceAreaCategory", rank);
        localFailures += checkCellsToMPs<polyMPO::MeshF_IceVolumeCategory,
                                         polyMPO::MPF_IceVolumeCategory>(mpMesh, "iceVolumeCategory", rank);
        localFailures += checkCellsToMPs<polyMPO::MeshF_SnowVolumeCategory,
                                         polyMPO::MPF_SnowVolumeCategory>(mpMesh, "snowVolumeCategory", rank);

        int globalFailures = 0;

        MPI_Allreduce(&localFailures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (rank == 0) {

            if (globalFailures == 0) {
                std::cout
                    << "Particle-cell mapping test PASSED."
                    << std::endl;
            }
            else {
                std::cerr
                    << "Particle-cell mapping test FAILED with "
                    << globalFailures
                    << " errors."
                    << std::endl;
            }
        }

        if (globalFailures != 0) {
            testResult = 1;
        }

        // Do not delete mesh or MPs here.
        // mpMesh owns the Mesh and MaterialPoints objects.
    }

    Kokkos::finalize();
    MPI_Finalize();

    return testResult;
}
