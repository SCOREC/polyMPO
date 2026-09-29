#include "pmpo_MPMesh.hpp"

#include <mpi.h>
#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#ifdef POLYMPO_HAS_NETCDF
#include <netcdf.h>
#endif

// Unit test for particle-cell mapping of the category fields on an MPAS
// spherical centroidal Voronoi (SCVT) mesh:
//   MPF_IceAreaCategory    <-> MeshF_IceAreaCategory
//   MPF_IceVolumeCategory  <-> MeshF_IceVolumeCategory
//   MPF_SnowVolumeCategory <-> MeshF_SnowVolumeCategory
// Each field has nIceCategories components.
//
// Mesh: read from an MPAS mesh file (test/sample_mpas_meshes/
//       spherical_cvt_642elms.nc): vertex coordinates, cell centers
//       (SCVT generators), verticesOnCell and cellsOnCell. Cell areas are
//       computed from the vertices.
//
// MPs: built by this test so the number of MPs per cell is controlled.
// Cells are split into four groups by (cell mod 4):
//
//   Group A (0): 2 MPs per vertex, spread inside the cell, all carry base
//   Group B (1): 2 MPs per vertex, alternating base and base + 0.5
//   Group C (2): no MPs
//   Group D (3): exactly 1 MP, carries base
//
//   base(f, k, cell) = f * 1e7 + (k + 1) * 1e4 + cell
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
//   2a: every cell set to 1.5 * (k + 1)
//       -> every MP exactly 1.5 * (k + 1)
//   2b: every cell set to (k + 1) * 1e4 + cell (varies between neighbours)
//       -> every MP within the min/max of its cell and neighbour cells
//          (gradient limiter, uses the real SCVT neighbours)
//
// Usage: ./testParticleCellMap <path to spherical MPAS mesh .nc file>

namespace {

constexpr double SPREAD = 0.5;
constexpr double REL_TOL = 1.0e-12;
constexpr int MPS_PER_VTX = 2;
constexpr int MAX_PRINT = 10;

int cellGroup(const int elm) { return elm % 4; }

KOKKOS_INLINE_FUNCTION
double baseValue(const int field, const int cat, const int elm)
{
    return field * 1.0e7 + (cat + 1) * 1.0e4 + elm;
}

KOKKOS_INLINE_FUNCTION
double constantCellValue(const int cat)
{
    return 1.5 * (cat + 1);
}

KOKKOS_INLINE_FUNCTION
double varyingCellValue(const int cat, const int elm)
{
    return (cat + 1) * 1.0e4 + elm;
}

//--------------------------------------------------------------------------
// MPAS mesh file
//--------------------------------------------------------------------------

struct MPASMesh {
    int nCells = 0;
    int nVertices = 0;
    int maxEdges = 0;
    double sphereRadius = 0.0;
    std::vector<double> xVertex, yVertex, zVertex;
    std::vector<double> xCell, yCell, zCell;
    std::vector<int> nEdgesOnCell;
    std::vector<int> verticesOnCell;   // [nCells][maxEdges], 1-based
    std::vector<int> cellsOnCell;      // [nCells][maxEdges], 1-based
};

#ifdef POLYMPO_HAS_NETCDF

void ncCheck(const int status, const std::string& what)
{
    if (status != NC_NOERR) {
        std::cerr << "NetCDF error (" << what << "): " << nc_strerror(status) << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

int ncDimLen(const int ncid, const char* name)
{
    int dimid;
    size_t len;
    ncCheck(nc_inq_dimid(ncid, name, &dimid), name);
    ncCheck(nc_inq_dimlen(ncid, dimid, &len), name);
    return static_cast<int>(len);
}

void ncReadDouble(const int ncid, const char* name, std::vector<double>& data, const size_t size)
{
    int varid;
    data.resize(size);
    ncCheck(nc_inq_varid(ncid, name, &varid), name);
    ncCheck(nc_get_var_double(ncid, varid, data.data()), name);
}

void ncReadInt(const int ncid, const char* name, std::vector<int>& data, const size_t size)
{
    int varid;
    data.resize(size);
    ncCheck(nc_inq_varid(ncid, name, &varid), name);
    ncCheck(nc_get_var_int(ncid, varid, data.data()), name);
}

MPASMesh readMPASMesh(const std::string& filename)
{
    MPASMesh m;
    int ncid;
    ncCheck(nc_open(filename.c_str(), NC_NOWRITE, &ncid), filename);

    size_t attLen = 0;
    ncCheck(nc_inq_attlen(ncid, NC_GLOBAL, "on_a_sphere", &attLen), "on_a_sphere");
    std::string onSphere(attLen, ' ');
    ncCheck(nc_get_att_text(ncid, NC_GLOBAL, "on_a_sphere", &onSphere[0]), "on_a_sphere");
    if (onSphere.compare(0, 3, "YES") != 0) {
        std::cerr << "Mesh file is not spherical (on_a_sphere = " << onSphere << ")" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ncCheck(nc_get_att_double(ncid, NC_GLOBAL, "sphere_radius", &m.sphereRadius), "sphere_radius");

    m.nCells = ncDimLen(ncid, "nCells");
    m.nVertices = ncDimLen(ncid, "nVertices");
    m.maxEdges = ncDimLen(ncid, "maxEdges");

    ncReadDouble(ncid, "xVertex", m.xVertex, m.nVertices);
    ncReadDouble(ncid, "yVertex", m.yVertex, m.nVertices);
    ncReadDouble(ncid, "zVertex", m.zVertex, m.nVertices);
    ncReadDouble(ncid, "xCell", m.xCell, m.nCells);
    ncReadDouble(ncid, "yCell", m.yCell, m.nCells);
    ncReadDouble(ncid, "zCell", m.zCell, m.nCells);
    ncReadInt(ncid, "nEdgesOnCell", m.nEdgesOnCell, m.nCells);
    ncReadInt(ncid, "verticesOnCell", m.verticesOnCell, static_cast<size_t>(m.nCells) * m.maxEdges);
    ncReadInt(ncid, "cellsOnCell", m.cellsOnCell, static_cast<size_t>(m.nCells) * m.maxEdges);

    ncCheck(nc_close(ncid), "close");
    return m;
}

#endif

//--------------------------------------------------------------------------
// polyMPO mesh from the MPAS mesh
//--------------------------------------------------------------------------

polyMPO::Mesh* createMesh(const MPASMesh& m)
{
    polyMPO::MeshFView<polyMPO::MeshF_VtxCoords> vtxCoords("scvtVtxCoords", m.nVertices);
    auto vtxCoordsHost = Kokkos::create_mirror_view(vtxCoords);
    for (int v = 0; v < m.nVertices; v++) {
        vtxCoordsHost(v, 0) = m.xVertex[v];
        vtxCoordsHost(v, 1) = m.yVertex[v];
        vtxCoordsHost(v, 2) = m.zVertex[v];
    }
    Kokkos::deep_copy(vtxCoords, vtxCoordsHost);

    polyMPO::IntVtx2ElmView elm2Vtx("scvtElm2VtxConn", m.nCells);
    polyMPO::IntElm2ElmView elm2Elm("scvtElm2ElmConn", m.nCells);
    auto elm2VtxHost = Kokkos::create_mirror_view(elm2Vtx);
    auto elm2ElmHost = Kokkos::create_mirror_view(elm2Elm);

    const int maxConn = static_cast<int>(elm2VtxHost.extent(1)) - 1;
    for (int c = 0; c < m.nCells; c++) {
        const int nv = m.nEdgesOnCell[c];
        if (nv > maxConn) {
            std::cerr << "Cell " << c << " has " << nv << " edges, polyMPO supports "
                      << maxConn << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        elm2VtxHost(c, 0) = nv;
        elm2ElmHost(c, 0) = nv;
        for (int j = 0; j < nv; j++) {
            elm2VtxHost(c, j + 1) = m.verticesOnCell[c * m.maxEdges + j];
            elm2ElmHost(c, j + 1) = m.cellsOnCell[c * m.maxEdges + j];
        }
    }
    Kokkos::deep_copy(elm2Vtx, elm2VtxHost);
    Kokkos::deep_copy(elm2Elm, elm2ElmHost);

    return new polyMPO::Mesh(polyMPO::mesh_general_polygonal,
                             polyMPO::geom_spherical_surf,
                             m.sphereRadius,
                             m.nVertices,
                             m.nCells,
                             vtxCoords,
                             elm2Vtx,
                             elm2Elm);
}

// Cell centers from the file (SCVT generators) and cell areas computed as
// the sum of triangles (center, v_i, v_i+1).
void setCellGeometry(polyMPO::Mesh* mesh, const MPASMesh& m)
{
    const int nElms = mesh->getNumElements();
    auto elm2Vtx = mesh->getElm2VtxConn();
    auto vtxCoords = mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elmCenter = mesh->getMeshField<polyMPO::MeshF_ElmCenterXYZ>();
    auto cellArea = mesh->getMeshField<polyMPO::MeshF_CellArea>();

    auto elmCenterHost = Kokkos::create_mirror_view(elmCenter);
    for (int c = 0; c < nElms; c++) {
        elmCenterHost(c, 0) = m.xCell[c];
        elmCenterHost(c, 1) = m.yCell[c];
        elmCenterHost(c, 2) = m.zCell[c];
    }
    Kokkos::deep_copy(elmCenter, elmCenterHost);

    Kokkos::parallel_for("setScvtCellArea", nElms, KOKKOS_LAMBDA(const int elm) {
        const int nv = elm2Vtx(elm, 0);
        double area = 0.0;
        for (int i = 1; i <= nv; i++) {
            const int a = elm2Vtx(elm, i) - 1;
            const int b = elm2Vtx(elm, (i % nv) + 1) - 1;
            double u[3], w[3];
            for (int d = 0; d < 3; d++) {
                u[d] = vtxCoords(a, d) - elmCenter(elm, d);
                w[d] = vtxCoords(b, d) - elmCenter(elm, d);
            }
            const double cx = u[1]*w[2] - u[2]*w[1];
            const double cy = u[2]*w[0] - u[0]*w[2];
            const double cz = u[0]*w[1] - u[1]*w[0];
            area += 0.5 * Kokkos::sqrt(cx*cx + cy*cy + cz*cz);
        }
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
polyMPO::MaterialPoints* createTestMPs(polyMPO::Mesh* mesh, const MPASMesh& m)
{
    const int nElms = m.nCells;
    const double radius = m.sphereRadius;

    std::vector<double> px, py, pz;
    std::vector<int> owner;
    std::vector<int> mpsPerCell(nElms, 0);

    auto addMP = [&](const int elm, double x, double y, double z) {
        const double norm = std::sqrt(x*x + y*y + z*z);
        px.push_back(radius * x / norm);
        py.push_back(radius * y / norm);
        pz.push_back(radius * z / norm);
        owner.push_back(elm);
        mpsPerCell[elm]++;
    };

    for (int elm = 0; elm < nElms; elm++) {
        const double c[3] = {m.xCell[elm], m.yCell[elm], m.zCell[elm]};
        const int nv = m.nEdgesOnCell[elm];
        const int group = cellGroup(elm);

        if (group == 2) {
            continue;                                           // C: no MPs
        }
        if (group == 3) {                                       // D: one MP
            const int v = m.verticesOnCell[elm * m.maxEdges] - 1;
            addMP(elm, 0.5 * (c[0] + m.xVertex[v]),
                       0.5 * (c[1] + m.yVertex[v]),
                       0.5 * (c[2] + m.zVertex[v]));
            continue;
        }
        for (int j = 0; j < nv; j++) {                          // A, B: spread inside
            const int v = m.verticesOnCell[elm * m.maxEdges + j] - 1;
            for (int s = 1; s <= MPS_PER_VTX; s++) {
                const double t = static_cast<double>(s) / (MPS_PER_VTX + 1);
                addMP(elm, (1.0 - t) * c[0] + t * m.xVertex[v],
                           (1.0 - t) * c[1] + t * m.yVertex[v],
                           (1.0 - t) * c[2] + t * m.zVertex[v]);
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

//--------------------------------------------------------------------------
// Part 1: MPs -> cells
//--------------------------------------------------------------------------

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

template <polyMPO::MeshFieldIndex mfIndex>
int checkCellValues(polyMPO::MPMesh& mpMesh, const int field, const char* name, const int rank)
{
    auto meshField = mpMesh.p_mesh->getMeshField<mfIndex>();
    auto fieldHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), meshField);

    const int nElms = mpMesh.p_mesh->getNumElements();
    int failures = 0;
    int groupFailures[4] = {0, 0, 0, 0};

    for (int elm = 0; elm < nElms; elm++) {
        for (int k = 0; k < nIceCategories; k++) {

            const double value = fieldHost(elm, k);
            const double base = baseValue(field, k, elm);
            const double tol = REL_TOL * std::abs(base);
            const int group = cellGroup(elm);

            bool ok = std::isfinite(value);
            const char* expected = "";

            switch (group) {
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
                ++groupFailures[group];
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

    std::cout
        << "Rank " << rank << ": MPs -> cells " << name
        << ": failures A/B/C/D = "
        << groupFailures[0] << "/" << groupFailures[1] << "/"
        << groupFailures[2] << "/" << groupFailures[3]
        << std::endl;

    return failures;
}

//--------------------------------------------------------------------------
// Part 2: cells -> MPs
//--------------------------------------------------------------------------

// Min/max of each cell and its valid neighbours for the varying cell field.
void neighbourBounds(const MPASMesh& m,
                     Kokkos::View<double**>& lo, Kokkos::View<double**>& hi)
{
    const int nElms = m.nCells;
    lo = Kokkos::View<double**>("limiterLo", nElms, nIceCategories);
    hi = Kokkos::View<double**>("limiterHi", nElms, nIceCategories);
    auto loHost = Kokkos::create_mirror_view(lo);
    auto hiHost = Kokkos::create_mirror_view(hi);

    for (int elm = 0; elm < nElms; elm++) {
        for (int k = 0; k < nIceCategories; k++) {
            double vmin = varyingCellValue(k, elm);
            double vmax = vmin;
            for (int j = 0; j < m.nEdgesOnCell[elm]; j++) {
                const int nb = m.cellsOnCell[elm * m.maxEdges + j] - 1;
                if (nb < 0 || nb >= nElms) continue;
                const double v = varyingCellValue(k, nb);
                vmin = std::min(vmin, v);
                vmax = std::max(vmax, v);
            }
            loHost(elm, k) = vmin;
            hiHost(elm, k) = vmax;
        }
    }
    Kokkos::deep_copy(lo, loHost);
    Kokkos::deep_copy(hi, hiHost);
}

template <polyMPO::MeshFieldIndex mfIndex, polyMPO::MaterialPointSlice mpSlice>
int checkCellsToMPs(polyMPO::MPMesh& mpMesh, const char* name, const int rank,
                    const Kokkos::View<double**>& lo, const Kokkos::View<double**>& hi)
{
    auto p_mesh = mpMesh.p_mesh;
    auto p_MPs = mpMesh.p_MPs;
    const int nElms = p_mesh->getNumElements();
    auto meshField = p_mesh->getMeshField<mfIndex>();
    auto mpField = p_MPs->getData<mpSlice>();

    // 2a: constant cell field -> every MP exactly that value
    Kokkos::parallel_for("setConstantCellValues", nElms, KOKKOS_LAMBDA(const int elm) {
        for (int k = 0; k < nIceCategories; k++) {
            meshField(elm, k) = constantCellValue(k);
        }
    });
    Kokkos::fence();
    mpMesh.mapCellsToMPs<mfIndex>();
    Kokkos::fence();

    Kokkos::View<int> badConstant("badConstant");
    auto checkConstant = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if (mask) {
            for (int k = 0; k < nIceCategories; k++) {
                const double expected = constantCellValue(k);
                if (!(Kokkos::fabs(mpField(mp, k) - expected) <= REL_TOL * expected)) {
                    Kokkos::atomic_increment(&badConstant());
                }
            }
        }
    };
    p_MPs->parallel_for(checkConstant, "checkConstantCellToMP");

    // 2b: varying cell field -> every MP within its cell/neighbour min/max
    Kokkos::parallel_for("setVaryingCellValues", nElms, KOKKOS_LAMBDA(const int elm) {
        for (int k = 0; k < nIceCategories; k++) {
            meshField(elm, k) = varyingCellValue(k, elm);
        }
    });
    Kokkos::fence();
    mpMesh.mapCellsToMPs<mfIndex>();
    Kokkos::fence();

    Kokkos::View<int> badBounded("badBounded");
    auto checkBounded = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if (mask) {
            for (int k = 0; k < nIceCategories; k++) {
                const double v = mpField(mp, k);
                const double tol = REL_TOL * Kokkos::fabs(hi(elm, k));
                if (!(v >= lo(elm, k) - tol && v <= hi(elm, k) + tol)) {
                    Kokkos::atomic_increment(&badBounded());
                }
            }
        }
    };
    p_MPs->parallel_for(checkBounded, "checkBoundedCellToMP");
    Kokkos::fence();

    auto badConstantHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), badConstant);
    auto badBoundedHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), badBounded);

    std::cout
        << "Rank " << rank << ": cells -> MPs " << name
        << ": constant failures = " << badConstantHost()
        << ", out-of-bounds (limiter) failures = " << badBoundedHost()
        << std::endl;

    return badConstantHost() + badBoundedHost();
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

#ifndef POLYMPO_HAS_NETCDF
        if (rank == 0) {
            std::cerr << "testParticleCellMap requires NetCDF; skipping." << std::endl;
        }
        testResult = 77;
#else
        if (argc < 2) {
            if (rank == 0) {
                std::cerr << "Usage: " << argv[0] << " <path to spherical MPAS mesh .nc file>" << std::endl;
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        const std::string meshFile = argv[1];

        if (rank == 0) {
            std::cout
                << "Particle-cell mapping test running with "
                << size << " MPI ranks on " << meshFile
                << std::endl;
        }

        // SCVT mesh from the MPAS mesh file.
        // The mapping is local to each rank, so every rank runs the same test.

        const MPASMesh mpasMesh = readMPASMesh(meshFile);

        polyMPO::Mesh* mesh = createMesh(mpasMesh);
        setCellGeometry(mesh, mpasMesh);

        polyMPO::MaterialPoints* p_MPs = createTestMPs(mesh, mpasMesh);

        polyMPO::MPMesh mpMesh(mesh, p_MPs);
        mpMesh.p_MPs->setMPIComm(MPI_COMM_WORLD);

        // mapCellsToMPs uses the gnomonic projection of cells and MPs.
        mesh->setGnomonicProjection(mesh->getRotatedFlag());

        const int nElms = mesh->getNumElements();
        const bool spherical = (mesh->getGeomType() == polyMPO::geom_spherical_surf);

        std::cout
            << "Rank " << rank
            << ": geometry = " << (spherical ? "spherical (SCVT)" : "NOT spherical")
            << ", radius = " << mesh->getSphereRadius()
            << ", cells = " << nElms
            << ", vertices = " << mesh->getNumVertices()
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

        int cellsPerGroup[4] = {0, 0, 0, 0};
        int countFailures = 0;
        for (int elm = 0; elm < nElms; elm++) {
            ++cellsPerGroup[cellGroup(elm)];
            if (nMPsHost(elm) != expectedMPs(elm, mpasMesh.nEdgesOnCell[elm])) {
                ++countFailures;
            }
        }
        std::cout
            << "Rank " << rank
            << ": cells per group A/B/C/D = "
            << cellsPerGroup[0] << "/" << cellsPerGroup[1] << "/"
            << cellsPerGroup[2] << "/" << cellsPerGroup[3]
            << ", MP count mismatches = " << countFailures
            << std::endl;
        localFailures += countFailures;

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

        Kokkos::View<double**> lo, hi;
        neighbourBounds(mpasMesh, lo, hi);

        localFailures += checkCellsToMPs<polyMPO::MeshF_IceAreaCategory,
                                         polyMPO::MPF_IceAreaCategory>(mpMesh, "iceAreaCategory", rank, lo, hi);
        localFailures += checkCellsToMPs<polyMPO::MeshF_IceVolumeCategory,
                                         polyMPO::MPF_IceVolumeCategory>(mpMesh, "iceVolumeCategory", rank, lo, hi);
        localFailures += checkCellsToMPs<polyMPO::MeshF_SnowVolumeCategory,
                                         polyMPO::MPF_SnowVolumeCategory>(mpMesh, "snowVolumeCategory", rank, lo, hi);

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
#endif
    }

    Kokkos::finalize();
    MPI_Finalize();

    return testResult;
}
