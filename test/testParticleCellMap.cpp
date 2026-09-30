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

// Partitioned (4-rank) unit test for particle-cell mapping on spherical centroidal Voronoi (SCVT) mesh.
//
// Partition: the SCVT mesh (spherical_cvt_642elms.nc) is split into 4 longitude wedges; rank r owns the cells whose center lies in wedge r.
//
// MP layout and values (based on the global cell ID, so the partitioned and serial runs see identical inputs):
//   Group A (gid mod 4 = 0): 2 MPs per vertex, all carry base
//   Group B (gid mod 4 = 1): 2 MPs per vertex, base or base + 0.5
//                            (+0.5 for MPs on one side of the cell center)
//   Group C (gid mod 4 = 2): no MPs
//   Group D (gid mod 4 = 3): exactly 1 MP, carries base
//   base(f, k, gid) = f * 1e7 + (k + 1) * 1e4 + gid
//
//   Two important parts of this test:
//   Part 1: MPs -> cells: every owned cell value equals the serial value.
//   Part 2: cells -> MPs: cells (owned + halo) set to
//           (k + 1) * 1e4 + gid (varies between neighbours); the MPs of
//           every owned cell must match the serial MPs of that cell
//           (sum, min and max per category).
// Halo cells get no MP contribution (MPs live on their owning rank); the test reports them but does not require a value.

namespace {

constexpr int NUM_PARTS = 4;
constexpr int NUM_FIELDS = 3;
constexpr double SPREAD = 0.5;
constexpr double REL_TOL = 1.0e-12;
constexpr int MPS_PER_VTX = 2;
constexpr int MAX_PRINT = 10;

KOKKOS_INLINE_FUNCTION
int cellGroup(const int gid) { return gid % 4; }

KOKKOS_INLINE_FUNCTION
double baseValue(const int field, const int cat, const int gid)
{
    return field * 1.0e7 + (cat + 1) * 1.0e4 + gid;
}

KOKKOS_INLINE_FUNCTION
double varyingCellValue(const int cat, const int gid)
{
    return (cat + 1) * 1.0e4 + gid;
}

bool nearlyEqual(const double a, const double b)
{
    return std::isfinite(a) && std::isfinite(b) &&
           std::abs(a - b) <= REL_TOL * std::max(1.0, std::abs(b));
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
// Local (rank) mesh: owned cells first, then halo cells
//--------------------------------------------------------------------------

struct LocalMesh {
    int nCells = 0;
    int nOwned = 0;
    int nVertices = 0;
    int maxEdges = 0;
    double sphereRadius = 0.0;
    std::vector<double> xVertex, yVertex, zVertex;
    std::vector<double> xCell, yCell, zCell;
    std::vector<int> nEdgesOnCell;
    std::vector<int> verticesOnCell;   // [nCells][maxEdges], 1-based local
    std::vector<int> cellsOnCell;      // [nCells][maxEdges], 1-based local, invalid = nCells + 1
    std::vector<int> globalId;         // local cell -> global cell
    std::vector<bool> isBoundary;      // owned cell with a neighbour owned by another rank
};

// owner[g] = owning rank of global cell g. serial = true: all cells owned.
LocalMesh buildLocalMesh(const MPASMesh& m, const std::vector<int>& owner,
                         const int rank, const bool serial)
{
    LocalMesh L;
    L.maxEdges = m.maxEdges;
    L.sphereRadius = m.sphereRadius;

    std::vector<int> g2l(m.nCells, -1);
    std::vector<int> cells;

    for (int g = 0; g < m.nCells; g++) {
        if (serial || owner[g] == rank) {
            g2l[g] = static_cast<int>(cells.size());
            cells.push_back(g);
        }
    }
    L.nOwned = static_cast<int>(cells.size());

    if (!serial) {
        for (int i = 0; i < L.nOwned; i++) {
            const int g = cells[i];
            for (int j = 0; j < m.nEdgesOnCell[g]; j++) {
                const int nb = m.cellsOnCell[g * m.maxEdges + j] - 1;
                if (nb < 0 || nb >= m.nCells || g2l[nb] >= 0) continue;
                g2l[nb] = static_cast<int>(cells.size());
                cells.push_back(nb);
            }
        }
    }
    L.nCells = static_cast<int>(cells.size());

    std::vector<int> v2l(m.nVertices, -1);
    L.nEdgesOnCell.resize(L.nCells);
    L.verticesOnCell.assign(static_cast<size_t>(L.nCells) * L.maxEdges, 0);
    L.cellsOnCell.assign(static_cast<size_t>(L.nCells) * L.maxEdges, L.nCells + 1);
    L.isBoundary.assign(L.nCells, false);

    for (int l = 0; l < L.nCells; l++) {
        const int g = cells[l];
        L.globalId.push_back(g);
        L.xCell.push_back(m.xCell[g]);
        L.yCell.push_back(m.yCell[g]);
        L.zCell.push_back(m.zCell[g]);
        L.nEdgesOnCell[l] = m.nEdgesOnCell[g];

        for (int j = 0; j < m.nEdgesOnCell[g]; j++) {
            const int v = m.verticesOnCell[g * m.maxEdges + j] - 1;
            if (v2l[v] < 0) {
                v2l[v] = static_cast<int>(L.xVertex.size());
                L.xVertex.push_back(m.xVertex[v]);
                L.yVertex.push_back(m.yVertex[v]);
                L.zVertex.push_back(m.zVertex[v]);
            }
            L.verticesOnCell[l * L.maxEdges + j] = v2l[v] + 1;

            const int nb = m.cellsOnCell[g * m.maxEdges + j] - 1;
            if (nb >= 0 && nb < m.nCells && g2l[nb] >= 0) {
                L.cellsOnCell[l * L.maxEdges + j] = g2l[nb] + 1;
            }
            if (l < L.nOwned && nb >= 0 && nb < m.nCells && !serial && owner[nb] != rank) {
                L.isBoundary[l] = true;
            }
        }
    }
    L.nVertices = static_cast<int>(L.xVertex.size());
    return L;
}

//--------------------------------------------------------------------------
// polyMPO mesh and MPs from a local mesh
//--------------------------------------------------------------------------

polyMPO::Mesh* createMesh(const LocalMesh& L)
{
    polyMPO::MeshFView<polyMPO::MeshF_VtxCoords> vtxCoords("localVtxCoords", L.nVertices);
    auto vtxCoordsHost = Kokkos::create_mirror_view(vtxCoords);
    for (int v = 0; v < L.nVertices; v++) {
        vtxCoordsHost(v, 0) = L.xVertex[v];
        vtxCoordsHost(v, 1) = L.yVertex[v];
        vtxCoordsHost(v, 2) = L.zVertex[v];
    }
    Kokkos::deep_copy(vtxCoords, vtxCoordsHost);

    polyMPO::IntVtx2ElmView elm2Vtx("localElm2VtxConn", L.nCells);
    polyMPO::IntElm2ElmView elm2Elm("localElm2ElmConn", L.nCells);
    auto elm2VtxHost = Kokkos::create_mirror_view(elm2Vtx);
    auto elm2ElmHost = Kokkos::create_mirror_view(elm2Elm);

    const int maxConn = static_cast<int>(elm2VtxHost.extent(1)) - 1;
    for (int c = 0; c < L.nCells; c++) {
        const int nv = L.nEdgesOnCell[c];
        if (nv > maxConn) {
            std::cerr << "Cell with " << nv << " edges, polyMPO supports " << maxConn << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        elm2VtxHost(c, 0) = nv;
        elm2ElmHost(c, 0) = nv;
        for (int j = 0; j < nv; j++) {
            elm2VtxHost(c, j + 1) = L.verticesOnCell[c * L.maxEdges + j];
            elm2ElmHost(c, j + 1) = L.cellsOnCell[c * L.maxEdges + j];
        }
    }
    Kokkos::deep_copy(elm2Vtx, elm2VtxHost);
    Kokkos::deep_copy(elm2Elm, elm2ElmHost);

    return new polyMPO::Mesh(polyMPO::mesh_general_polygonal,
                             polyMPO::geom_spherical_surf,
                             L.sphereRadius,
                             L.nVertices,
                             L.nCells,
                             vtxCoords,
                             elm2Vtx,
                             elm2Elm);
}

// Cell centers from the file (SCVT generators) and cell areas computed as
// the sum of triangles (center, v_i, v_i+1).
void setCellGeometry(polyMPO::Mesh* mesh, const LocalMesh& L)
{
    const int nElms = mesh->getNumElements();
    auto elm2Vtx = mesh->getElm2VtxConn();
    auto vtxCoords = mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elmCenter = mesh->getMeshField<polyMPO::MeshF_ElmCenterXYZ>();
    auto cellArea = mesh->getMeshField<polyMPO::MeshF_CellArea>();

    auto elmCenterHost = Kokkos::create_mirror_view(elmCenter);
    for (int c = 0; c < nElms; c++) {
        elmCenterHost(c, 0) = L.xCell[c];
        elmCenterHost(c, 1) = L.yCell[c];
        elmCenterHost(c, 2) = L.zCell[c];
    }
    Kokkos::deep_copy(elmCenter, elmCenterHost);

    Kokkos::parallel_for("setLocalCellArea", nElms, KOKKOS_LAMBDA(const int elm) {
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

// MPs only in owned cells, laid out by the global cell ID.
polyMPO::MaterialPoints* createMPs(polyMPO::Mesh* mesh, const LocalMesh& L)
{
    const int nElms = L.nCells;
    const double radius = L.sphereRadius;

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

    for (int elm = 0; elm < L.nOwned; elm++) {
        const double c[3] = {L.xCell[elm], L.yCell[elm], L.zCell[elm]};
        const int nv = L.nEdgesOnCell[elm];
        const int group = cellGroup(L.globalId[elm]);

        if (group == 2) {
            continue;                                           // C: no MPs
        }
        if (group == 3) {                                       // D: one MP
            const int v = L.verticesOnCell[elm * L.maxEdges] - 1;
            addMP(elm, 0.5 * (c[0] + L.xVertex[v]),
                       0.5 * (c[1] + L.yVertex[v]),
                       0.5 * (c[2] + L.zVertex[v]));
            continue;
        }
        for (int j = 0; j < nv; j++) {                          // A, B: spread inside
            const int v = L.verticesOnCell[elm * L.maxEdges + j] - 1;
            for (int s = 1; s <= MPS_PER_VTX; s++) {
                const double t = static_cast<double>(s) / (MPS_PER_VTX + 1);
                addMP(elm, (1.0 - t) * c[0] + t * L.xVertex[v],
                           (1.0 - t) * c[1] + t * L.yVertex[v],
                           (1.0 - t) * c[2] + t * L.zVertex[v]);
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
// One run of the mapping on a local mesh
//--------------------------------------------------------------------------

struct Results {
    int numMPs = 0;
    std::vector<int> mpCount;                    // [nCells]
    std::vector<double> cell[NUM_FIELDS];        // MPs -> cells, [nCells * nIceCategories]
    std::vector<double> mpSum[NUM_FIELDS];       // cells -> MPs, per cell and category
    std::vector<double> mpMin[NUM_FIELDS];
    std::vector<double> mpMax[NUM_FIELDS];
};

template <polyMPO::MaterialPointSlice mpSlice>
void setMPValues(polyMPO::MPMesh& mpMesh, const int field, const polyMPO::IntView& gid)
{
    auto p_MPs = mpMesh.p_MPs;
    auto mpField = p_MPs->getData<mpSlice>();
    auto mpPos = p_MPs->getData<polyMPO::MPF_Cur_Pos_XYZ>();
    auto elmCenter = mpMesh.p_mesh->getMeshField<polyMPO::MeshF_ElmCenterXYZ>();

    auto setValues = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if (mask) {
            const int g = gid(elm);
            // group B: +0.5 for MPs on one side of the cell center (position
            // based, so it is identical in the serial and partitioned runs)
            const double side = (mpPos(mp, 0) + mpPos(mp, 1) + mpPos(mp, 2))
                              - (elmCenter(elm, 0) + elmCenter(elm, 1) + elmCenter(elm, 2));
            const bool addSpread = (cellGroup(g) == 1) && (side > 0.0);
            for (int k = 0; k < nIceCategories; k++) {
                mpField(mp, k) = baseValue(field, k, g) + (addSpread ? SPREAD : 0.0);
            }
        }
    };
    p_MPs->parallel_for(setValues, "setCategoryMPValues");
    Kokkos::fence();
}

template <polyMPO::MeshFieldIndex mfIndex>
std::vector<double> readCellValues(polyMPO::MPMesh& mpMesh)
{
    auto meshField = mpMesh.p_mesh->getMeshField<mfIndex>();
    auto fieldHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), meshField);
    const int nElms = mpMesh.p_mesh->getNumElements();
    std::vector<double> values(static_cast<size_t>(nElms) * nIceCategories);
    for (int elm = 0; elm < nElms; elm++) {
        for (int k = 0; k < nIceCategories; k++) {
            values[elm * nIceCategories + k] = fieldHost(elm, k);
        }
    }
    return values;
}

template <polyMPO::MeshFieldIndex mfIndex, polyMPO::MaterialPointSlice mpSlice>
void cellsToMPs(polyMPO::MPMesh& mpMesh, const polyMPO::IntView& gid, Results& r, const int f)
{
    auto p_MPs = mpMesh.p_MPs;
    const int nElms = mpMesh.p_mesh->getNumElements();
    auto meshField = mpMesh.p_mesh->getMeshField<mfIndex>();

    // Owned and halo cells, as MPAS provides them.
    Kokkos::parallel_for("setVaryingCellValues", nElms, KOKKOS_LAMBDA(const int elm) {
        for (int k = 0; k < nIceCategories; k++) {
            meshField(elm, k) = varyingCellValue(k, gid(elm));
        }
    });
    Kokkos::fence();
    mpMesh.mapCellsToMPs<mfIndex>();
    Kokkos::fence();

    Kokkos::View<double**> sum("mpSum", nElms, nIceCategories);
    Kokkos::View<double**> vmin("mpMin", nElms, nIceCategories);
    Kokkos::View<double**> vmax("mpMax", nElms, nIceCategories);
    Kokkos::deep_copy(vmin, 1.0e300);
    Kokkos::deep_copy(vmax, -1.0e300);

    auto mpField = p_MPs->getData<mpSlice>();
    auto aggregate = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if (mask) {
            for (int k = 0; k < nIceCategories; k++) {
                const double v = mpField(mp, k);
                Kokkos::atomic_add(&sum(elm, k), v);
                Kokkos::atomic_fetch_min(&vmin(elm, k), v);
                Kokkos::atomic_fetch_max(&vmax(elm, k), v);
            }
        }
    };
    p_MPs->parallel_for(aggregate, "aggregateMPValues");
    Kokkos::fence();

    auto sumHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sum);
    auto minHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vmin);
    auto maxHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vmax);

    const size_t n = static_cast<size_t>(nElms) * nIceCategories;
    r.mpSum[f].resize(n);
    r.mpMin[f].resize(n);
    r.mpMax[f].resize(n);
    for (int elm = 0; elm < nElms; elm++) {
        for (int k = 0; k < nIceCategories; k++) {
            r.mpSum[f][elm * nIceCategories + k] = sumHost(elm, k);
            r.mpMin[f][elm * nIceCategories + k] = minHost(elm, k);
            r.mpMax[f][elm * nIceCategories + k] = maxHost(elm, k);
        }
    }
}

Results runMapping(const LocalMesh& L)
{
    Results r;

    polyMPO::Mesh* mesh = createMesh(L);
    setCellGeometry(mesh, L);
    polyMPO::MaterialPoints* p_MPs = createMPs(mesh, L);

    polyMPO::MPMesh mpMesh(mesh, p_MPs);
    mpMesh.p_MPs->setMPIComm(MPI_COMM_WORLD);
    mesh->setGnomonicProjection(mesh->getRotatedFlag());

    r.numMPs = mpMesh.p_MPs->getCount();

    polyMPO::IntView gid("globalCellId", L.nCells);
    auto gidHost = Kokkos::create_mirror_view(gid);
    for (int l = 0; l < L.nCells; l++) gidHost(l) = L.globalId[l];
    Kokkos::deep_copy(gid, gidHost);

    // uniform MP area
    auto mpArea = mpMesh.p_MPs->getData<polyMPO::MPF_Area>();
    auto setArea = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if (mask) mpArea(mp, 0) = 1.0;
    };
    mpMesh.p_MPs->parallel_for(setArea, "setUniformMPArea");

    // MP count per cell
    Kokkos::View<int*> nMPs("nMPsPerCell", L.nCells);
    auto countMPs = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if (mask) Kokkos::atomic_increment(&nMPs(elm));
    };
    mpMesh.p_MPs->parallel_for(countMPs, "countMPsPerCell");
    Kokkos::fence();
    auto nMPsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nMPs);
    r.mpCount.assign(nMPsHost.data(), nMPsHost.data() + L.nCells);

    // MPs -> cells: map all three fields, then read all three
    setMPValues<polyMPO::MPF_IceAreaCategory>(mpMesh, 1, gid);
    mpMesh.mapMPsToCells<polyMPO::MPF_IceAreaCategory>();
    setMPValues<polyMPO::MPF_IceVolumeCategory>(mpMesh, 2, gid);
    mpMesh.mapMPsToCells<polyMPO::MPF_IceVolumeCategory>();
    setMPValues<polyMPO::MPF_SnowVolumeCategory>(mpMesh, 3, gid);
    mpMesh.mapMPsToCells<polyMPO::MPF_SnowVolumeCategory>();
    Kokkos::fence();

    r.cell[0] = readCellValues<polyMPO::MeshF_IceAreaCategory>(mpMesh);
    r.cell[1] = readCellValues<polyMPO::MeshF_IceVolumeCategory>(mpMesh);
    r.cell[2] = readCellValues<polyMPO::MeshF_SnowVolumeCategory>(mpMesh);

    // cells -> MPs
    cellsToMPs<polyMPO::MeshF_IceAreaCategory, polyMPO::MPF_IceAreaCategory>(mpMesh, gid, r, 0);
    cellsToMPs<polyMPO::MeshF_IceVolumeCategory, polyMPO::MPF_IceVolumeCategory>(mpMesh, gid, r, 1);
    cellsToMPs<polyMPO::MeshF_SnowVolumeCategory, polyMPO::MPF_SnowVolumeCategory>(mpMesh, gid, r, 2);

    // mpMesh owns and deletes the Mesh and MaterialPoints
    return r;
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
            std::cerr << "testParticleCellMapPartitioned requires NetCDF; skipping." << std::endl;
        }
        testResult = 77;
#else
        if (size != NUM_PARTS) {
            if (rank == 0) {
                std::cerr
                    << "This test requires exactly " << NUM_PARTS
                    << " MPI ranks (got " << size << "); skipping."
                    << std::endl;
            }
            testResult = 77;
        }
        else {
      #ifdef PMPO_TEST_SCVT_MESH
              const std::string meshFile = (argc >= 2) ? argv[1] : PMPO_TEST_SCVT_MESH;
      #else
             if (argc < 2) {
                if (rank == 0) std::cerr << "Usage: " << argv[0] << " <mesh.nc>" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            const std::string meshFile = argv[1];
      #endif
 
            if (rank == 0) {
                std::cout
                    << "Partitioned particle-cell mapping test running with "
                    << size << " MPI ranks on " << meshFile
                    << std::endl;
            }

            const MPASMesh m = readMPASMesh(meshFile);

            // Partition: 4 longitude wedges of the cell centers.
            const double pi = 4.0 * std::atan(1.0);
            std::vector<int> owner(m.nCells);
            for (int g = 0; g < m.nCells; g++) {
                const double lon = std::atan2(m.yCell[g], m.xCell[g]);          // [-pi, pi]
                int p = static_cast<int>((lon + pi) / (2.0 * pi) * NUM_PARTS);
                owner[g] = std::min(std::max(p, 0), NUM_PARTS - 1);
            }

            const LocalMesh serialMesh = buildLocalMesh(m, owner, rank, true);
            const LocalMesh localMesh = buildLocalMesh(m, owner, rank, false);

            const Results serial = runMapping(serialMesh);
            const Results part = runMapping(localMesh);

            int nBoundary = 0;
            for (int l = 0; l < localMesh.nOwned; l++) {
                if (localMesh.isBoundary[l]) ++nBoundary;
            }

            std::cout
                << "Rank " << rank
                << ": owned cells = " << localMesh.nOwned
                << " (partition-boundary " << nBoundary << ")"
                << ", halo cells = " << localMesh.nCells - localMesh.nOwned
                << ", MPs = " << part.numMPs
                << std::endl;

            int localFailures = 0;

            // Part 0: totals over all ranks, MP counts per owned cell.
            int ownedTotal = 0;
            int mpsTotal = 0;
            MPI_Allreduce(&localMesh.nOwned, &ownedTotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&part.numMPs, &mpsTotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (rank == 0) {
                std::cout
                    << "All ranks: owned cells = " << ownedTotal << " (mesh " << m.nCells << ")"
                    << ", MPs = " << mpsTotal << " (serial " << serial.numMPs << ")"
                    << std::endl;
                if (ownedTotal != m.nCells) ++localFailures;
                if (mpsTotal != serial.numMPs) ++localFailures;
            }

            int countMismatch = 0;
            for (int l = 0; l < localMesh.nOwned; l++) {
                if (part.mpCount[l] != serial.mpCount[localMesh.globalId[l]]) ++countMismatch;
            }
            localFailures += countMismatch;

            // Part 1: MPs -> cells, owned cells vs serial.
            int p1Interior = 0, p1Boundary = 0, printed = 0;
            for (int f = 0; f < NUM_FIELDS; f++) {
                for (int l = 0; l < localMesh.nOwned; l++) {
                    const int g = localMesh.globalId[l];
                    for (int k = 0; k < nIceCategories; k++) {
                        const double p = part.cell[f][l * nIceCategories + k];
                        const double s = serial.cell[f][g * nIceCategories + k];
                        if (!nearlyEqual(p, s)) {
                            (localMesh.isBoundary[l] ? p1Boundary : p1Interior)++;
                            if (printed++ < MAX_PRINT) {
                                std::cerr << "Rank " << rank << ": MPs -> cells field " << f + 1
                                          << " cell " << g << " cat " << k
                                          << ": partitioned = " << p << ", serial = " << s << std::endl;
                            }
                        }
                    }
                }
            }
            localFailures += p1Interior + p1Boundary;

            // Halo cells: no MPs on this rank (reported only).
            int haloNonZero = 0;
            for (int f = 0; f < NUM_FIELDS; f++) {
                for (int l = localMesh.nOwned; l < localMesh.nCells; l++) {
                    for (int k = 0; k < nIceCategories; k++) {
                        if (part.cell[f][l * nIceCategories + k] != 0.0) ++haloNonZero;
                    }
                }
            }

            // Part 2: cells -> MPs, per owned cell vs serial (sum, min, max).
            int p2Interior = 0, p2Boundary = 0;
            printed = 0;
            for (int f = 0; f < NUM_FIELDS; f++) {
                for (int l = 0; l < localMesh.nOwned; l++) {
                    const int g = localMesh.globalId[l];
                    if (part.mpCount[l] == 0) continue;
                    for (int k = 0; k < nIceCategories; k++) {
                        const size_t il = static_cast<size_t>(l) * nIceCategories + k;
                        const size_t ig = static_cast<size_t>(g) * nIceCategories + k;
                        const bool ok = nearlyEqual(part.mpSum[f][il], serial.mpSum[f][ig]) &&
                                        nearlyEqual(part.mpMin[f][il], serial.mpMin[f][ig]) &&
                                        nearlyEqual(part.mpMax[f][il], serial.mpMax[f][ig]);
                        if (!ok) {
                            (localMesh.isBoundary[l] ? p2Boundary : p2Interior)++;
                            if (printed++ < MAX_PRINT) {
                                std::cerr << "Rank " << rank << ": cells -> MPs field " << f + 1
                                          << " cell " << g << " cat " << k
                                          << ": partitioned sum/min/max = " << part.mpSum[f][il]
                                          << "/" << part.mpMin[f][il] << "/" << part.mpMax[f][il]
                                          << ", serial = " << serial.mpSum[f][ig]
                                          << "/" << serial.mpMin[f][ig] << "/" << serial.mpMax[f][ig]
                                          << std::endl;
                            }
                        }
                    }
                }
            }
            localFailures += p2Interior + p2Boundary;

            std::cout
                << "Rank " << rank
                << ": MP count mismatches = " << countMismatch
                << "; MPs -> cells mismatches (interior/boundary) = " << p1Interior << "/" << p1Boundary
                << "; cells -> MPs mismatches (interior/boundary) = " << p2Interior << "/" << p2Boundary
                << "; halo cell values != 0 (info) = " << haloNonZero
                << std::endl;

            int globalFailures = 0;
            MPI_Allreduce(&localFailures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            if (rank == 0) {
                if (globalFailures == 0) {
                    std::cout << "Partitioned particle-cell mapping test PASSED." << std::endl;
                }
                else {
                    std::cerr << "Partitioned particle-cell mapping test FAILED with "
                              << globalFailures << " errors." << std::endl;
                }
            }
            if (globalFailures != 0) {
                testResult = 1;
            }
        }
#endif
    }

    Kokkos::finalize();
    MPI_Finalize();

    return testResult;
}
