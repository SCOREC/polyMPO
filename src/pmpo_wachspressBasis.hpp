#ifndef POLYMPO_WACHSPRESS_BASIS_H
#define POLYMPO_WACHSPRESS_BASIS_H

#include "pmpo_utils.hpp"
#include "pmpo_MPMesh.hpp"

namespace polyMPO{

/** \brief calculate the basis and gradient of Basis for a give MP with its element Vtxs
 *
 *  \details based on the 4.1 section from:
 *  https://www.mn.uio.no/math/english/people/aca/michaelf/papers/gbc.pdf
 *
 *  \param MP: single Material Point
 *
 *  \param numVtxs: number of vertices of vtxCoords
 *
 *  \param vtxCoords: vertices of its corresponding element
 *
 *  \param basis && gadBasis: hold the return values
 *
 *  \return basis and gradient of basis
 * */
//TODO:Change this to support 3d
KOKKOS_INLINE_FUNCTION
void getBasisAndGradByAreaGblForm2d(Vec2d MP,
                                  int numVtxs,
                                  Vec2d* vtxCoords,
                                  double* basis,
                                  Vec2d* gradBasis){
    Vec2d e[maxVtxsPerElm + 1];
    Vec2d p[maxVtxsPerElm];
    double w[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        e[i + 1] = vtxCoords[i + 1] - vtxCoords[i];
        p[i] = vtxCoords[i] - MP;
    }
    e[0] = e[numVtxs];

    double c[maxVtxsPerElm];
    double a[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        c[i] = e[i].cross(e[i + 1]);
        a[i] = p[i].cross(e[i + 1]);
    }
    double wSum = 0.0;

    double wdx[maxVtxsPerElm];
    double wdy[maxVtxsPerElm];
    initArray(wdx, maxVtxsPerElm, 0.0);
    initArray(wdy, maxVtxsPerElm, 0.0);
    double wdxSum = 0.0;
    double wdySum = 0.0;
    for (int i = 0; i < numVtxs; i++){
        double aProduct = 1.0;
        for (int j = 0; j < numVtxs - 2; j++){
            int index1 = (j + i + 1) % numVtxs;
            aProduct *= a[index1];

            double productX = 1.0;
            double productY = 1.0;
            for (int k = 0; k < j; k++){
                int index2 = (i + k + 1) % numVtxs;
                productX *= a[index2];
                productY *= a[index2];
            }
            productX *= -(e[index1+1][1]);
            productY *= e[index1+1][0];
            for (int k = j + 1; k < numVtxs - 2; k++){
                int index2 = (i + k + 1) % numVtxs;
                productX *= a[index2];
                productY *= a[index2];
            }
            wdx[i] += productX;
            wdy[i] += productY;
        }
        wdx[i] *= c[i];
        wdy[i] *= c[i];
        wdxSum += wdx[i];
        wdySum += wdy[i];
        w[i] = c[i] * aProduct;
        wSum += w[i];
    }

    double wSumInv = 1.0 / wSum;
    for (int i = 0; i < numVtxs; i++){
        basis[i] = w[i] * wSumInv;
        gradBasis[i] = Vec2d(wdx[i] * wSumInv - w[i] * wSumInv * wSumInv * wdxSum, wdy[i] * wSumInv - w[i] * wSumInv * wSumInv * wdySum);
    }
}

KOKKOS_INLINE_FUNCTION
void getBasisAndGradByAreaGblForm3d(Vec3d MP,
                                  int numVtxs,
                                  Vec3d* vtxCoords,
                                  double* basis,
                                  Vec3d* gradBasis){
    Vec3d e[maxVtxsPerElm + 1];
    Vec3d p[maxVtxsPerElm];
    double w[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        e[i + 1] = vtxCoords[i + 1] - vtxCoords[i];
        p[i] = vtxCoords[i] - MP;
    }
    e[0] = e[numVtxs];

    double c[maxVtxsPerElm];
    double a[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        c[i] = e[i].cross(e[i + 1]).magnitude();
        a[i] = p[i].cross(e[i + 1]).magnitude();
    }
    double wSum = 0.0;

    double wdx[maxVtxsPerElm];
    double wdy[maxVtxsPerElm];
    double wdz[maxVtxsPerElm];
    initArray(wdx, maxVtxsPerElm, 0.0);
    initArray(wdy, maxVtxsPerElm, 0.0);
    initArray(wdz, maxVtxsPerElm, 0.0);
    double wdxSum = 0.0;
    double wdySum = 0.0;
    double wdzSum = 0.0;
    for (int i = 0; i < numVtxs; i++){
        double aProduct = 1.0;
        for (int j = 0; j < numVtxs - 2; j++){
            int index1 = (j + i + 1) % numVtxs;
            aProduct *= a[index1];

            double productX = 1.0;
            double productY = 1.0;
            double productZ = 1.0;
            for (int k = 0; k < j; k++){
                int index2 = (i + k + 1) % numVtxs;
                productX *= a[index2];
                productY *= a[index2];
                productZ *= a[index2];
            }

            // when j = k, find gradient of A_k
            double c1 = e[index1+1][0];
            double c2 = e[index1+1][1];
            double c3 = e[index1+1][2];
            double f1 = c2 * (-p[index1][2]) - c3 * (-p[index1][1]);
            double f2 = c3 * (-p[index1][0]) - c1 * (-p[index1][2]);
            double f3 = c1 * (-p[index1][1]) - c2 * (-p[index1][0]);
            productX *= (f2 * c3 - f3 * c2)/a[index1];
            productY *= (f3 * c1 - f1 * c3)/a[index1];
            productZ *= (f1 * c2 - f2 * c1)/a[index1];

            for (int k = j + 1; k < numVtxs - 2; k++){
                int index2 = (i + k + 1) % numVtxs;
                productX *= a[index2];
                productY *= a[index2];
                productZ *= a[index2];
            }
            wdx[i] += productX;
            wdy[i] += productY;
            wdz[i] += productZ;
        }
        wdx[i] *= c[i];
        wdy[i] *= c[i];
        wdz[i] *= c[i];
        wdxSum += wdx[i];
        wdySum += wdy[i];
        wdzSum += wdz[i];
        w[i] = c[i] * aProduct;
        wSum += w[i];
    }

    double wSumInv = 1.0 / wSum;
    for (int i = 0; i < numVtxs; i++){
        basis[i] = w[i] * wSumInv;
        gradBasis[i] = Vec3d(wdx[i] * wSumInv - w[i] * wSumInv * wSumInv * wdxSum,
                             wdy[i] * wSumInv - w[i] * wSumInv * wSumInv * wdySum,
                             wdz[i] * wSumInv - w[i] * wSumInv * wSumInv * wdzSum);
    }
}

KOKKOS_INLINE_FUNCTION
void calcBasis(int numVtxs, double* a, double* c, double* basis){
    double w[maxVtxsPerElm];
    double wSum = 0.0;
    for (int i = 0; i < numVtxs; i++){
        double aProduct = 1.0;
        for (int j = 0; j < numVtxs - 2; j++){
            int index1 = (j + i + 1) % numVtxs;
            aProduct *= a[index1];
        }
        w[i] = c[i] * aProduct;
        wSum += w[i];
    }

    double wSumInv = 1.0 / wSum;
    for (int i = 0; i < numVtxs; i++){
        basis[i] = w[i] * wSumInv;
    }
}

KOKKOS_INLINE_FUNCTION
void getBasisByAreaGblForm(Vec2d MP, int numVtxs, Vec2d* vtxCoords, double* basis) {
    Vec2d e[maxVtxsPerElm + 1];
    Vec2d p[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        e[i + 1] = vtxCoords[i + 1] - vtxCoords[i];
        p[i] = vtxCoords[i] - MP;
    }
    e[0] = e[numVtxs];

    double c[maxVtxsPerElm];
    double a[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        c[i] = e[i].cross(e[i + 1]);
        a[i] = p[i].cross(e[i + 1]);
    }

    calcBasis(numVtxs, a, c, basis);
}

KOKKOS_INLINE_FUNCTION
void getBasisByAreaGblForm3d(Vec3d MP, int numVtxs, Vec3d* vtxCoords, double* basis) {
    Vec3d e[maxVtxsPerElm + 1];
    Vec3d p[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        e[i + 1] = vtxCoords[i + 1] - vtxCoords[i];
        p[i] = vtxCoords[i] - MP;
    }
    e[0] = e[numVtxs];

    double c[maxVtxsPerElm];
    double a[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        c[i] = (e[i].cross(e[i + 1])).magnitude();
        a[i] = (p[i].cross(e[i + 1])).magnitude();
    }

    calcBasis(numVtxs, a, c, basis);
}

//3d
KOKKOS_INLINE_FUNCTION
void getBasisByAreaGblFormSpherical(Vec3d MP, int numVtxs, Vec3d* v,
                                    double radius, double* basis) {
    double c[maxVtxsPerElm];
    double a[maxVtxsPerElm];
    for (int i = 1; i < numVtxs; i++){
        //c = vi-1 vi vi+1
        //a = x    vi vi+1
        c[i] = sphericalTriangleArea(v[i-1],v[i],v[i+1],radius);
        a[i] = sphericalTriangleArea(v[i],v[i+1],MP,radius);
    }
    c[0] = sphericalTriangleArea(v[numVtxs-1],v[0],v[1],radius);
    a[0] = sphericalTriangleArea(v[0],v[1],MP,radius);

    calcBasis(numVtxs, a, c, basis);
}

KOKKOS_INLINE_FUNCTION
void getBasisByAreaGblFormSpherical2(Vec3d MP, int numVtxs, Vec3d* v,
                                    double radius, double* basis) {
    double c[maxVtxsPerElm];
    double a[maxVtxsPerElm];
    for (int i = 1; i < numVtxs; i++){
        //c = vi-1 vi vi+1
        //a = x    vi vi+1
        c[i] = sphericalTriangleArea2(v[i-1],v[i],v[i+1],radius);
        a[i] = sphericalTriangleArea2(v[i],v[i+1],MP,radius);
    }
    c[0] = sphericalTriangleArea2(v[numVtxs-1],v[0],v[1],radius);
    a[0] = sphericalTriangleArea2(v[0],v[1],MP,radius);
    
    calcBasis(numVtxs, a, c, basis);
}

/*
KOKKOS_INLINE_FUNCTION
void getBasisByAreaGblForm_1(Vec2d MP, int numVtxs, Vec2d* vtxCoords, double* basis) {
    double denominator, product;
    Vec2d v1, v2;

    denominator = 0.0;
    for (int i = 1; i <= numVtxs; i++) {
        v1 = vtxCoords[i] - vtxCoords[i-1];
        v2 = vtxCoords[i+1] - vtxCoords[i];
        product = 0.5 * v1.cross(v2);
        for (int k = 1; k <= i-2; k++) {
            v1 = vtxCoords[k] - MP;
            v2 = vtxCoords[k+1] - MP;
            product *= 0.5 * v1.cross(v2);
        }
        for (int k = i+1; k <= std::min(i-2,1) + numVtxs; k++) {
            v1 = vtxCoords[k] - MP;
            v2 = vtxCoords[k+1] - MP;
            product *= 0.5 * v1.cross(v2);
        }
        basis[i-1] = product;
        denominator += product;
    }
    for (int i = 0; i < numVtxs; i++) {
        basis[i] /= denominator;
    }
} 
*/

// spherical interpolation of values from mesh vertices to MPsi
template <MeshFieldIndex meshFieldIndex>
void sphericalInterpolation(MPMesh& mpMesh){
    Kokkos::Timer timer;
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    int numVtxs = p_mesh->getNumVertices();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();
   
    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();
    double radius = p_mesh->getSphereRadius();
    PMT_ALWAYS_ASSERT(radius >0);
    constexpr MaterialPointSlice mpfIndex = meshFieldIndexToMPSlice<meshFieldIndex>;
    auto mpField = p_MPs->getData<mpfIndex>();
    
    const int numEntries = mpSliceToNumEntries<mpfIndex>();
    //check field correspondence
    auto meshField = p_mesh->getMeshField<meshFieldIndex>(); 

    auto interpolation = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if(mask) { //if material point is 'active'/'enabled'
            Vec3d position3d(MPsPosition(mp,0),MPsPosition(mp,1),MPsPosition(mp,2));
            // formating 
            Vec3d v3d[maxVtxsPerElm+1];
            int numVtx = elm2VtxConn(elm,0);
            for(int i = 1; i<=numVtx; i++){
                v3d[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);
                v3d[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
                v3d[i-1][2] = vtxCoords(elm2VtxConn(elm,i)-1,2);
            }
            v3d[numVtx][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
            v3d[numVtx][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
            v3d[numVtx][2] = vtxCoords(elm2VtxConn(elm,1)-1,2);
            
            double basisByArea3d[maxVtxsPerElm] = {0.0};
            initArray(basisByArea3d,maxVtxsPerElm,0.0);

            // calc basis
            getBasisByAreaGblFormSpherical(position3d, numVtx, v3d, radius, basisByArea3d);
            
            // interpolation step
            for(int entry=0; entry<numEntries; entry++){
                double mpValue = 0.0;
                for(int i=1; i<= numVtx; i++){
                    mpValue += meshField(elm2VtxConn(elm,i)-1,entry)*basisByArea3d[i-1];
                }
                mpField(mp,entry) = mpValue;
            }
        }
    };
    p_MPs->parallel_for(interpolation, "interpolation");
    pumipic::RecordTime("PolyMPO_sphericalInterpolation", timer.seconds());
}


KOKKOS_INLINE_FUNCTION
void compute2DplanarTriangleArea(int numVtx, 
     const Kokkos::View<double[maxVtxsPerElm][2], Kokkos::LayoutStride, Kokkos::MemoryTraits<Kokkos::Unmanaged>>& gnom_vtx_subview, 
     double mpProjX, double mpProjY, double* basis){

  // Temporary storage
  double vertCoords[2][maxVtxsPerElm + 1];
  for (int i = 0; i < numVtx; ++i) {
    vertCoords[0][i] = gnom_vtx_subview(i, 0);
    vertCoords[1][i] = gnom_vtx_subview(i, 1);
  }
  vertCoords[0][numVtx] = vertCoords[0][0];
  vertCoords[1][numVtx] = vertCoords[1][0];
  
  //Helper lambda for signed triangle area
  auto triArea = [&](const double p1[2], const double p2[2], const double p3[2]) -> double {
    return 0.5 * (p1[0] * (p2[1] - p3[1]) - p2[0] * (p1[1] - p3[1]) + p3[0] * (p1[1] - p2[1]));
  };
  
  // Compute areaV and areaXV
  double areaV[maxVtxsPerElm];
  double areaXV[maxVtxsPerElm];
  double xy[2] = { mpProjX, mpProjY };
  
  //Special case
  double p1[2] = { vertCoords[0][numVtx - 1], vertCoords[1][numVtx - 1] };
  double p2[2] = { vertCoords[0][0], vertCoords[1][0] };
  double p3[2] = { vertCoords[0][1], vertCoords[1][1] };
  areaV[0] = triArea(p1, p2, p3);
  double q1[2] = { vertCoords[0][0], vertCoords[1][0] };
  double q2[2] = { xy[0], xy[1] };
  double q3[2] = { vertCoords[0][1], vertCoords[1][1] };
  areaXV[0] = triArea(q1, q2, q3);
  
  for (int i = 1; i < numVtx; ++i) {
    double p1[2] = { vertCoords[0][i - 1], vertCoords[1][i - 1] };
    double p2[2] = { vertCoords[0][i], vertCoords[1][i] };
    double p3[2] = { vertCoords[0][i + 1], vertCoords[1][i + 1] };
    areaV[i] = triArea(p1, p2, p3);
    double q1[2] = { vertCoords[0][i], vertCoords[1][i] };
    double q2[2] = { xy[0], xy[1] };
    double q3[2] = { vertCoords[0][i + 1], vertCoords[1][i + 1] };
    areaXV[i] = triArea(q1, q2, q3);
  }
  
  // Compute Wachspress-like weights
  double denominator = 0.0;
  for (int i = 0; i < numVtx; ++i){
    double product = areaV[i];
    for (int j = 0; j < numVtx - 2; ++j) {
      int ind1 = (i + j + 1) % numVtx;
      product *= areaXV[ind1];
    }
    basis[i] = product;
    denominator += product;
  }

  // Normalize
  for (int i = 0; i < numVtx; ++i){
    basis[i] /= denominator;
    //printf("i %d basis %.15e \n", i, basis[i]);
  }
}

inline void sphericalInterpolationDispVelIncr(MPMesh& mpMesh){
    Kokkos::Timer timer;
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    int numVtxs = p_mesh->getNumVertices();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();
    
    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();
    double radius = p_mesh->getSphereRadius();
    PMT_ALWAYS_ASSERT(radius > 0);
 
    constexpr MeshFieldIndex meshFieldIndex1 = polyMPO::MeshF_RotLatLonIncr;
    constexpr MeshFieldIndex meshFieldIndex2 = polyMPO::MeshF_OnSurfVeloIncr;
    
    auto meshField1 = p_mesh->getMeshField<meshFieldIndex1>();
    auto meshField2 = p_mesh->getMeshField<meshFieldIndex2>();
    
    constexpr MaterialPointSlice mpfIndex1 = meshFieldIndexToMPSlice<meshFieldIndex1>;
    constexpr MaterialPointSlice mpfIndex2 = meshFieldIndexToMPSlice<meshFieldIndex2>;

    const int numEntries1 = mpSliceToNumEntries<mpfIndex1>();
    const int numEntries2 = mpSliceToNumEntries<mpfIndex2>();

    auto mpField1 = p_MPs->getData<mpfIndex1>();
    auto mpField2 = p_MPs->getData<mpfIndex2>();

    // Field required for calculting gnomonic projection of MPs
    bool use3DArea = false;
    bool isRotated = p_mesh->getRotatedFlag();
    auto gnomProjVtx = p_mesh->getMeshField<polyMPO::MeshF_VtxGnomProj>();
    auto gnomProjElmCenter = p_mesh->getMeshField<polyMPO::MeshF_ElmCenterGnomProj>();
      
    auto interpolation = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
      if(mask) {

        double basisByArea[maxVtxsPerElm] = {0.0};
        initArray(basisByArea, maxVtxsPerElm, 0.0);
        int numVtx = elm2VtxConn(elm, 0);
        Vec3d position3d(MPsPosition(mp, 0), MPsPosition(mp, 1), MPsPosition(mp, 2));

        if(use3DArea){
          Vec3d v3d[maxVtxsPerElm + 1];
          for (int i = 1; i <= numVtx; i++) {
            v3d[i-1][0] = vtxCoords(elm2VtxConn(elm, i) - 1, 0);
            v3d[i-1][1] = vtxCoords(elm2VtxConn(elm, i) - 1, 1);
            v3d[i-1][2] = vtxCoords(elm2VtxConn(elm, i) - 1, 2);
          }
          v3d[numVtx][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
          v3d[numVtx][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
          v3d[numVtx][2] = vtxCoords(elm2VtxConn(elm,1)-1,2);        
          getBasisByAreaGblFormSpherical(position3d, numVtx, v3d, radius, basisByArea);
        }
        else{ //if using gnomonic Projection for weights
          double mpProjX, mpProjY;
          auto gnomProjElmCenter_sub = Kokkos::subview(gnomProjElmCenter, elm, Kokkos::ALL);
          if(isRotated){
            position3d[0] = -MPsPosition(mp, 2);
            position3d[2] = MPsPosition(mp, 0);
          }
          computeGnomonicProjectionAtPoint(position3d, gnomProjElmCenter_sub, mpProjX, mpProjY);
          auto gnom_vtx_subview = Kokkos::subview(gnomProjVtx, elm, Kokkos::ALL, Kokkos::ALL); 
          compute2DplanarTriangleArea(numVtx, gnom_vtx_subview, mpProjX, mpProjY, basisByArea);
        }
        
        for(int entry=0; entry<numEntries1; entry++){
          double mpValue = 0.0;
          for(int i=1; i<= numVtx; i++)
            mpValue += meshField1(elm2VtxConn(elm,i)-1,entry)*basisByArea[i-1];
          mpField1(mp,entry) = mpValue;
        }
        
        for(int entry=0; entry<numEntries2; entry++){
          double mpValue = 0.0;
          for(int i=1; i<= numVtx; i++)
            mpValue += meshField2(elm2VtxConn(elm,i)-1,entry)*basisByArea[i-1];
          mpField2(mp,entry) = mpValue;
        } 
      }//mask
    };//lambda
    p_MPs->parallel_for(interpolation, "sphericalInterpolationMultiField");
    pumipic::RecordTime("PolyMPO_sphericalInterpolationDispVelIncr", timer.seconds());
}

} //namespace polyMPO end
#endif
