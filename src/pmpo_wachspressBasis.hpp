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
void getBasisAndGradByAreaGblFormSpherical(Vec3d MP,
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
    double v[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        c[i] = e[i].cross(e[i + 1]).magnitude();
        a[i] = p[i].cross(e[i + 1]).magnitude();
        v[i] = vtxCoords[i].dot(MP);
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
        w[i] = c[i] * aProduct;
        wdxSum += wdx[i] * v[i] + w[i] * vtxCoords[i][0];
        wdySum += wdy[i] * v[i] + w[i] * vtxCoords[i][1];
        wdzSum += wdz[i] * v[i] + w[i] * vtxCoords[i][2];
        wSum += w[i] * v[i];
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
void calcSphericalBasis(int numVtxs, double* a, double* c, double* vtxDotMP, double* basis){
    double w[maxVtxsPerElm];
    double wSum = 0.0;
    for (int i = 0; i < numVtxs; i++){
        double aProduct = 1.0;
        for (int j = 0; j < numVtxs - 2; j++){
            int index1 = (j + i + 1) % numVtxs;
            aProduct *= a[index1];
        }
        w[i] = c[i] * aProduct;
        wSum += w[i] * vtxDotMP[i];
    }

    double wSumInv = 1.0 / wSum;
    for (int i = 0; i < numVtxs; i++){
        basis[i] = w[i] * wSumInv;
    }
}


KOKKOS_INLINE_FUNCTION
void getBasisByAreaGblForm2d(Vec2d MP, int numVtxs, Vec2d* vtxCoords, double* basis) {
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


KOKKOS_INLINE_FUNCTION
void getBasisByAreaGblFormSpherical3(Vec3d MP, int numVtxs, Vec3d* vtxCoords, double* basis) {
    Vec3d e[maxVtxsPerElm + 1];
    Vec3d p[maxVtxsPerElm];
    for (int i = 0; i < numVtxs; i++){
        e[i + 1] = vtxCoords[i + 1] - vtxCoords[i];
        p[i] = vtxCoords[i] - MP;
    }
    e[0] = e[numVtxs];

    double c[maxVtxsPerElm];
    double a[maxVtxsPerElm];
    double vDotMP[maxVtxsPerElm]; // v_j * x
    for (int i = 0; i < numVtxs; i++){
        c[i] = (e[i].cross(e[i + 1])).magnitude();
        a[i] = (p[i].cross(e[i + 1])).magnitude();
        vDotMP[i] = vtxCoords[i].dot(MP);
    }

    calcSphericalBasis(numVtxs, a, c, vDotMP, basis);
}

// spherical interpolation of values from mesh vertices to MPsi
template <MeshFieldIndex mfIndex, MaterialPointSlice mpfIndex>
void sphericalInterpolation(MPMesh& mpMesh){
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    int numVtxs = p_mesh->getNumVertices();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();
   
    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();
    double radius = p_mesh->getSphereRadius();
    PMT_ALWAYS_ASSERT(radius >0);
    auto mpField = p_MPs->getData<mpfIndex>();
    
    const int numEntries = mpSlice2MeshFieldIndex.at(mpfIndex).first;
    //check field correspondence
    const MeshFieldIndex meshFieldIndex = mpSlice2MeshFieldIndex.at(mpfIndex).second;
    PMT_ALWAYS_ASSERT(meshFieldIndex == mfIndex);
    auto meshField = p_mesh->getMeshField<mfIndex>(); 

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
            getBasisByAreaGblFormSpherical2(position3d, numVtx, v3d, radius, basisByArea3d);
            
            // interpolation step
            for(int entry=0; entry<numEntries; entry++){
                double mpValue = 0.0;
                for(int i=0; i<= numVtx; i++){
                    mpValue += meshField(elm2VtxConn(elm,i),entry)*basisByArea3d[i];
                }
                mpField(mp,entry) = mpValue;
            }
        }
    };
    p_MPs->parallel_for(interpolation, "interpolation");
}

} //namespace polyMPO end
#endif
