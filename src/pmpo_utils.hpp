#ifndef POLYMPO_UTILS_H
#define POLYMPO_UTILS_H

#include <stdlib.h>
#include <iostream>

#include <Kokkos_Core.hpp>

//assert borrowed from scorec/core/pcu/pcu_util.h
void PMT_Assert_Fail(const char* msg) __attribute__ ((noreturn));

#define PMT_ALWAYS_ASSERT(cond)                   \
  do {                                            \
    if (! (cond)) {                               \
      char omsg[2048];                            \
      sprintf(omsg, "%s failed at %s + %d \n",    \
        #cond, __FILE__, __LINE__);               \
      PMT_Assert_Fail(omsg);                      \
    }                                             \
  } while (0)

namespace polyMPO{

#define vec2d_nEntries 2
typedef double vec2d_t[vec2d_nEntries];
#define vec3d_nEntries 3
typedef double vec3d_t[vec3d_nEntries];

class Vec2d;
class Vec3d;
using Vec2dView = Kokkos::View<Vec2d*>;
using DoubleView = Kokkos::View<double*>;
using IntView = Kokkos::View<int*>;
using BoolView = Kokkos::View<bool*>;

class Vec2d {
  private:
    vec2d_t coords_;

  public:
    //constructors
    KOKKOS_INLINE_FUNCTION
    Vec2d():coords_{0.0, 0.0}{}
    
    ~Vec2d() = default;
    
    KOKKOS_INLINE_FUNCTION
    Vec2d(double x1, double x2):coords_{x1,x2}{}
    
    KOKKOS_INLINE_FUNCTION
    Vec2d(double x[2]):coords_{x[0], x[1]}{}

    //operators 
    KOKKOS_INLINE_FUNCTION
    double operator[](int i) const { return coords_[i]; }
    
    KOKKOS_INLINE_FUNCTION
    double &operator[](int i) { return coords_[i]; }

    KOKKOS_INLINE_FUNCTION
    Vec2d operator-() { return Vec2d(-coords_[0], -coords_[1]); }
    
    KOKKOS_INLINE_FUNCTION
    Vec2d operator+(Vec2d v) { return Vec2d(coords_[0] + v[0],
                                                 coords_[1] + v[1]); }
    
    KOKKOS_INLINE_FUNCTION
    Vec2d operator-(Vec2d v) { return Vec2d(coords_[0] - v[0],
                                                 coords_[1] - v[1]); }
    
    KOKKOS_INLINE_FUNCTION
    Vec2d operator*(double scalar) const { return Vec2d(coords_[0]*scalar,
                                                            coords_[1]*scalar); }
    
    KOKKOS_INLINE_FUNCTION
    double dot(Vec2d v) { return  coords_[0]*v[0]+coords_[1]*v[1]; }
    
    KOKKOS_INLINE_FUNCTION
    double cross(Vec2d v) { return (coords_[0]*v[1] - coords_[1]*v[0]); }
    
    KOKKOS_INLINE_FUNCTION
    double magnitude() {return sqrt(coords_[0]*coords_[0] +
                                    coords_[1]*coords_[1]); }
};

class Vec3d{
  private:
    vec3d_t coords_;

  public:
    //constructors
    KOKKOS_INLINE_FUNCTION
    Vec3d():coords_{0.0, 0.0, 0.0}{}
    
    ~Vec3d() = default;
    
    KOKKOS_INLINE_FUNCTION
    Vec3d(double x1, double x2, double x3):coords_{x1,x2,x3}{}
    
    KOKKOS_INLINE_FUNCTION
    Vec3d(double x[3]):coords_{x[0], x[1], x[2]}{}

    //operators
    KOKKOS_INLINE_FUNCTION
    double operator[](int i) const { return coords_[i]; }
    
    KOKKOS_INLINE_FUNCTION
    double &operator[](int i) { return coords_[i]; }
    
    KOKKOS_INLINE_FUNCTION
    Vec3d operator+(const Vec3d& v) const {
        return Vec3d(coords_[0] + v.coords_[0], coords_[1] + v.coords_[1], coords_[2] + v.coords_[2]);
    }

    KOKKOS_INLINE_FUNCTION
    Vec3d operator-(const Vec3d& v) const {
        return Vec3d(coords_[0] - v.coords_[0], coords_[1] - v.coords_[1], coords_[2] - v.coords_[2]);
    }
    
    KOKKOS_INLINE_FUNCTION
    Vec3d operator-() { return Vec3d(-coords_[0], -coords_[1], -coords_[2]); }

    KOKKOS_INLINE_FUNCTION
    Vec3d operator*(double scalar) const {
        return Vec3d(coords_[0] * scalar, coords_[1] * scalar, coords_[2] * scalar);
    }

    KOKKOS_INLINE_FUNCTION
    double dot(const Vec3d& v) const {
        return coords_[0] * v.coords_[0] + coords_[1] * v.coords_[1] + coords_[2] * v.coords_[2];
    }

    KOKKOS_INLINE_FUNCTION
    Vec3d cross(const Vec3d& v) const {
        return Vec3d(coords_[1] * v.coords_[2] - coords_[2] * v.coords_[1],
                     coords_[2] * v.coords_[0] - coords_[0] * v.coords_[2],
                     coords_[0] * v.coords_[1] - coords_[1] * v.coords_[0]);
    }

    KOKKOS_INLINE_FUNCTION
    double magnitude() const {
        return std::sqrt(coords_[0] * coords_[0] + coords_[1] * coords_[1] + coords_[2] * coords_[2]);
    }
};

KOKKOS_INLINE_FUNCTION
void initArray(Vec2d* arr, int n, Vec2d fill){
    for(int i=0; i<n; i++){
        arr[i] = fill;
    }
}

KOKKOS_INLINE_FUNCTION
void initArray(Vec3d* arr, int n, Vec3d fill){
    for(int i=0; i<n; i++){
        arr[i] = fill;
    }
}

KOKKOS_INLINE_FUNCTION
void initArray(double* arr, int n, double fill){
    for(int i=0; i<n; i++){
        arr[i] = fill;
    }
}


KOKKOS_INLINE_FUNCTION
double arcLength(Vec3d &a, Vec3d &b){
    Vec3d c = b-a;
    double r = a.magnitude();
    return r * 2.0 * std::asin(c.magnitude() / (2.0 * r)) ;
}

// implement from: https://github.com/MPAS-Seaice-MPM-Project/MPAS-Seaice-MPM/blob/9bb3b4f3a09ac59fbc93b8ab317ca0fa6a6f18bb/components/mpas-framework/src/operators/mpas_geometry_utils.F#L651-L671
//TODO:check for more effcient calculation
KOKKOS_INLINE_FUNCTION
double sphericalTriangleArea(Vec3d &a, Vec3d &b, Vec3d &c, double radius){
    double ab, bc, ca, semiperim, tanqe;
    Vec3d ablen, aclen, dlen;

    ab = arcLength(a,b) / radius;
    bc = arcLength(b,c) / radius;
    ca = arcLength(c,a) / radius;
    semiperim = 0.5 * (ab + bc + ca);

    tanqe = Kokkos::sqrt(Kokkos::fmax(0.0, Kokkos::tan(0.5 * semiperim) * 
                                           Kokkos::tan(0.5 * (semiperim - ab)) *
                                           Kokkos::tan(0.5 * (semiperim - bc)) * 
                                           Kokkos::tan(0.5 * (semiperim - ca))));

    return 4.0 * radius * radius * Kokkos::atan(tanqe);
}

//implement from: https://www.maa.org/sites/default/files/Eriksson14108673.pdf 
KOKKOS_INLINE_FUNCTION
double sphericalTriangleArea2(Vec3d &a, Vec3d &b, Vec3d &c, double radius){
    double inv_radius = 1.0/radius;
    Vec3d a_unit = a*inv_radius;
    Vec3d b_unit = b*inv_radius;
    Vec3d c_unit = c*inv_radius;
    double tripleProduct = a_unit.dot(b_unit.cross(c_unit));
    double tangent = tripleProduct / (1 + b_unit.dot(c_unit) + c_unit.dot(a_unit) + a_unit.dot(b_unit));

    return Kokkos::fabs(2.0 * radius * radius * Kokkos::atan(tangent));
}

//this is a lazy comparison and shouldn't be relied on beyond simple testing
bool isEqual(double a, double b, double tol=1e-9);

}//namespace polyMPO end

#endif

