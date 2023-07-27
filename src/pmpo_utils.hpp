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

//TODO: discuss the sign determination of the area and the center of the origin
//double check the inline function
//check for more effcient calculation
KOKKOS_INLINE_FUNCTION
double sphereTriangleArea(Vec3d &a, Vec3d &b, Vec3d &c, double radius){
    double ab, bc, ca, semiperim, tanqe;
    Vec3d ablen, aclen, dlen;
    //PMT_ALWAYS_ASSERT();

    ab = arcLength(a,b) / radius;
    bc = arcLength(b,c) / radius;
    ca = arcLength(c,a) / radius;
    semiperim = 0.5 * (ab + bc + ca);

    tanqe = std::sqrt(std::tan(0.5 * semiperim) * 
                      std::tan(0.5 * (semiperim - ab)) *
                      std::tan(0.5 * (semiperim - bc)) * 
                      std::tan(0.5 * (semiperim - ca)));
    tanqe = tanqe > 0.0 ? tanqe : 0.0;

    double triangleArea = 4.0 * radius * radius * std::atan(tanqe);

    ablen = b-a;
    aclen = c-a;

    dlen[0] =  (ablen[1] * aclen[2]) - (ablen[2] * aclen[1]);
    dlen[1] = -((ablen[0] * aclen[2]) - (ablen[2] * aclen[0]));
    dlen[2] =  (ablen[0] * aclen[1]) - (ablen[1] * aclen[0]);

    if ((dlen[0] * a[0] + dlen[1] * a[1] + dlen[2] * a[2]) < 0.0) {
        triangleArea = -triangleArea;
    }

    return triangleArea;
}

//this is a lazy comparison and shouldn't be relied on beyond simple testing
bool isEqual(double a, double b, double tol=1e-9);

}//namespace polyMPO end

#endif

