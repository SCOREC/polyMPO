#ifndef POLYMPMTEST_UTILS_H
#define POLYMPMTEST_UTILS_H

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

namespace polyMpmTest{

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
    
    KOKKOS_FUNCTION
    void operator=(const double &scalar);
    //Vec2d &operator+=(const Vec2d &v);
   
};

KOKKOS_INLINE_FUNCTION
void Vec2d::operator=(const double &scalar){
    coords_[0] = scalar;
    coords_[1] = scalar;
}

KOKKOS_INLINE_FUNCTION
void initArray(Vec2d* arr, int n, Vec2d fill){
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
};

//this is a lazy comparison and shouldn't be relied on beyond simple testing
bool isEqual(double a, double b, double tol=1e-9);

}//namespace polyMpmTest end

#endif

