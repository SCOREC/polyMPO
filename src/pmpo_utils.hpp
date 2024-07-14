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
typedef double doubleSclr_t[1];
#define vec4d_nEntries 4
typedef double vec4d_t[vec4d_nEntries];


class Vec2d;
class Vec3d;
class Vec4d;
class Matrix;
using Vec2dView = Kokkos::View<Vec2d*>;
using Vec3dView = Kokkos::View<Vec3d*>;
using Vec4dView = Kokkos::View<Vec4d*>;
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

class Vec4d{
  private:
    vec4d_t coords_;

  public:
    //constructors
    KOKKOS_INLINE_FUNCTION
    Vec4d():coords_{0.0, 0.0, 0.0, 0.0}{}
    
    ~Vec4d() = default;
    
    KOKKOS_INLINE_FUNCTION
    Vec4d(double x1, double x2, double x3, double x4):coords_{x1,x2,x3,x4}{}
    
    KOKKOS_INLINE_FUNCTION
    Vec4d(double x[4]):coords_{x[0], x[1], x[2], x[3]}{}

    //operators
    KOKKOS_INLINE_FUNCTION
    double operator[](int i) const { return coords_[i]; }
    
    KOKKOS_INLINE_FUNCTION
    double &operator[](int i) { return coords_[i]; }
    
    KOKKOS_INLINE_FUNCTION
    Vec4d operator+(const Vec4d& v) const {
        return Vec4d(coords_[0] + v.coords_[0], coords_[1] + v.coords_[1], 
                     coords_[2] + v.coords_[2], coords_[3] + v.coords_[3]);
    }

    KOKKOS_INLINE_FUNCTION
    Vec4d operator-(const Vec4d& v) const {
        return Vec4d(coords_[0] - v.coords_[0], coords_[1] - v.coords_[1], 
                     coords_[2] - v.coords_[2], coords_[3] - v.coords_[3]);
    }
    
    KOKKOS_INLINE_FUNCTION
    Vec4d operator-() { return Vec4d(-coords_[0], -coords_[1], -coords_[2], -coords_[3]); }

    KOKKOS_INLINE_FUNCTION
    Vec4d operator*(double scalar) const {
        return Vec4d(coords_[0] * scalar, coords_[1] * scalar, coords_[2] * scalar, coords_[3] * scalar);
    }

    KOKKOS_INLINE_FUNCTION
    double dot(const Vec4d& v) const {
        return coords_[0] * v.coords_[0] + coords_[1] * v.coords_[1] + coords_[2] * v.coords_[2] + coords_[3] * v.coords_[3];
    }

    KOKKOS_INLINE_FUNCTION
    double magnitude() const {
        return std::sqrt(coords_[0] * coords_[0] + coords_[1] * coords_[1] + coords_[2] * coords_[2] + coords_[3] * coords_[3]);
    }
};

class Matrix {
  //2x2, 3x3, or 4x4
  //stored as such: 
  /*
    | v0[0] v0[1] v0[2] v0[3] |
    | v1[0] v1[1] v1[2] v1[3] |
    | v2[0] v2[1] v2[2] v2[3] |
    | v3[0] v3[1] v3[2] v3[3] |
  */
  private:
    std::vector<std::vector<double>> data_;
    int rows_;
    int cols_;

  public:
  
  KOKKOS_INLINE_FUNCTION
  Matrix(Vec4d v0, Vec4d v1, Vec4d v2, Vec4d v3){
    data_ = std::vector<std::vector<double>>(4, std::vector<double>(4, 0));
    rows_ = 4;
    cols_ = 4;
    
    for (int i=0; i<4; i++){
      data_[0][i] = v0[i];
      data_[1][i] = v1[i];
      data_[2][i] = v2[i];
      data_[3][i] = v3[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  std::vector<double> solve_helper(std::vector<double> b){
    std::vector<double> x(rows_, 0);
    for(int i=0; i<rows_; i++){
      x[i] = 0;
    }
    for (int i=0; i<rows_; i++){
      double pivot = data_[i][i];
      for( int j = i+1; j<rows_; j++){
        double ratio = data_[j][i] / pivot;
        for(int k = i; k<cols_; k++){
          data_[j][k] -= ratio * data_[i][k];
        }
        b[j] -= ratio * b[i];
      }
    }
    for(int i = rows_-1; i>=0; i--){
      double sum = 0;
      for(int j = i+1; j<cols_; j++){
        sum += data_[i][j] * x[j];
      }
      x[i] = (b[i] - sum) / data_[i][i];
    }
    return x;
  }

  Vec4d solve(Vec4d b){
    std::vector<double> x = solve_helper(b);
    return Vec4d(x[0], x[1], x[2], x[3]);
  }
  
  KOKKOS_INLINE_FUNCTION
  double operator()(int i, int j) const { return data_[i][j];}
  

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
void initArray(Vec4d* arr, int n, Vec4d fill){
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

    return 4.0 * Kokkos::atan(tanqe);
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

    return 2.0 * Kokkos::atan(tangent);
}

//this is a lazy comparison and shouldn't be relied on beyond simple testing
bool isEqual(double a, double b, double tol=1e-9);

}//namespace polyMPO end

#endif

