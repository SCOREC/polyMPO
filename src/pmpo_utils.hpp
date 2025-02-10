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

const double EPSILON = std::numeric_limits<double>::epsilon();

#define vec2d_nEntries 2
typedef double vec2d_t[vec2d_nEntries];
#define vec3d_nEntries 3
typedef double vec3d_t[vec3d_nEntries];
#define vec4d_nEntries 4
typedef double vec4d_t[vec4d_nEntries];

typedef double doubleSclr_t[1];

class Vec2d;
class Vec3d;
class Vec4d;
class Matrix4d;

using Vec2dView = Kokkos::View<Vec2d*>;
using Vec3dView = Kokkos::View<Vec3d*>;
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

class Matrix4d {
    
  private:
  
    double data_[4][4];
   
  public:
			        
    KOKKOS_INLINE_FUNCTION
    Matrix4d(Vec4d v0, Vec4d v1, Vec4d v2, Vec4d v3){
  
      for (int i=0; i<4; i++){
        data_[0][i] = v0[i];
        data_[1][i] = v1[i];
        data_[2][i] = v2[i];
        data_[3][i] = v3[i];
      }
    }
 
    //Retrieval
    KOKKOS_INLINE_FUNCTION
    double& operator()(int i, int j) { return data_[i][j];}

    //Scalar division operator overloading
    KOKKOS_INLINE_FUNCTION
    Matrix4d operator/(double scalar) const {
      Matrix4d result(Vec4d(0, 0, 0, 0), Vec4d(0, 0, 0, 0), Vec4d(0, 0, 0, 0), Vec4d(0, 0, 0, 0));
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          result(i, j) = data_[i][j] / scalar;
        }
      }
      return result;
    }

    // Function to calculate the Frobenius norm
    KOKKOS_INLINE_FUNCTION
    double frobeniusNorm() const {
      double sum = 0.0; // Initialize sum of squares
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          sum += data_[i][j]* data_[i][j]; // Sum of squares
        }
      }
      return std::sqrt(sum); // Return the square root of the sum
    }

    // Function to calculate the trace of the matrix
    KOKKOS_INLINE_FUNCTION
    double trace() const {
      double traceSum = 0.0; // Initialize sum of diagonal elements
      for (int i = 0; i < 4; i++) {
        traceSum += data_[i][i]; // Sum diagonal elements
      }
      return traceSum;
    }
    
    //Function to regularize the matrix
    KOKKOS_INLINE_FUNCTION
    void addToDiag(double eps) {
      for (int i = 1; i < 4; i++) {
        data_[i][i]=data_[i][i]+eps;
      }
    }

    KOKKOS_INLINE_FUNCTION
    void scaleFirstRowAndColumn(double factor) {
      // Scale the first row
      for (int j = 0; j < 4; j++) {
        data_[0][j] *= factor;
      }
      // Scale the first column
      for (int i = 0; i < 4; i++) {
        data_[i][0] *= factor;
      }
    }
				        
};


KOKKOS_INLINE_FUNCTION
void initArray(Vec2d* arr, int n, Vec2d fill){
    for(int i=0; i<n; i++){
        arr[i] = fill;
    }
}


KOKKOS_INLINE_FUNCTION
void CholeskySolve4d_UnitRHS(Matrix4d& A, double* x){

    double a_00=A(0,0);
    if (A(0,0)==0){
        x[0]=0;
	x[1]=0;
	x[2]=0;
	x[3]=0;
	return;
    }
    
    A(0,0) = std::sqrt(A(0,0));
    A(0,1) /= A(0,0);
    A(0,2) /= A(0,0);
    A(0,3) /= A(0,0);

    double diag = A(1,1) - A(0,1)*A(0,1);
    if(diag>EPSILON){
         A(1,1)=std::sqrt(diag);
	 A(1,2)=(A(1,2)-A(0,1)*A(0,2))/A(1,1);
	 A(1,3)=(A(1,3)-A(0,1)*A(0,3))/A(1,1);
    }
    else{
        A(1,1)=0.0;
	A(1,2)=0.0;
	A(1,3)=0.0;
    }
    
    diag = A(2,2) - A(0,2)*A(0,2)-A(1,2)*A(1,2);
    if(diag>EPSILON){
        A(2,2)=std::sqrt(diag);
        A(2,3)=(A(2,3)-A(0,2)*A(0,3)-A(1,2)*A(1,3))/A(2,2);
    }
    else{
       A(2,2)=0.0;
       A(2,3)=0.0;
    }

    diag=A(3,3)-A(0,3)*A(0,3)-A(1,3)*A(1,3)-A(2,3)*A(2,3);
    if(diag>EPSILON)
       A(3,3)=std::sqrt(diag);
    else 
       A(3,3)=0.0;

    if(abs(A(1,1))<EPSILON || abs(A(2,2))<EPSILON || abs(A(3,3))<EPSILON){
       x[0]=1.0/a_00;
       x[1]=0.0;
       x[2]=0.0;
       x[3]=0.0;
    }
    else{
       x[0]= 1.0/A(0,0);
       x[1]= -A(0,1)*x[0]/A(1,1);                
       x[2]= -(A(0,2)*x[0]+A(1,2)*x[1])/A(2,2);  
       x[3]= -(A(0,3)*x[0]+A(1,3)*x[1]+A(2,3)*x[2])/A(3,3); 

       x[3] = x[3]/A(3,3);
       x[2] = ( x[2] - A(2,3)*x[3] )/A(2,2);
       x[1] = ( x[1] - A(1,2)*x[2] - A(1,3)*x[3])/A(1,1);
       x[0] = ( x[0] - A(0,1)*x[1] - A(0,2)*x[2] - A(0,3)*x[3])/A(0,0);
    } 
}

KOKKOS_INLINE_FUNCTION
void QRDecomp4d(Matrix4d& A, Matrix4d& Q, Matrix4d& R){
  int m_size=4;
  for (int i=0; i<m_size; i++){
      Vec4d A_column;
      Vec4d u_column;
      for (int l=0; l<m_size; l++){
         A_column[l] = A(l,i);
         u_column[l] = A(l,i);
      }

      for (int j=0; j<i; j++){
         Vec4d q_column;
         for (int l=0; l<m_size; l++)
           q_column[l] = Q(l,j);
         double r=q_column.dot(A_column);
         R(j,i) = r;
         for (int k=0; k<m_size; k++) u_column[k] -= r*Q(k,j);
      }
      double uNorm = std::sqrt(u_column.dot(u_column));
      R(i,i)=uNorm;
      for (int k=0; k<m_size; k++)
          Q(k,i)= u_column[k]/uNorm;
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
double sphericalTriangleArea(Vec3d& a, Vec3d& b, Vec3d& c, double radius){
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

KOKKOS_INLINE_FUNCTION
Vec3d xyz_from_lat_lon(double lat, double lon, double r){
  Vec3d xyz;
  xyz[0] = r*cos(lat)*cos(lon);
  xyz[1] = r*cos(lat)*sin(lon);
  xyz[2] = r*sin(lat);
  return xyz;
}

//Rotate 90 degrees about Y axis
KOKKOS_INLINE_FUNCTION
Vec3d grid_rotation_backward(Vec3d& xyz_input){
  Vec3d xyz_output;
  xyz_output[0] =  xyz_input[2];
  xyz_output[1] =  xyz_input[1];
  xyz_output[2] = -xyz_input[0];
  return xyz_output;
}

KOKKOS_INLINE_FUNCTION
void lat_lon_from_xyz(double& lat, double& lon, Vec3d& xyz, double r){
  lon = Kokkos::atan2(xyz[1], xyz[0]);
  lat = Kokkos::asin(xyz[2]/r);
} 

}//namespace polyMPO end

#endif

