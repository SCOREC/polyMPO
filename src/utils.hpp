#ifndef POLYMPMTEST_UTILS_H
#define POLYMPMTEST_UTILS_H


namespace polyMpmTest{

using VectorView = Kokkos::View<Vector*>;
using DoubleView = Kokkos::View<double*>;
using IntView = Kokkos::View<int*>;
using BoolView = Kokkos::View<bool*>;

class Vector {
  private:
    double coords_[2];

  public:
    Vector():coords_{0.0, 0.0}{}
    ~Vector() = default;
    Vector(double x1, double x2):coords_{x1,x2}{}
    Vector(double x[2]):coords_{x[0], x[1]}{}

    double operator[](int i) const;
    double &operator[](int i);

    Vector operator-();
    Vector operator+(Vector v);
    Vector operator-(Vector v);
    Vector operator*(double scalar) const;
    Vector operator*(Vector v) const;
    Vector &operator+=(const Vector &v);
    void operator=(const double &scalar);
    double dot(Vector v);
    double cross(Vector v);
    double magnitude();
   
    


};


void initArray(Vector arr*, int n, Vector fill){
}

void initArray(double arr*, int n, double fill)}{
}

}//namespace polyMpmTest end

#endif

