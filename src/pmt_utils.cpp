#include "pmt_utils.hpp"
#include <cmath>
#include <limits>

void PMT_Assert_Fail(const char* msg) {
  fprintf(stderr, "%s", msg);
  abort();
}

bool polyMpmTest::isEqual(double a, double b, double tol) {
  return std::fabs(a - b) < tol;
}
