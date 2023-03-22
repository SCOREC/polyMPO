#include "pmt_utils.hpp"

void PMT_Assert_Fail(const char* msg) {
  fprintf(stderr, "%s", msg);
  abort();
}
