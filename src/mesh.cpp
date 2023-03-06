#include "mesh.hpp"

Mesh::Mesh() {
  coords = new int[10];
}

Mesh::~Mesh() {
  delete [] coords;
}
