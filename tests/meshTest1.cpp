#include "../src/tools/math.hpp"
#include "../src/tools/mesh.hpp"
#include "../src/assembly.hpp"
#include <cmath>
using namespace std;
using ld = long double;

bool inDomain(ld x, ld y) {
  return (x >= 0) && (x <= 8.5) && (y <= 0) && (y >= -2 * sqrt(3));
}

int main() {

  Mesh2D mesh = Mesh2D(inDomain,
      {-0.1, 0.1},
      {9.0, -4.1},
      1.0);
  mesh.print();
}
