#include "../src/tools/math.hpp"
#include "../src/tools/mesh.hpp"
#include "../src/assembly.hpp"
#include <cmath>
#include <iostream>
using namespace std;

int main() {
  auto inDomain =  [&] (ld x, ld y) {
    return (x >= 0) && (x <= 1) && (y >= 0) && (y <= 1);
  };
  auto f = [&] (ld x, ld y) {
    return 2 * (M_PI * M_PI) * sin(M_PI * x) * sin(M_PI * y);
  };
  ld M[2][2];
  M[0][0] = 1;
  M[0][1] = 0;
  M[1][0] = 0;
  M[1][1] = 1;
  Mesh2D mesh;
  ull numOfElementsToGraph;
  Vector sol = solve(
      inDomain,
      {0, 1},
      {1, 0},
      0.0025,
      M,
      {0, 0},
      0.0,
      f,
      100,
      mesh,
      numOfElementsToGraph);
  cout << sol.size << " " << numOfElementsToGraph << "\n";
  for (ull i = 0; i < sol.size; i++) {
    cout << get<0>(mesh.getNodeCoordinates[i]) 
      << " "
      << get<1>(mesh.getNodeCoordinates[i])
      << " "
      << sol[i]
      << "\n";
  } 
  for (auto [n0, n1, n2] : mesh.elements) {
    if (n0 < mesh.boundaryNodeCutOff && n1 < mesh.boundaryNodeCutOff && n2 < mesh.boundaryNodeCutOff) {
      cout << n0 << " " << n1 << " " << n2 << "\n";
    }
  }

}
