#include "../src/tools/math.hpp"
#include "../src/tools/mesh.hpp"
#include "../src/assembly.hpp"
#include <cmath>
#include <iostream>
using namespace std;

int main() {
  auto inDomain =  [&] (ld x, ld y) {
    return x*x + y*y <= 1;
  };
  auto f = [&] (ld x, ld y) {
    return 4;
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
      {-1, 1},
      {1, -1},
      0.01,
      M,
      {0, 0},
      0.0,
      f,
      1000,
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
