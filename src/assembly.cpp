#include <iostream>
#include <functional>
#include "tools/mesh.hpp"
#include "tools/math.hpp"

using namespace std;
using ld = long double;
using ull = unsigned long long;

inline void apply(ld** A, ld* x) {
  ld tmp = A[0][0] * (*x) + A[0][1] * (*(x + 1));
  *(x + 1) = A[1][0] * (*x) + A[1][1] * (*(x + 1));
  *x = tmp;
}

SparseMatrix<6> generateStiffnessMatrix(ld A[2][2], Mesh2D& mesh) {
  SparseMatrix<6> stiffMat = SparseMatrix<6>(mesh.neighbours);
  ld gradPhi0[2] = {-1, 1};
  ld gradPhi1[2] = {1, 0};
  ld gradPhi2[2] = {0, 1};
  for (auto [node0, node1, node2] : mesh.elements) {
    tuple<ld, ld> xy0 = mesh.getNodeCoordinates[node0];
    tuple<ld, ld> xy1 = mesh.getNodeCoordinates[node1];
    tuple<ld, ld> xy2 = mesh.getNodeCoordinates[node2];
    ld x0 = get<0>(xy0);
    ld x1 = get<0>(xy1);
    ld x2 = get<0>(xy2);
    ld y0 = get<1>(xy0);
    ld y1 = get<1>(xy1);
    ld y2 = get<1>(xy2);
    ld detJ = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0); 
    //K is transpose of inverse of J multiplied by detJ
    ld K[2][2] = {{y2 - y0, y0 - y1}, {x0 - x2, x1 - x0}};
     
  }

  return stiffMat;
}

template<ull n>
Vector<n> generateFreeTerm(
    function<ld(ld, ld)> f, 
    function<ld(ld, ld)> R,
    Vector<2> a,
    ld F) {

}

