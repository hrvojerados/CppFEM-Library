#pragma once
#include <iostream>
#include <functional>
#include <tuple>
#include <utility>
#include "tools/mesh.hpp"
#include "tools/math.hpp"

using namespace std;
using ld = long double;
using ull = unsigned long long;

inline pair<ld, ld> applyMat(ld A[2][2], pair<ld, ld> x) {
  pair<ld, ld> result; 
  result.first = A[0][0] * x.first + A[0][1] * x.second;
  result.second = A[1][0] * x.first + A[1][1] * x.second;
  return result;
}

inline ld dot(pair<ld, ld> x, pair<ld, ld> y) {
  return x.first * y.first + x.second * y.second;
}

inline ld phi0(ld x, ld y) {
  return 1 - x - y;
}
inline ld phi1(ld x, ld y) {
  return x;
}
inline ld phi2(ld x, ld y) {
  return y;
}

SparseMatrix<6> generateStiffnessMatrix(
    ld M[2][2],
    pair <ld, ld> a,
    ld C5,
    Mesh2D& mesh) {
  SparseMatrix<6> stiffMat = SparseMatrix<6>(mesh.neighbours);
  pair<ld, ld> gradPhi0 = {-1, -1};
  pair<ld, ld> gradPhi1 = {1, 0};
  pair<ld, ld> gradPhi2 = {0, 1};
  ld gaussNodes[3][3] = {{0, 0, 1.0 / 3.0},{1, 0, 1.0 / 3.0}, {0, 1, 1.0 / 3.0}};
  for (auto [node0, node1, node2] : mesh.elements) {
    tuple<ld, ld> xy0 = mesh.getNodeCoordinates[node0];
    tuple<ld, ld> xy1 = mesh.getNodeCoordinates[node1];
    tuple<ld, ld> xy2 = mesh.getNodeCoordinates[node2];
    bool isEdge0 = node0 >= mesh.boundaryNodeCutOff;
    bool isEdge1 = node1 >= mesh.boundaryNodeCutOff;
    bool isEdge2 = node2 >= mesh.boundaryNodeCutOff;
    ld x0 = get<0>(xy0);
    ld x1 = get<0>(xy1);
    ld x2 = get<0>(xy2);
    ld y0 = get<1>(xy0);
    ld y1 = get<1>(xy1);
    ld y2 = get<1>(xy2);
    ld detJ = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)); 
    //K is transpose of inverse of J multiplied by detJ
    ld K[2][2] = {{y2 - y0, y0 - y1}, {x0 - x2, x1 - x0}};
    ld grad0grad1 = (dot(
          applyMat(M, applyMat(K, gradPhi0)),
          applyMat(K, gradPhi1)
          )) / (2 * detJ);
    ld grad1grad0 = (dot(
          applyMat(M, applyMat(K, gradPhi1)),
          applyMat(K, gradPhi0)
          )) / (2 * detJ);
    ld grad1grad2 = (dot(
          applyMat(M, applyMat(K, gradPhi1)),
          applyMat(K, gradPhi2)
          )) / (2 * detJ);
    ld grad2grad1 = (dot(
          applyMat(M, applyMat(K, gradPhi2)),
          applyMat(K, gradPhi1)
          )) / (2 * detJ);
    ld grad2grad0 = (dot(
          applyMat(M, applyMat(K, gradPhi2)),
          applyMat(K, gradPhi0)
          )) / (2 * detJ);
    ld grad0grad2 = (dot(
          applyMat(M, applyMat(K, gradPhi0)),
          applyMat(K, gradPhi2)
          )) / (2 * detJ);
    ld grad0grad0 = (dot(
          applyMat(M, applyMat(K, gradPhi0)),
          applyMat(K, gradPhi0)
          )) / (2 * detJ);
    ld grad1grad1 = (dot(
          applyMat(M, applyMat(K, gradPhi1)),
          applyMat(K, gradPhi1)
          )) / (2 * detJ);
    ld grad2grad2 = (dot(
          applyMat(M, applyMat(K, gradPhi2)),
          applyMat(K, gradPhi2)
          )) / (2 * detJ);

    ld grad0Phi1 = 0;
    for (int i = 0; i < 3; i++) {
      grad0Phi1 += gaussNodes[i][2] * phi1(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad0Phi1 *= dot(applyMat(K, gradPhi0), a);

    ld grad1Phi0 = 0;
    for (int i = 0; i < 3; i++) {
      grad1Phi0 += gaussNodes[i][2] * phi0(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad1Phi0 *= dot(applyMat(K, gradPhi1), a);

    ld grad1Phi2 = 0;
    for (int i = 0; i < 3; i++) {
      grad1Phi2 += gaussNodes[i][2] * phi2(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad1Phi2 *= dot(applyMat(K, gradPhi1), a);

    ld grad2Phi1 = 0;
    for (int i = 0; i < 3; i++) {
      grad2Phi1 += gaussNodes[i][2] * phi1(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad2Phi1 *= dot(applyMat(K, gradPhi2), a);

    ld grad0Phi2 = 0;
    for (int i = 0; i < 3; i++) {
      grad0Phi2 += gaussNodes[i][2] * phi2(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad0Phi2 *= dot(applyMat(K, gradPhi0), a);

    ld grad2Phi0 = 0;
    for (int i = 0; i < 3; i++) {
      grad2Phi0 += gaussNodes[i][2] * phi0(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad2Phi0 *= dot(applyMat(K, gradPhi2), a);

    ld grad0Phi0 = 0;
    for (int i = 0; i < 3; i++) {
      grad0Phi0 += gaussNodes[i][2] * phi0(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad0Phi0 *= dot(applyMat(K, gradPhi0), a);

    ld grad1Phi1 = 0;
    for (int i = 0; i < 3; i++) {
      grad1Phi1 += gaussNodes[i][2] * phi1(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad1Phi1 *= dot(applyMat(K, gradPhi1), a);

    ld grad2Phi2 = 0;
    for (int i = 0; i < 3; i++) {
      grad2Phi2 += gaussNodes[i][2] * phi2(gaussNodes[i][0], gaussNodes[i][1]);
    }
    grad2Phi2 *= dot(applyMat(K, gradPhi2), a);

    ld phi0phi1 = 0;
    for (int i = 0; i < 3; i++) {
      phi0phi1 += gaussNodes[i][2] * 
        phi0(gaussNodes[i][0], gaussNodes[i][1]) *
        phi1(gaussNodes[i][0], gaussNodes[i][1]);
    }
    phi0phi1 *= (C5 * detJ);

    ld phi1phi2 = 0;
    for (int i = 0; i < 3; i++) {
      phi1phi2 += gaussNodes[i][2] * 
        phi1(gaussNodes[i][0], gaussNodes[i][1]) *
        phi2(gaussNodes[i][0], gaussNodes[i][1]);
    }
    phi1phi2 *= (C5 * detJ);

    ld phi0phi2 = 0;
    for (int i = 0; i < 3; i++) {
      phi0phi2 += gaussNodes[i][2] * 
        phi0(gaussNodes[i][0], gaussNodes[i][1]) *
        phi2(gaussNodes[i][0], gaussNodes[i][1]);
    }
    phi0phi2 *= (C5 * detJ);

    ld phi0phi0 = 0;
    for (int i = 0; i < 3; i++) {
      phi0phi0 += gaussNodes[i][2] * 
        phi0(gaussNodes[i][0], gaussNodes[i][1]) *
        phi0(gaussNodes[i][0], gaussNodes[i][1]);
    }
    phi0phi0 *= (C5 * detJ);

    ld phi1phi1 = 0;
    for (int i = 0; i < 3; i++) {
      phi1phi1 += gaussNodes[i][2] * 
        phi1(gaussNodes[i][0], gaussNodes[i][1]) *
        phi1(gaussNodes[i][0], gaussNodes[i][1]);
    }
    phi1phi1 *= (C5 * detJ);

    ld phi2phi2 = 0;
    for (int i = 0; i < 3; i++) {
      phi2phi2 += gaussNodes[i][2] * 
        phi2(gaussNodes[i][0], gaussNodes[i][1]) *
        phi2(gaussNodes[i][0], gaussNodes[i][1]);
    }
    phi2phi2 *= (C5 * detJ);

    if (!isEdge0 && !isEdge1) {
      stiffMat.inc(
          node0,
          node1,
          grad1grad0 + grad1Phi0 + phi0phi1);
      stiffMat.inc(
          node1,
          node0,
          grad0grad1 + grad0Phi1 + phi0phi1);
    }
    if (!isEdge1 && !isEdge2) {
      stiffMat.inc(
          node1,
          node2,
          grad2grad1 + grad2Phi1 + phi1phi2);
      stiffMat.inc(
          node2,
          node1,
          grad1grad2 + grad1Phi2 + phi1phi2);
    }
    if (!isEdge0 && !isEdge2) {
      stiffMat.inc(
          node0,
          node2,
          grad2grad0 + grad2Phi0 + phi0phi2);
      stiffMat.inc(
          node2,
          node0,
          grad0grad2 + grad0Phi2 + phi0phi2);
    }
    if (!isEdge0) {
      stiffMat.inc(
          node0,
          node0,
          grad0grad0 + grad0Phi0 + phi0phi0);
    }
    if (!isEdge1) {
      stiffMat.inc(
          node1,
          node1,
          grad1grad1 + grad1Phi1 + phi1phi1);
    }
    if (!isEdge2) {
      stiffMat.inc(
          node2,
          node2,
          grad2grad2 + grad2Phi2 + phi2phi2);
    }
  }

  return stiffMat;
}

Vector generateFreeTerm(
    function<ld(ld, ld)> f, 
    Mesh2D mesh) {
  ull n = mesh.boundaryNodeCutOff;
  Vector result = Vector(n);
  ld gaussNodes[3][3] = {{0, 0, 1.0 / 3.0},{1, 0, 1.0 / 3.0}, {0, 1, 1.0 / 3.0}};
  for (auto [node0, node1, node2] : mesh.elements) {
    tuple<ld, ld> xy0 = mesh.getNodeCoordinates[node0];
    tuple<ld, ld> xy1 = mesh.getNodeCoordinates[node1];
    tuple<ld, ld> xy2 = mesh.getNodeCoordinates[node2];
    bool isEdge0 = node0 >= mesh.boundaryNodeCutOff;
    bool isEdge1 = node1 >= mesh.boundaryNodeCutOff;
    bool isEdge2 = node2 >= mesh.boundaryNodeCutOff;
    ld x0 = get<0>(xy0);
    ld x1 = get<0>(xy1);
    ld x2 = get<0>(xy2);
    ld y0 = get<1>(xy0);
    ld y1 = get<1>(xy1);
    ld y2 = get<1>(xy2);
    ld detJ = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)); 
    auto subsX = [&] (ld xHat, ld yHat) {
      return (x1 - x0) * xHat + (x2 - x0) * yHat + x0; 
    }; 
    auto subsY = [&] (ld xHat, ld yHat) {
      return (y1 - y0) * xHat + (y2 - y0) * yHat + y0; 
    }; 

    if (!isEdge0) {
      ld b0 = 0;
      for (int i = 0; i < 3; i++) {
        ld x = subsX(gaussNodes[i][0], gaussNodes[i][1]);
        ld y = subsY(gaussNodes[i][0], gaussNodes[i][1]);
        b0 += gaussNodes[i][2] * phi0(gaussNodes[i][0], gaussNodes[i][1]) * f(x, y) * detJ;
      }
      result[node0] += b0;
    }
    if (!isEdge1) {
      ld b1 = 0;
      for (int i = 0; i < 3; i++) {
        ld x = subsX(gaussNodes[i][0], gaussNodes[i][1]);
        ld y = subsY(gaussNodes[i][0], gaussNodes[i][1]);
        b1 += gaussNodes[i][2] * phi1(gaussNodes[i][0], gaussNodes[i][1]) * f(x, y) * detJ;
      }
      result[node1] += b1;
    }
    if (!isEdge2) {
      ld b2 = 0;
      for (int i = 0; i < 3; i++) {
        ld x = subsX(gaussNodes[i][0], gaussNodes[i][1]);
        ld y = subsY(gaussNodes[i][0], gaussNodes[i][1]);
        b2 += gaussNodes[i][2] * phi2(gaussNodes[i][0], gaussNodes[i][1]) * f(x, y) * detJ;
      }
      result[node2] += b2;
    }
  }
  return result;
}

Vector solve(
    function<bool(ld, ld)> inDomain,
    tuple<ld, ld> topLeft,
    tuple<ld, ld> bottomRight,
    ld precision,
    ld M[2][2],
    pair<ld, ld> a,
    ld C5,
    function<ld(ld, ld)> f,
    ull numOfIterations,
    Mesh2D& returnedMesh) {
  Mesh2D mesh = Mesh2D(inDomain, topLeft, bottomRight, precision);
  vector<tuple<ull, ull, ull>> newElements;
  for (auto [n0, n1, n2] : mesh.elements) {
    if (n0 < mesh.boundaryNodeCutOff && n1 < mesh.boundaryNodeCutOff && n2 < mesh.boundaryNodeCutOff) {
      newElements.push_back({n0, n1, n2});
    }
  }
  std::swap(mesh.elements, newElements);
  returnedMesh = mesh;

  SparseMatrix<6> A =  generateStiffnessMatrix(M, a, C5, mesh);
  Vector b = generateFreeTerm(f, mesh);
  Vector x = Vector(b.size);
  solveGaussSeidel<6>(x, A, b, numOfIterations);
  return x;
}
