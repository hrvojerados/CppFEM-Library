//#include <bits/stdc++.h>
#pragma once

#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <stdexcept>

using namespace std;
using ll = long long;
using u = unsigned;
using ull = unsigned long long;
using ld = long double;


class Vector {
public:
  ld* V;
  ull size;
  Vector(ull size) {
    this->size = size;
    V = new ld[size];
    for (ull i = 0; i < size; i++) {
       V[i] = 0;
    }
  }
  Vector(initializer_list<ld> init) {
    size = init.size();
    V = new ld[size];
    u i = 0;
    for (ld val : init) {
        if (i >= size) break;
        V[i++] = val;
    }
    while (i < size) V[i++] = 0;
}

  Vector(const Vector& other) {
    size = other.size;
    V = new ld[size];
    for (u i = 0; i < size; i++) {
      V[i] = other.V[i];
    }
  }

  Vector(Vector&& other) noexcept {
    V = other.V;
    other.V = nullptr;
  }

  Vector& operator=(const Vector& other) {
    if (this == &other)
      return *this;
    delete[] V;
    V = new ld[size];
    for (u i = 0; i < size; i++) {
      V[i] = other[i];
    }
    return *this;
  }

  Vector& operator=(Vector&& other) noexcept {
    if (this != &other) {
      delete[] V;
      V = other.V;
      other.V = nullptr;
    }
    return *this;
  }

  ~Vector() {
    delete[] V;
  }

  Vector operator+(const Vector& other) const {
    Vector result = Vector(size);
    for (u i = 0; i < size; i++) {
      result.V[i] = V[i] + other.V[i];
    } 
    return result;
  }

  Vector operator*(const ld scalar) const {
    Vector result = Vector(size);
    for (u i = 0; i < size; i++) {
      result.V[i] = V[i] * scalar;
    }
    return result;
  }

  const ld& operator[](u i) const {
    return V[i];
  }

  ld& operator[](u i) {
    return V[i];
  }

  void print() {
    cerr << "v = (";
    for (u i = 0; i < size; i++) {
      cerr << V[i] << ", ";
    }
    cerr << ")\n";
  }

};

ld dot(const Vector& V1, const Vector& V2) {
  ld result = 0;
  ull size = V1.size;
  for (u i = 0; i < size; i++) {
    result += V1[i] * V2[i];
  }
  return result;
}

Vector operator*(const ld scalar, const Vector& vec) {
  Vector result = Vector(vec.size);
  for (ull i = 0; i < vec.size; i++) 
    result[i] = scalar * result[i];
  return result;  
}

template<u meshBranchFactor>
class SparseMatrix {
//CSR format matrix
//NB only square matrix
public:
  ull numOfElements;
  ull numOfRows;
  ld* A;
  ull* C;
  ull* R;

  SparseMatrix( 
      unordered_map<ull, set<ull>> neigbours) {
      if (neigbours.size() == 0) { 
        throw invalid_argument("SparseMatrix requires non-empty neighbour map.");
      }
      numOfRows = neigbours.size();
      R = new ull[numOfRows + 1]; //+1 for easier loop at the end
      R[0] = 0;
      numOfElements = neigbours[0].size() + 1;
      ull i;
      for (i = 1; i < numOfRows; i++) {
        R[i] = numOfElements; 
        numOfElements += neigbours[i].size() + 1;
      }
      R[i] = numOfElements;
      A = new ld[numOfElements];
      C = new ull[numOfElements];
      for (ull i = 0; i < numOfElements; i++) {
        A[i] = 0;
      }
      
      for (ull i = 0; i < numOfRows; i++) {
        ull j = R[i];
        bool isIInserted = false;
        for (ull n : neigbours[i]) {
          if (!isIInserted && n > i) {
            C[j] = i;
            isIInserted = true;
            j++;
          }
          C[j] = n;
          j++;
        }
        if (!isIInserted) {
          C[j] = i;
        }
      }
  }

  inline ld get(ull i, ull j) {
    //for now linear search
    ull it = R[i];
    while (it < R[i + 1] - 1 && C[it] < j) it++;
    if (C[it] == j) {
      return A[it];
    }
    return 0;
  }

  inline void setVal(ull i, ull j, ld val) { 
    //for now linear search
    ull it = R[i];
    while (it < R[i + 1] - 1 && C[it] < j) it++;
    if (C[it] == j) {
      A[it] = val;
    }
  }

  inline void inc(ull i, ull j, ld val) {
    //for now linear search
    ull it = R[i];
    while (it < R[i + 1] - 1 && C[it] < j) it++;
    if (C[it] == j) {
      A[it] += val;
    }
  }

  inline void print() {
    cerr << "A:\n";
    for (ull i = 0; i < numOfElements; i++) {
      cerr << A[i] << " ";
    }
    cerr << "\n";

    cerr << "C:\n";
    for (ull i = 0; i < numOfElements; i++) {
      cerr << C[i] << " ";
    }
    cerr << "\n";

    cerr << "R:\n";
    for (ull i = 0; i < numOfRows; i++) {
      cerr << R[i] << " ";
    }
    cerr << "\n";

    for (ull i = 0; i < numOfRows; i++) {
      for (ull j = 0; j < numOfRows; j++) {
        cerr << get(i, j) << " ";
      }
      cerr << "\n";
    }
  }

  ~SparseMatrix() {//TODO no rule of 5
    delete[] A;
    delete[] C;
    delete[] R;
  }
};

pair<ld, ld> grad(ld x, ld y, function<ld(ld, ld)> f, ld eps) {
  pair<ld, ld> result;
  result.first = (f(x + eps, y) - f(x, y)) / eps;
  result.second = (f(x, y + eps) - f(x, y)) / eps;
  return result;
} 

template <u meshBranchFactor>
void solveGaussSeidel(
    Vector &x,
    SparseMatrix<meshBranchFactor> &M,
    Vector &b,
    u numOfIterations) {
  ull n = b.size;
  if (M.numOfRows != n) {
    cerr << M.numOfRows << " " << n << "\n";
    throw invalid_argument("dimension error in Gauss-Seidel solver");
  }
  for (u i = 0; i < n; i++) x[i] = 0;

  for (u k = 1; k <= numOfIterations; k++) {
    for (u i = 0; i < M.numOfRows; i++) {
      ld sumOfLessThan_i = 0;
      u j;
      for (j = M.R[i];
          j < M.R[i + 1] && M.C[j] < i;
          j++) {
        sumOfLessThan_i += M.A[j] * x[M.C[j]];
      }
      ld sumOfGreaterThan_i = 0;
      ld Mii = 0;
      if (M.C[j] == i) {
        Mii = M.A[j];
        j++;
      } 
      for (; j < M.R[i + 1]; j++) {
        sumOfGreaterThan_i += M.A[j] * x[M.C[j]];
      }
      x[i] = 
        (b[i] - sumOfLessThan_i - sumOfGreaterThan_i) / Mii; 
    }
  }
}
template <u meshBranchFactor>
void solveCG(
    Vector &x,
    SparseMatrix<meshBranchFactor> &M,
    Vector &b,
    u numOfIterations) 
{
  ull n = b.size;
  if (M.numOfRows != n) {
    cerr << M.numOfRows << " " << n << "\n";
    throw invalid_argument("dimension error in CG solver");
  }
  Vector r(n), p(n), Ap(n);
  for (ull i = 0; i < n; ++i) 
    x[i] = 0;
  for (ull i = 0; i < n; ++i) 
    r[i] = b[i];
  for (ull i = 0; i < n; ++i) 
    p[i] = r[i];
  long double rsold = dot(r, r);
  for (ull k = 0; k < numOfIterations; k++) {
    for (ull i = 0; i < n; ++i) {
      ld sum = 0;
      for (ull j = M.R[i]; j < M.R[i+1]; j++)
        sum += M.A[j] * p[M.C[j]];
      Ap[i] = sum;
    }
    ld pAp = dot(p, Ap);
    if (pAp == 0) break; 
    ld alpha = rsold / pAp;
    for (ull i = 0; i < n; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * Ap[i];
    }
    ld rsnew = dot(r, r);
    if (sqrt(rsnew) < 1e-10L) 
      break;

    ld beta = rsnew / rsold;
    for (ull i = 0; i < n; ++i) {
      p[i] = r[i] + beta * p[i];
    }
    rsold = rsnew;
  }
}


