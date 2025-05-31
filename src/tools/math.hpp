//#include <bits/stdc++.h>
#pragma once

#include <unordered_map>
#include <algorithm>
#include <boost/container/small_vector.hpp>
#include <cstddef>
#include <cstdio>
#include <stdexcept>

using namespace std;
using ll = long long;
using u = unsigned;
using ull = unsigned long long;
using ld = long double;


template<u size>
class Vector {
public:
  ld* V;

  Vector() {
    V = new ld[size];
  }
  Vector(initializer_list<ld> init) {
    V = new ld[size];
    u i = 0;
    for (ld val : init) {
        if (i >= size) break;
        V[i++] = val;
    }
    while (i < size) V[i++] = 0;
}

  Vector(const Vector& other) {
    V = new ld[size];
    for (u i = 0; i < size; i++) {
      V[i] = other.V[i];
    }
  }

  Vector(Vector&& other) noexcept {
    V = other.V;
    other.V = nullptr;
  }

  Vector<size>& operator=(const Vector& other) {
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

  Vector<size> operator+(const Vector& other) const {
    Vector<size>  result;
    for (u i = 0; i < size; i++) {
      result.V[i] = V[i] + other.V[i];
    } 
    return result;
  }

  Vector<size> operator*(const ld scalar) const {
    Vector<size>  result;
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
    printf("v = (");
    for (u i = 0; i < size; i++) {
      printf("%Lf, ", V[i]);
    }
    printf(")\n");
  }

};

template<u size>
ld dot(const Vector<size>& V1, const Vector<size>& V2) {
  ld result = 0;
  for (u i = 0; i < size; i++) {
    result += V1[i] * V2[i];
  }
  return result;
}

template<u size>
Vector<size> operator*(const ld scalar, const Vector<size>& vec) {
  return vec * scalar;  
}

template<u meshBranchFactor>
class SparseMatrix {
//CSR format matrix
public:
  ull numOfElements;
  ull numOfRows;
  ld* A;
  ull* C;
  ull* R;

  SparseMatrix(
      unordered_map<ull, boost::container::small_vector<ull, meshBranchFactor>> neigbours) {
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
        sort(neigbours[i].begin(), neigbours[i].end());
        bool isIInserted = false;
        ull j = R[i];
        for (
            u it = 0;
            j < R[i + 1] && it < neigbours[i].size();
            j++, it++) {
          if (neigbours[i][it] > i && !isIInserted) {
            C[j] = i;
            isIInserted = true;
            j++;
          }
          C[j] = neigbours[i][it];
        }
        if (!isIInserted) {
          C[j] = i;
        } 
      }
  }

  inline ld get(ull i, ull j) {
    //for now linear search
    ull it = R[i];
    while (C[it] < j) it++;
    if (C[it] == j) {
      return A[it];
    }
    return 0;
  }

  inline void set(ull i, ull j, ld val) { 
    //for now linear search
    ull it = R[i];
    while (C[it] < j) it++;
    if (C[it] == j) {
      A[it] = val;
    }
  }

  inline void inc(ull i, ull j, ld val) {
    //for now linear search
    ull it = R[i];
    while (C[it] < j) it++;
    if (C[it] == j) {
      A[it] += val;
    }
  }

  inline void print() {
    printf("A:\n");
    for (ull i = 0; i < numOfElements; i++) {
      printf("%Lf ", A[i]);
    }
    printf("\n");

    printf("C:\n");
    for (ull i = 0; i < numOfElements; i++) {
      printf("%llu ", C[i]);
    }
    printf("\n");

    printf("R:\n");
    for (ull i = 0; i < numOfRows; i++) {
      printf("%llu ", R[i]);
    }
    printf("\n");

    for (ull i = 0; i < numOfRows; i++) {
      for (ull j = 0; j < numOfRows; j++) {
        printf("%Lf ", get(i, j));
      }
      printf("\n");
    }
  }

  ~SparseMatrix() {//TODO no rule of 5
    delete[] A;
    delete[] C;
    delete[] R;
  }
};


