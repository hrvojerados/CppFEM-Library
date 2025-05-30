#include <bits/stdc++.h>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <sys/types.h>
#include <unordered_map>
#include <boost/container/small_vector.hpp>
using namespace std;
using ll = long long;
using u = unsigned;
using ull = unsigned long long;
using ld = long double;

class Mesh2D {
public:
  tuple<ll, ll, ll> *Elements;
  ll NumberOfElements;
  unordered_map<ll, tuple<double, double>> getNodeCoordinates;
};
class Mesh3D {};
//CSR format matrix
template<u meshBranchFactor>
class SparseMatrix {
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

  ~SparseMatrix() {
    delete[] A;
    delete[] C;
    delete[] R;
  }
};

int main() {
  unordered_map<ull, boost::container::small_vector<ull, 6>> neighbours = {
    {0, {2, 4}},
    {1, {2, 3}},
    {2, {0, 1, 4}},
    {3, {1, 4}},
    {4, {0, 2, 3}}
  };
  SparseMatrix<6>* A = new SparseMatrix<6>(neighbours);
  //A->print();
  A->set(4, 3, 0.69);
  A->set(0, 2, 3.14);
  A->print();
  printf("%Lf\n", A->get(0, 1));
}
