#include "../src/tools/math.hpp"
#include "../src/tools/mesh.hpp"

#include <iostream>
#include <boost/container/small_vector.hpp>
#include <unordered_map>

using namespace std;
int main() {
  // unordered_map<ull, boost::container::small_vector<ull, 6>> neighbours = {
  // {0, {2, 4}},
  // {1, {2, 3}},
  // {2, {0, 1, 4}},
  // {3, {1, 4}},
  // {4, {0, 2, 3}}
  // };
  // SparseMatrix<6>* A = new SparseMatrix<6>(neighbours);
  // //A->print();
  // A->set(0, 0, 2);
  // A->set(0, 2, 1);
  // A->set(0, 4, 5);
  // A->set(1, 1, 1);
  // A->set(1, 2, 1);
  // A->set(1, 3, 1);
  // A->set(2, 0, 1.2);
  // A->set(2, 1, 5.5);
  // A->set(2, 2, 3);
  // A->set(2, 4, 2);
  // A->set(3, 1, 1);
  // A->set(3, 3, 3);
  // A->set(3, 4, 2);
  // A->set(4, 0, 2.2);
  // A->set(4, 2, 3.1);
  // A->set(4, 3, 1);
  // A->set(4, 4, 9);
  // A->print();
  unordered_map<ull, set<ull>> neighbours = {
    {0, {}},
    {1, {}},
    {2, {}},
    {3, {}},
    {4, {}},
  };
  
  SparseMatrix<6>* A = new SparseMatrix<6>(neighbours);
  A->setVal(0, 0, 0.9);
  A->setVal(1, 1, 0.9);
  A->setVal(2, 2, 0.9);
  A->setVal(3, 3, 0.9);
  A->setVal(4, 4, 0.9);

  Vector b = {3.5, 2.15, 6.28, 6.69, 10};
  Vector x = Vector(5);
  A->print();
  solveGaussSeidel(x, *A, b, 100);
  x.print();
  
  printf("%Lf\n", A->get(0, 1));
}
