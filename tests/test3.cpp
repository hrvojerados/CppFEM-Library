#include "../src/tools/math.hpp"
#include "../src/tools/mesh.hpp"
#include "../src/assembly.hpp"

using namespace std;
int main() {
  unordered_map<ull, set<ull>> neighbours = {
  {0, {2, 4}},
  {1, {2, 3}},
  {2, {0, 1, 4}},
  {3, {1, 4}},
  {4, {0, 2, 3}}
  };
  SparseMatrix<6>* A = new SparseMatrix<6>(neighbours);
  //A->print();
  A->setVal(0, 0, 2);
  A->setVal(0, 2, 1);
  A->setVal(0, 4, 5);
  A->setVal(1, 1, 1);
  A->setVal(1, 2, 1);
  A->setVal(1, 3, 1);
  A->setVal(2, 0, 1.2);
  A->setVal(2, 1, 5.5);
  A->setVal(2, 2, 3);
  A->setVal(2, 4, 2);
  A->setVal(3, 1, 1);
  A->setVal(3, 3, 3);
  A->setVal(3, 4, 2);
  A->setVal(4, 0, 2.2);
  A->setVal(4, 2, 3.1);
  A->setVal(4, 3, 1);
  A->setVal(4, 4, 9);
  A->print();
}
