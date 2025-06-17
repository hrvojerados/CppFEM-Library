#include "../src/tools/math.hpp"
#include "../src/tools/mesh.hpp"

#include <unordered_map>

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
  A->setVal(4, 3, 0.69);
  A->setVal(0, 2, 3.14);
  A->print();
  printf("%Lf\n", A->get(0, 1));
  Vector v = {1.2, 3.14, 2.71};
  Vector u = (2 * v);

  Vector x = u + v;
  x.print();
  v.print();
  u.print();
  Mesh2D mesh;

  
}
