#include "../src/tools/math.hpp"
#include "../src/tools/mesh.hpp"

using namespace std;

int main() {
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
