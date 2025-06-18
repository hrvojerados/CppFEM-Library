#pragma once
#include <cmath>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <utility>
#include <vector>
#include <functional>
#include <math.h>
#include <iostream>
using ll = long long;
using ull = unsigned long long;
using ld = long double;

using namespace std;

void printTriple(tuple<ld, ld, ld> t) {
  cerr << "(" 
    << get<0>(t) 
    << ", "
    << get<1>(t)
    << ", "
    << get<2>(t)
    << ")\n";
}
struct TupleHash3ull {
  size_t operator()(const tuple<ull, ull, ull>& t) const {
    size_t h1 = hash<ull>{}(get<0>(t));
    size_t h2 = hash<ull>{}(get<1>(t));
    size_t h3 = hash<ull>{}(get<2>(t));
    size_t seed = h1;
    seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

struct tupleHash2 {
    template <typename T1, typename T2>
    size_t operator()(const tuple<T1, T2>& t) const {
        size_t h1 = hash<T1>{}(get<0>(t));
        size_t h2 = hash<T2>{}(get<1>(t));
        return h1 ^ (h2 << 1);
    }
};

class Mesh2D {
public:
  ull numberOfNodes;
  ull boundaryNodeCutOff;
  ld h;
  vector<tuple<ull, ull, ull>> elements;
  unordered_map<ull, tuple<ld, ld>> getNodeCoordinates;
  unordered_map<tuple<ld, ld>, ull, tupleHash2> getNodeIndex;
  unordered_map<ull, set<ull>> neighbours;
  
  Mesh2D() = default;

  Mesh2D(function<bool(ld, ld)> inDomain,
      tuple<ld, ld> topLeft,
      tuple<ld, ld> bottomRight,
      ld precision) {
    numberOfNodes = 0;
    h = precision;
    get<0>(topLeft) -= 2 * h;
    get<0>(bottomRight) += 2 * h ;
    get<1>(topLeft) += 2 * h ;
    get<1>(bottomRight) -= 2 * h;
    ld xLeft = get<0>(topLeft) + (h / 3);
    ld xRight = get<0>(bottomRight);
    ld yTop = get<1>(topLeft);
    ld yDown = get<1>(bottomRight);
    
    ld x = xLeft;
    ull N = 0;
    while (x < xRight) {
      N++;
      x += h;
    }
    N--;
    h = (xRight - xLeft) / N;
    ld down = (sqrt(3) / 2) * h;
    ld shiftRight = h / 2;
    unordered_map<ull, bool> domainMap;
    bool indent = false;
    x = xLeft;
    ld y = yTop;
    ull ind = 0;
    while (y >= yDown) {
      ull cnt = 0;
      while (cnt != N + 1) {
        tuple<ld, ld> p = {x, y};
        // cerr << x << " " << y << " " << ind << "\n"; 
        if (inDomain(x, y)) {   
          numberOfNodes++;
          get<0>(getNodeCoordinates[ind]) = x;
          get<1>(getNodeCoordinates[ind]) = y;
          getNodeIndex[p] = ind;
          domainMap[ind] = true;
        } 
        ind++;
        cnt++;
        x += h;
      }
      indent = !indent;
      if (!indent) x = xLeft;
      else x = xLeft + shiftRight;
      y -= down;
    }
    unordered_map<ull, bool> isBoundary;
    for (auto[node, _] : domainMap) {
      if (node % (N + 1) == N)
        continue;
      if (node > N) {
        //up
        if (node % (2 * (N + 1)) <= N) {
          //   p
          //  / \
          // /   \
          // 0---1
          ull p = node - N - 1;  
          bool isPIn = (domainMap.find(p) != domainMap.end());
          bool isNode1In = (domainMap.find(node + 1) != domainMap.end());
          if (isPIn && isNode1In) {
            tuple<ull, ull, ull> key = {p, node, node + 1};
            elements.emplace_back(key);
          } else {
            isBoundary[node] = true;
          }
        } else if (node % (2 * (N + 1)) > N) {
          ull p = node - N;
          bool isPIn = (domainMap.find(p) != domainMap.end());
          bool isNode1In = (domainMap.find(node + 1) != domainMap.end());
          if (isPIn && isNode1In) {
            tuple<ull, ull, ull> key = {p, node, node + 1};
            elements.emplace_back(key);
          } else {
            isBoundary[node] = true;
          }
        }
      }
      if (node < ind - 1 - N) {
        //down
        if (node % (2 * (N + 1)) <= N) {
          ull p = node + N + 1;  
          bool isPIn = (domainMap.find(p) != domainMap.end());
          bool isNode1In = (domainMap.find(node + 1) != domainMap.end());
          if (isPIn && isNode1In) {
            tuple<ull, ull, ull> key = {p, node, node + 1};
            elements.emplace_back(key);
          } else {
            isBoundary[node] = true;
          }
        } else if (node % (2 * (N + 1)) > N) {
          ull p = node + N + 2;
          bool isPIn = (domainMap.find(p) != domainMap.end());
          bool isNode1In = (domainMap.find(node + 1) != domainMap.end());
          if (isPIn && isNode1In) {
            tuple<ull, ull, ull> key = {p, node, node + 1};
            elements.emplace_back(key);
          } else {
            isBoundary[node] = true;
          }
        }

      }
    }

    boundaryNodeCutOff = numberOfNodes - isBoundary.size();
    unordered_map<ull, ull> substitution;
    ull i = 0;
    for (auto [node, _] : domainMap) {
      if (!isBoundary[node]) {
        // cerr << node << "\n";
        substitution[node] = i;
        i++;
      }
    } 
    if (i != boundaryNodeCutOff) {
      cerr << "hmmm i != boundaryNodeCutOff!!!\n";
      cerr
        << "i = " 
        << i
        << " != "
        << boundaryNodeCutOff
        << "\n";
    }
    for (auto [node, _ ] : isBoundary) {
      // cerr << node << "\n";
      if (isBoundary[node]) {
        substitution[node]  = i;
        i++;
      }
    }
    if (i!= numberOfNodes) {
      cerr << "hmmm i != numberOfNodes!!!\n";
      cerr 
        << "i = "
        << i
        << " numberOfNodes = "
        << numberOfNodes
        << "\n";
    }

    unordered_map<ull, tuple<ld, ld>> newGetNodeCoordinates;    
    for (auto [key, val] : getNodeCoordinates) {
      newGetNodeCoordinates[substitution[key]] = val;
    }
    std::swap(getNodeCoordinates, newGetNodeCoordinates);

    for (auto [key, val] : getNodeIndex) {
      getNodeIndex[key] = substitution[val];
    }
    //elements and neighbours
    for (ull i = 0; i < elements.size(); i++) {
      get<0>(elements[i]) = substitution[get<0>(elements[i])]; 
      get<1>(elements[i]) = substitution[get<1>(elements[i])]; 
      get<2>(elements[i]) = substitution[get<2>(elements[i])]; 
      ull e0 = get<0>(elements[i]);
      ull e1 = get<1>(elements[i]);
      ull e2 = get<2>(elements[i]);
      if (e0 < boundaryNodeCutOff) {
        if (e1 < boundaryNodeCutOff) {
          neighbours[e0].insert(e1);
        }
        if (e2 < boundaryNodeCutOff) {
          neighbours[e0].insert(e2);
        } 
      }
      if (e1 < boundaryNodeCutOff) {
        if (e0 < boundaryNodeCutOff) {
          neighbours[e1].insert(e0);
        }
        if (e2 < boundaryNodeCutOff) {
          neighbours[e1].insert(e2);
        }
      }
      if (e2 < boundaryNodeCutOff) {
        if (e0 < boundaryNodeCutOff) {
          neighbours[e2].insert(e0);
        }
        if (e1 < boundaryNodeCutOff) {
          neighbours[e2].insert(e1);
        }
      }
    }
  }

  void checkMesh() {
    auto isEquilateral = [&] (tuple<ld, ld> p1, tuple<ld, ld> p2, tuple<ld, ld> p3) {
      auto dist = [&] (tuple<ld, ld> a, tuple<ld, ld> b) {
        ld dx = get<0>(a) - get<0>(b);
        ld dy = get<1>(a) - get<1>(b);
        return hypot(dx, dy);
      };
      ld d1 = dist(p1, p2);
      ld d2 = dist(p2, p3);
      ld d3 = dist(p3, p1);
      const ld EPS = 1e-9;
      return abs(d1 - d2) < EPS &&
        abs(d2 - d3) < EPS &&
        abs(d3 - d1) < EPS;
    };
    for (auto [n0, n1, n2] : elements) {
      if (!isEquilateral(getNodeCoordinates[n0], getNodeCoordinates[n1], getNodeCoordinates[n2])) {
        cerr << "element is not euqilateral\n";
        cerr << n0 << " " << n1 << " " << n2 << "\n";
      }
    }
  }
  void print() {
    cout << "Number of nodes: " << numberOfNodes << "\n";
    cout << "boundaryNodeCutOff: " << boundaryNodeCutOff << "\n";
    cout << "precision: " << h << "\n";
    cout << "elements:\n";
    for (auto tup : elements) {
      cout << "[" << get<0>(tup) 
        << ", " << get<1>(tup)
        << ", " << get<2>(tup)
        << "],\n";
    }
    cout << "Index to coordinates and vice versa:\n";
    for (auto [ind, coordinates] : getNodeCoordinates) {
      cout << ind 
        << " -> (" << get<0>(coordinates)
        << ", " << get<1>(coordinates)
        << ") & " 
        << "(" << get<0>(coordinates)
        << ", " << get<1>(coordinates)
        << ") -> " 
        << getNodeIndex[coordinates] << "\n";
    }
    cout << "neighbours:\n";
    for (ull i = 0; i < numberOfNodes; i++){
      if (neighbours.find(i) == neighbours.end()) { 
        continue;
      }
      cout << i << ": "; 
      for (ull n : neighbours[i]) {
        cout << n << " ";
      }
      cout << "\n";
    }
  }
};


class Mesh3D {
  //One day...
};

