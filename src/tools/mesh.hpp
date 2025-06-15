#pragma once
#include <cmath>
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
  unordered_map<ull, vector<ull>> neighbours;
  
  Mesh2D() = default;
  Mesh2D(function<bool(ld, ld)> inDomain,
      tuple<ld, ld> topLeft,
      tuple<ld, ld> bottomRight,
      ld precision) {
    h = precision;
    ld xLeft = get<0>(topLeft);
    ld xRight = get<0>(bottomRight);
    ld yTop = get<1>(topLeft);
    ld yDown = get<1>(bottomRight);
    
    ld x = xLeft;
    ull N = 0;
    while (x < xRight) {
      N++;
      x += h;
    }
    h = (xRight - xLeft) / N;
    ld down = (sqrt(3) / 2) * h;

    ld shiftRight = h / 2;

    unordered_map<ull, bool> domainMap;
    x = xLeft;
    ld y = yTop;
    ull ind = 0;
    while (y >= yDown) {
      while (x <= xRight) {
        tuple<ld, ld> p = {x, y};
        if (inDomain(x, y)) {   
          numberOfNodes++;
          get<0>(getNodeCoordinates[ind]) = x;
          get<1>(getNodeCoordinates[ind]) = y;
          getNodeIndex[p] = ind;
          domainMap[ind] = true;
        } 
        ind++;
        x += h;
      }
      x = xLeft;
      y -= down;
    }
    /*    2_____3
     *   /       \
     *  1        4
     *  \       /
     *   6_____5
     *   pretpostavka da je broj nodeova manji od 2^63
     */
    //uvijek ih sortirane ubacuj
    // p2 < p3 < p1 < node < p4 < p6 < p5
    unordered_map<ull, bool> isBoundary;
    unordered_map<tuple<ull, ull, ull>, bool, TupleHash3ull> addedElements;
    for (auto[node, _] : domainMap) {
      ull p1 = node - 1;
      ull p4 = node + 1;
      ull p2;
      if (node % (N + 1) <= N) p2 = node - N - 2;
      else p2 = node - N - 1;
      ull p3 = p2 + 1;
      ull p6;
      if (node % (N + 1) <= N) p6 = node + N;
      else p6 = node + N + 1;
      ull p5 = p6 + 1;

      bool isP1In = (domainMap.find(p1) != domainMap.end());
      bool isP2In = (domainMap.find(p2) != domainMap.end());
      bool isP3In = (domainMap.find(p3) != domainMap.end());
      bool isP4In = (domainMap.find(p4) != domainMap.end());
      bool isP5In = (domainMap.find(p5) != domainMap.end());
      bool isP6In = (domainMap.find(p6) != domainMap.end());

      if (isP1In && isP2In) {
        tuple<ull, ull, ull> key = {p2, p1, node};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }
      if (isP2In && isP3In) {
        tuple<ull, ull, ull> key = {p2, p3, node};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }

      if (isP3In && isP4In) {
        tuple<ull, ull, ull> key = {p3, node, p4};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }

      if (isP4In && isP5In) {
        tuple<ull, ull, ull> key = {node, p4, p5};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }

      if (isP5In && isP6In) {
        tuple<ull, ull, ull> key = {node, p6, p5};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }

      if (isP6In && isP1In) {
        tuple<ull, ull, ull> key = {p1, node, p6};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }
    }

    boundaryNodeCutOff = isBoundary.size();
    unordered_map<ull, ull> substitution;
    ull i = 0;
    for (auto [node, _ ] : isBoundary) {
      substitution[node]  = i;
      i++;
    }
    for (auto [node, _] : domainMap) {
      if (!isBoundary[node]) {
        substitution[node] = i;
        i++;
      }
    } 
    if (i + 1 != numberOfNodes)
      cerr << "hmmm i + 1 != numberOfNodes!!!\n";

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
      neighbours[e0].push_back(e1);
      neighbours[e0].push_back(e2);
      neighbours[e1].push_back(e0);
      neighbours[e1].push_back(e2);
      neighbours[e2].push_back(e0);
      neighbours[e2].push_back(e1);
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

