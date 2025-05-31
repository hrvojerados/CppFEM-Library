#pragma once
#include <cmath>
#include <tuple>
#include <unordered_map>
#include <set>
#include <utility>
#include <vector>
#include <functional>
#include <math.h>
#include <boost/container/small_vector.hpp>
 
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
  vector<tuple<ll, ll, ll>> elements;
  unordered_map<ull, tuple<ld, ld>> getNodeCoordinates;
  unordered_map<tuple<ld, ld>, ull, tupleHash2> getNodeIndex;
  unordered_map<ull, set<ull>> neighbours;
  
  Mesh2D() = default;
  Mesh2D(function<bool(tuple<ld,ld>)> inDomain,
      tuple<ld, ld> topLeft,
      tuple<ld, ld> bottomRight,
      ld precision) {
    h = precision;
    ld xLeft = get<0>(topLeft);
    ld xRight = get<0>(bottomRight);
    ld yTop = get<1>(topLeft);
    ld yDown = get<1>(bottomRight);
    
    ld x = xLeft;
    ld N = 0;
    while (x < xRight) {
      N++;
      x += h;
    }
    h = (xRight - xLeft) / N;
    ld down = (sqrt(3) / 2) * h;

    ld shiftRight = h / 2;

    set<ld> domain;
    unordered_map<ull, bool> domainMap;
    x = xLeft;
    ld y = yTop;
    ull ind = 0;
    while (y >= yDown) {
      while (x <= xRight) {
        tuple<ld, ld> p = {x, y};
        if (inDomain(p)) {   
          numberOfNodes++;
          get<0>(getNodeCoordinates[ind]) = x;
          get<1>(getNodeCoordinates[ind]) = y;
          getNodeIndex[p] = ind;
          domain.insert(ind);
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
    for (ull node : domain) {
      ull p1 = node - 1;
      ull p2 = node - N - 1;
      ull p3 = node - N;
      ull p4 = node + 1;
      ull p5 = node + N + 1;
      ull p6 = node + N;
      if (domainMap[p1] && domainMap[p2]) {
        tuple<ull, ull, ull> key = {p2, p1, node};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }
      if (domainMap[p2] && domainMap[p3]) {
        tuple<ull, ull, ull> key = {p2, p3, node};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }

      if (domainMap[p3] && domainMap[p4]) {
        tuple<ull, ull, ull> key = {p3, node, p4};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }

      if (domainMap[p4] && domainMap[p5]) {
        tuple<ull, ull, ull> key = {node, p4, p5};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }

      if (domainMap[p5] && domainMap[p6]) {
        tuple<ull, ull, ull> key = {node, p6, p5};
        if (!addedElements[key]) {
          elements.emplace_back(key);
          addedElements[key] = true;
        }
      } else {
        isBoundary[node] = true;
      }

      if (domainMap[p6] && domainMap[p1]) {
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
    for (ull node : domain) {
      substitution[node] = i;
      i++;
      //TODO take boundary into consideration
    } 
    unordered_map<ull, tuple<ld, ld>> newGetNodeCoordinates;    
    for (auto [key, val] : getNodeCoordinates) {
      newGetNodeCoordinates[substitution[key]] = val;
    }
    getNodeCoordinates = move(newGetNodeCoordinates);

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
      neighbours[e0].insert(e1);
      neighbours[e0].insert(e2);
      neighbours[e1].insert(e0);
      neighbours[e1].insert(e2);
      neighbours[e2].insert(e0);
      neighbours[e2].insert(e1);
    }
  }
};


class Mesh3D {
  //One day...
};

