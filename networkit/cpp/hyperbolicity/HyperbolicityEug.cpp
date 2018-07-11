/*
 * Hyperbolicity.cpp
 *
 *  Created on: 07.05.2018
 *      Author: Eugenio Angriman
 */

#include "HyperbolicityEug.h"
#include "../auxiliary/BucketPQ.h"
#include "../components/BiconnectedComponents.h"

#include <time.h>

namespace NetworKit {

Hyperbolicity::Hyperbolicity(const Graph &G)
    : G(G), n(G.upperNodeIdBound()),
      pairs(Aux::PrioQueue<int64_t, uint64_t>((uint64_t)(n * (n - 1) / 2))),
      prunedAPSP(PrunedAPSP(G)), nonFarApartPairs(CSRMatrix(n)),
      distance(CSRMatrix(n)) {
  if (G.isDirected()) {
    throw std::runtime_error(
        "Hyperbolicity can be computed on unidrected graphs only.");
  }
  if (G.numberOfSelfLoops() > 0) {
    throw std::runtime_error(
        "Hyperbolicity cannot be computed on graphs with self-loops. Call the "
        "removeSelfLoops() methods before computing the hyperbolicity on this "
        "graph");
  }
}

void Hyperbolicity::init() {
  minPivotDist.assign(n, std::numeric_limits<count>::max());
  nearestPivot.assign(n, std::numeric_limits<count>::max());
  maxPivotDist.assign(n, 0);
  degreeInCmp.assign(n, 0);
  isInCmp.assign(n, 0);
  acc.assign(n, 0);
  ecc.assign(n, 0);
  distCentral.assign(n, 0);
  isPivot.assign(n, 0);
  valuable.assign(n, 0);
  seenBool.assign(n, 0);
  sourceDist.assign(n, 0);
  visited.assign(n, n);

  hasRun = false;
  maxHyp = 0;

  while (pairs.size() > 0) {
    pairs.extractMin();
  }
}

void Hyperbolicity::run() {
  init();

  prunedAPSP.run();
  BiconnectedComponents bcc(G);
  bcc.run();
  auto compSizes = bcc.getComponentSizes();
  auto comps = bcc.getComponents();

  std::multimap<count, count> sortedComps;
  for (auto iter = compSizes.begin(); iter != compSizes.end(); ++iter) {
    // Hyperbolicity can be computed in graphs with at least 4 nodes
    if (iter->second > 3) {
      sortedComps.insert(std::make_pair(iter->second, iter->first));
    }
  }

  if (sortedComps.size() == 0) {
    WARN("This graph has no biconnected components with at least 4 nodes. The "
         "hyperbolicity is 0.");
  }

  for (auto iter = sortedComps.rbegin(); iter != sortedComps.rend(); ++iter) {
    computeHypOfComponent(comps[iter->second]);
  }

  hasRun = true;
}

void Hyperbolicity::initHypComponent(const std::vector<node> &component) {
  std::fill(isInCmp.begin(), isInCmp.end(), 0);
  std::fill(degreeInCmp.begin(), degreeInCmp.end(), 0);
  std::fill(isPivot.begin(), isPivot.end(), 0);
  std::fill(valuable.begin(), valuable.end(), 0);
  std::fill(seenBool.begin(), seenBool.end(), 0);

  seen.clear();
  seen.reserve(component.size());
  mate.clear();
  mate.resize(n);

  nPivots = numberOfPivots(component.size());
  sortedComponent.resize(nPivots, 0);

  mostCentral = n;
  maxDeg = 0;

  for (count i = 0; i < component.size(); ++i) {
    isInCmp[component[i]] = true;
  }

  auto comparator = [](std::pair<node, count> p1, std::pair<node, count> p2) {
    return p1.second == p2.second ? p1.first < p2.first : p1.second > p2.second;
  };
  std::priority_queue<std::pair<node, count>,
                      std::vector<std::pair<node, count>>, decltype(comparator)>
      degSet(comparator);

  node u;
  count degree;
  for (count i = 0; i < component.size(); ++i) {
    u = component[i];
    degree = 0;
    G.forNeighborsOf(u, [&](node v) {
      if (isInCmp[v]) {
        ++degree;
      }
    });
    degreeInCmp[u] = degree;
    if (degree > maxDeg) {
      maxDeg = degree;
      mostCentral = u;
    } else if (degree == maxDeg && mostCentral > u) {
      mostCentral = u;
    }

    degSet.push(std::make_pair(u, degree));
    if (degSet.size() > nPivots) {
      degSet.pop();
    }
  }

  // GroupDegree gd(G, nPivots, true);
  // gd.run();
  // sortedComponent = gd.groupMaxDegree();
  for (count i = nPivots; i > 0; --i) {
    sortedComponent[i - 1] = degSet.top().first;
    degSet.pop();
  }

  for (count i = 0; i < nPivots; ++i) {
    node source = sortedComponent[i];
    bfsFrom(source, isInCmp);
    isPivot[source] = true;
  }

  clock_t t = -clock();
  storePairs(component, isInCmp);
  t += clock();
  INFO("Store pairs time = ", (double)t / CLOCKS_PER_SEC);
}

void Hyperbolicity::computeHypOfComponent(const std::vector<node> &component) {
  initHypComponent(component);
  INFO("Pivots = ", nPivots);
  INFO("Pairs = ", pairs.size());
  distance.sort();
  count numPairs = 0;
  delta = 0;
  node x, y;

  clock_t pairTime, deltaTime, first, second;
  std::pair<int64_t, uint64_t> top;
  std::pair<node, node> pair;

  while (pairs.size() > 0) {
    std::fill(acc.begin(), acc.end(), 0);
    clock_t t = clock();
    top = pairs.extractMin();
    pair = pairFromIndex(top.second);

    x = pair.first;
    y = pair.second;

    dist_xy = -top.first;

    // Checking distance upper bound
    if (dist_xy <= delta) {
      INFO("Distance ", -top.first, " <= ", delta, " aborting");
      break;
    }

    numPairs++;

    mate[x].push_back(y);
    mate[y].push_back(x);

    if (!seenBool[x]) {
      seen.push_back(x);
      seenBool[x] = true;
    }
    if (!seenBool[y]) {
      seen.push_back(y);
      seenBool[y] = true;
    }

    deltaPlusOne = delta + 1;
    condVal = dist_xy - delta;
    // TODO swap x and y if y is more central than x.

    if (3 * deltaPlusOne >= 2 * dist_xy) {
      condAcc1 = 3 * deltaPlusOne - 2 * dist_xy;
      condAcc2 = 2 * deltaPlusOne - dist_xy;

      for (node v : seen) {
        if (x == v || y == v) {
          continue;
        }

        clock_t t1 = clock();
        dist_yv = getDistanceUnknown(y, v);
        dist_xv = getDistanceUnknown(x, v);
        first += clock() - t1;
        if (2 * dist_xv >= deltaPlusOne && 2 * dist_yv >= deltaPlusOne) {
          ecc_v = eccUpperBound(v);
          if (2 * (ecc_v - dist_xv) >= condAcc1) {
            if (2 * (ecc_v - dist_yv) >= condAcc1) {
              if (2 * ecc_v >= condAcc2 + dist_xv + dist_yv) {
                acc[v] = true;
                if (2 * distCentral[v] + condVal > dist_xv + dist_yv) {
                  valuable[v] = true;
                }
              }
            }
          }
        }
      }
    } else {
      for (node v : component) {
        if (x != v && y != v) {
          clock_t t1 = clock();
          dist_yv = getDistanceUnknown(y, v);
          dist_xv = getDistanceUnknown(x, v);
          second += clock() - t1;
          if (2 * dist_xv >= deltaPlusOne && 2 * dist_yv >= deltaPlusOne) {
            ecc_v = eccUpperBound(v);
            if (2 * ecc_v + dist_xy >= 2 * deltaPlusOne + dist_xv + dist_yv) {
              acc[v] = true;
              if (2 * distCentral[v] + condVal > dist_xv + dist_yv) {
                valuable[v] = true;
              }
            }
          }
        }
      }
    }

    pairTime += clock() - t;
    t = clock();
    node v;
    for (count i = 0; i < component.size(); ++i) {
      v = component[i];
      if (x != v && y != v && valuable[v]) {
        updateDelta(x, y, v);
      }
    }
    deltaTime += clock() - t;
  }

  INFO("Pair time  = ", pairTime);
  INFO("Delta time = ", deltaTime);
  INFO("First      = ", first);
  INFO("Second     = ", second);
  INFO("Considered pairs = ", numPairs);
}

void Hyperbolicity::updateDelta(node x, node y, node v) {
  for (node w : mate[v]) {
    if (acc[w] && x != w && y != w) {
      dist_xv = getDistanceUnknown(x, v);
      dist_yv = getDistanceUnknown(y, v);
      dist_yw = getDistanceUnknown(y, w);
      dist_xw = getDistanceUnknown(x, w);
      dist_vw = getDistanceUnknown(v, w);
      s2 = dist_xv + dist_yw; // dist_xv + prunedAPSP.getDistance(y, w);
      s3 = dist_xw + dist_yv;
      // prunedAPSP.getDistance(x, w) + dist_yv;
      diff = (int64_t)dist_xy +
             dist_vw - //(int64_t)prunedAPSP.getDistance(v, w) -
             std::max(s2, s3);
      if (diff > (int64_t)delta) {
        delta = diff;
        if (delta > maxHyp) {
          maxHyp = delta;
          solX = x;
          solY = y;
          solV = v;
          solW = w;
        }
      }
    }
  }
}

void Hyperbolicity::bfsFrom(node source, std::vector<bool> &isInCmp) {
  count d;
  node u;
  node mostFar = source;
  // TODO We could save memory by indexing the component nodes!
  visited[source] = source;
  minPivotDist[source] = 0;
  sourceDist[source] = 0;
  std::queue<node> queue;
  queue.push(source);

  while (!queue.empty()) {
    u = queue.front();
    queue.pop();
    d = sourceDist[u] + 1;
    G.forNeighborsOf(u, [&](node v) {
      if (isInCmp[v]) {
        if (visited[v] != source) {
          sourceDist[v] = d;
          mostFar = v;
          if (source == mostCentral) {
            distCentral[v] = d;
          }
          maxPivotDist[v] = std::max(maxPivotDist[v], d);
          if (d < minPivotDist[v]) {
            minPivotDist[v] = d;
            nearestPivot[v] = source;
          }
          queue.push(v);
          visited[v] = source;

          setNonFarApart(u, v);
        } else if (sourceDist[v] == sourceDist[u] + 1) {
          setNonFarApart(u, v);
        }
      }
    });
  }

  ecc[source] = sourceDist[mostFar];
}

void Hyperbolicity::storePairs(const std::vector<node> &component,
                               const std::vector<bool> &isInCmp) {
  count d, eccLower;
  node u, v;
  for (count i = 0; i < component.size() - 1; ++i) {
    for (count j = i + 1; j < component.size(); ++j) {
      u = component[i];
      v = component[j];

      if (!G.hasEdge(u, v) && nonFarApartPairs(u, v) == 0) {
        eccLower = std::max(eccLowerBound(u), eccLowerBound(v));
        d = distUpperBound(u, v);
        if (d >= eccLower) {
          d = prunedAPSP.getDistance(u, v);
          if (d >= eccLower) {
            pairs.insert(-d, pairIndex(u, v));
            distance.setValue(u, v, d);
          }
        }
      }
    }
  }
}
} // namespace NetworKit
