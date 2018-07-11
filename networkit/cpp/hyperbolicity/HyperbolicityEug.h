/*
 * Hyperbolicity.h
 *
 *  Created on: 07.05.2018
 *      Author: Eugenio Angriman
 */

#ifndef HYPERBOLICITY_H_
#define HYPERBOLICITY_H_

#include "../algebraic/CSRMatrix.h"
#include "../auxiliary/PrioQueue.h"
#include "../base/Algorithm.h"
#include "../centrality/GroupDegree.h"
#include "../distance/Bisenius/PrunedAPSP.h"
#include "../graph/Graph.h"

#include <math.h>

namespace NetworKit {

/**
 * @ingroup distance
 */
class Hyperbolicity : public Algorithm {

public:
  Hyperbolicity(const Graph &G);
  void run() override;
  double getHyperbolicity();
  std::vector<node> getQuadruple();

protected:
  Graph G;
  const count n;
  // 'Compressed' keys. (0,1) = 0, (0,2) = 1, ... (1,2) = n
  Aux::PrioQueue<int64_t, uint64_t> pairs;
  PrunedAPSP prunedAPSP;
  std::vector<count> minPivotDist;
  std::vector<count> maxPivotDist;
  std::vector<count> degreeInCmp;
  std::vector<count> ecc;
  std::vector<count> distCentral;
  std::vector<node> nearestPivot;
  std::vector<bool> isInCmp;
  std::vector<bool> seenBool;
  std::vector<node> seen;
  std::vector<bool> valuable;
  std::vector<bool> isPivot;
  std::vector<node> visited;
  std::vector<bool> acc;
  std::vector<count> sourceDist;
  std::vector<std::vector<node>> mate;
  // std::vector<uint64_t> nonFarApartPairs;
  CSRMatrix nonFarApartPairs;
  std::vector<node> sortedComponent;
  CSRMatrix distance;
  uint64_t delta, deltaPlusOne, condAcc1, condAcc2, condVal, maxHyp;
  node mostCentral;
  count dist_xy, dist_xv, dist_yv, dist_xw, dist_yw, dist_vw, nPivots, ecc_v,
      maxDeg;
  node solX, solY, solV, solW;

  int64_t s2, s3, diff;

  void init();
  void initHypComponent(const std::vector<node> &component);
  count getDistanceUnknown(node u, node v);
  void computeHypOfComponent(const std::vector<node> &component);
  void bfsFrom(node source, std::vector<bool> &isInCmp);
  void storePairs(const std::vector<node> &component,
                  const std::vector<bool> &isInCmp);
  uint64_t pairIndex(node u, node v);
  std::pair<node, node> pairFromIndex(uint64_t idx);
  count numberOfPivots(count size);
  int64_t distUpperBound(node u, node v);
  count eccUpperBound(node u);
  count eccLowerBound(node u);
  void updateDelta(node x, node y, node v);
  void checkHasRun();
  void setNonFarApart(node u, node v);
};

inline void Hyperbolicity::setNonFarApart(node u, node v) {
  node _u = std::min(u, v);
  node _v = std::max(u, v);
  nonFarApartPairs.setValue(_u, _v, 1);
}

inline count Hyperbolicity::numberOfPivots(count size) {
  return (count)(size * 0.1 + 0.5); //(log2(size) + 0.5);
}

inline uint64_t Hyperbolicity::pairIndex(node u, node v) {
  assert(u != v);
  node from = std::min(u, v);
  node to = std::max(u, v) - 1;
  return (n - 1) * from + to;
}

inline std::pair<node, node> Hyperbolicity::pairFromIndex(uint64_t idx) {
  node from = idx / (n - 1);
  node to = idx - from * (n - 1) + 1;
  return std::make_pair(from, to);
}

inline count Hyperbolicity::getDistanceUnknown(node u, node v) {
  assert(u != v);

  if (G.hasEdge(u, v)) {
    return 1;
  }

  node _u = std::min(u, v);
  node _v = std::max(u, v);
  count d = distance(_u, _v);

  if (d == 0) {
    d = prunedAPSP.getDistance(u, v);
    distance.setValue(_u, _v, d);
  }

  return d;
}

inline int64_t Hyperbolicity::distUpperBound(node u, node v) {
  if (u == v) {
    return 0;
  }
  if (G.hasEdge(u, v)) {
    return 1;
  }

  return std::min(maxPivotDist[u] + minPivotDist[v],
                  minPivotDist[u] + maxPivotDist[v]);
}

inline count Hyperbolicity::eccUpperBound(node u) {
  return isPivot[u] ? ecc[u] : minPivotDist[u] + ecc[nearestPivot[u]];
}

inline count Hyperbolicity::eccLowerBound(node u) {
  return isPivot[u] ? ecc[u] : maxPivotDist[u];
}

inline double Hyperbolicity::getHyperbolicity() {
  checkHasRun();
  return ((double)maxHyp / 2.f);
}

inline std::vector<node> Hyperbolicity::getQuadruple() {
  checkHasRun();

  return {solX, solY, solV, solW};
}

inline void Hyperbolicity::checkHasRun() {
  if (!hasRun) {
    throw std::runtime_error("Run method has not been called.");
  }
}
} // namespace NetworKit

#endif /* ifndef HYPERBOLICITY_H_ */
