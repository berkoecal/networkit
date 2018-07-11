/*
 * CustomizedBFS.h
 *
 *  Created on: 18.02.2018
 *      Author: Berk Öcal
 */
#ifndef CUSTOMIZEDBFS_H
#define CUSTOMIZEDBFS_H

#include "SymmetricMatrix.h"
#include "../graph/Graph.h"
#include "../auxiliary/PrioQueue.h"
#include <queue>


namespace NetworKit{

enum class CentralNodeMethod{
		NATURAL,
		DEGREE,
		TOPCLOSENESS,
		RANDOM,
		GROUPDEGREE
};

/*
 * Class for computing upper bounds for distances of node pairs.
 * This class implements the notions of the paper
 * "Fast Shortest Path Distance Estimation in Large Networks" (Potamias et al.)
 */

class CustomizedAPSP{
public:
	/*
	 * The constructor requires a list of nodes it applies BFS for.
	 * If the list of @param centralNodes is a proper subset of all nodes,
	 *  the algorithm computes upper bounds for pairs not considered in the BFS searchs that have been run.
	 * Else, it computes exact distances.
	 */
	CustomizedAPSP(const Graph& G, CentralNodeMethod method = CentralNodeMethod::NATURAL);
	CustomizedAPSP(const Graph& G, count numberOfLandmarks, count maxComparisons, CentralNodeMethod method = CentralNodeMethod::DEGREE);
	virtual void run();
	SymMatrix<bool, node>& getRelevantPairs();
	SymMatrix<edgeweight, node> const& getExactDistances();
	SymMatrix<edgeweight, node>& getDistances();
	edgeweight getDistance(node u, node v);
	edgeweight distUpperBoundByComp(node u, node v);
	edgeweight distLowerBoundByComp(node u, node v);
	edgeweight distUpperBound(node u, node v);
	edgeweight distUpperBoundAdvanced(node u, node v);
	edgeweight distLowerBoundAdvanced(node u, node v);
	std::vector<edgeweight>& getEccentricityVec();
	std::vector<edgeweight>& getFarnessVec();
	edgeweight getMinFarness();
	edgeweight eccLowerBound(node u);
	node ithPivot(count i);

private:
	const Graph& G;
	SymMatrix<bool, node> relevantPairs;
	SymMatrix<edgeweight, node> distances;
	std::vector<edgeweight> ecc;
	std::vector<bool>  isPivot;
	std::vector<edgeweight> minPivotDist;
	std::vector<edgeweight> maxPivotDist;
	std::vector<node>  nearestPivot;
	std::vector<node>  farestPivot;
	std::vector<edgeweight> farnessToPivots;
	edgeweight minFarness;
	std::vector<node> pivots;
	const count numberOfLandmarks;
	const count maxComparisons;
	CentralNodeMethod method;
	void setBoundedDistances();
	void computeCentralNodes();
	bool hasRun;
};

inline SymMatrix<bool,node>& CustomizedAPSP::getRelevantPairs(){return relevantPairs;}
inline SymMatrix<edgeweight, node>& CustomizedAPSP::getDistances(){return distances;}
inline edgeweight CustomizedAPSP::getDistance(node u, node v){return distances.element(u,v);}
inline std::vector<edgeweight>& CustomizedAPSP::getEccentricityVec(){return ecc;}
inline std::vector<edgeweight>& CustomizedAPSP::getFarnessVec(){return farnessToPivots;}
inline edgeweight CustomizedAPSP::getMinFarness(){return minFarness;}
inline edgeweight CustomizedAPSP::eccLowerBound(node u){return isPivot[u] ? ecc[u] : maxPivotDist[u];}
inline node CustomizedAPSP::ithPivot(count i) {return pivots[i];}

}
#endif /*CUSTOMIZEDBFS_H */
