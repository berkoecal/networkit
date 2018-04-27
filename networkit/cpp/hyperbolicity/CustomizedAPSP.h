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
		RANDOM
};

//class NodeTupleWithDist
//{
//public:
//	NodeTupleWithDist(node index_a, node index_b, edgeweight dist_) :
//		indexA(index_a),
//		indexB(index_b),
//		dist(dist_){}
//
//	node index_a() const{return indexA;}
//	node index_b() const{return indexB;}
//	edgeweight distance() const {return dist;}
//
//private:
//	node indexA;
//	node indexB;
//	edgeweight dist;
//};
//
//struct CustomCompare{
//	bool operator() (NodeTupleWithDist const& firstPair, NodeTupleWithDist const& secondPair) const{
//		return firstPair.distance() >= secondPair.distance();
//	}
//};

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
	CustomizedAPSP(const Graph& G, count numberOfLandmarks, CentralNodeMethod method = CentralNodeMethod::DEGREE);
	virtual void run();
	SymMatrix<bool, node> getRelevantPairs() const;
	//std::set<NodeTupleWithDist, CustomCompare> const& getRelPairs() const;
	SymMatrix<edgeweight, node> const& getExactDistances() const;
	SymMatrix<edgeweight, node> getDistances() const;
	edgeweight getDistance(node u, node v);
	edgeweight getEccentricity(node v) const;
	//count numOfRelPairs;

private:
	const Graph& graph;
	SymMatrix<bool, node> relevantPairs;
	//std::set<NodeTupleWithDist, CustomCompare> relPairs;
	//Aux::PrioQueue<>
	SymMatrix<edgeweight, node> distances;
	std::vector<edgeweight> eccentricity;
	std::vector<node> nodesToConsider;
	const count numberOfLandmarks;
	CentralNodeMethod method;
	edgeweight computeBoundedDistance(node u, node v);
	void setBoundedDistances();
	void computeCentralNodes();

};

}
#endif /*CUSTOMIZEDBFS_H */
