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
#include <queue>

namespace NetworKit{

/*
 * Class for computing upper bounds for distances of node pairs.
 * This class implements the notions of the paper "Fast Shortest Path Distance Estimation in Large Networks" (Potamias et al.)
 */

class CustomizedBFS{
public:
	/*
	 * The constructor requires a list of nodes it applies BFS for.
	 * If the list of @param centralNodes is a proper subset of all nodes,
	 *  the algorithm computes upper bounds as in the paper.
	 * Else, it computes exact distances.
	 */
	CustomizedBFS(const Graph& G);
	CustomizedBFS(const Graph& G, const std::vector<node>& centralNodes);
	virtual void run();
	SymMatrix<bool, node> getRelevantPairs();
	SymMatrix<edgeweight, node> getExactDistances();
	SymMatrix<edgeweight, node> getBoundedDistances();
	edgeweight getBoundedDistance(node u, node v);
	edgeweight getEccentricity(node v);

private:
	const Graph& graph;
	SymMatrix<bool, node> relevantPairs;
	SymMatrix<edgeweight, node> boundedDistances;
	std::vector<edgeweight> eccentricity;
	std::vector<node> const nodesToConsider;
	edgeweight computeBoundedDistance(node u, node v);
	void setBoundedDistances();

};

}
#endif /*CUSTOMIZEDBFS_H */
