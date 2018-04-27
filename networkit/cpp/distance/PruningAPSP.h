/*
 * BFSPruned.h
 *       
 *   Created on: 01.12.2017
 *       Author: Berk Öcal
 */

#ifndef PruningAPSP_H
#define PruningAPSP_H

#include "../graph/Graph.h"
#include "../base/Algorithm.h"
namespace NetworKit{

class NodeWithDistance{
public:
	NodeWithDistance(node u, edgeweight distance): u(u), distance(distance){}
	node getNode() const{ return u;}
	edgeweight getDistance() const {return distance;}
private:
	node u;
	edgeweight distance;
};


class PruningAPSP : public Algorithm
{
public:
	/**
	* Constructs the PruningAPSP class for @a G and source node @a source.
	*
	* @param G The graph
	*/
	PruningAPSP(const Graph& G);
	
	/**
	* Breadth-first search from @a source.
	* @return Vector of unweighted distances from node @a source, i.e. the
	* length (number of edges) of the shortest path from @a source to any other node.
	*/
	virtual void run();
	virtual void runParallel();
	
	edgeweight getDistance (node v, node u) { if (!hasRun)  throw std::runtime_error("Call run method first"); return mergingComputation(v,u);};
	

	
private:
	const Graph& G;
	std::vector< std::vector<NodeWithDistance>> L;
	std::vector<std::vector<edgeweight>> distances;
	//std::vector<std::vector<node>> alreadyKnownDistances;
	void prunedBFS(node v);
	
	/**
	 * Internal querying step (i.e. preprocessing) of Akiba et al.
	 *
	 * @param T vector of distances
	 * @param v node
	 */
	edgeweight mergingComputation(std::vector<edgeweight> const& T, node u) const;
	
	/**
	 * Querying step of Akiba et al.
	 *
	 * @param u node
	 * @param v node
	 */
	edgeweight mergingComputation(node v, node u) const;
	std::vector<edgeweight> computeDistanceVector(node v) const;	
};

}
#endif // PruningAPSP_H
