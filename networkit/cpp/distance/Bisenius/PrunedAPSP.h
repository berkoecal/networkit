/*
 * AkibaPrunedLandmarks.h
 */

#ifndef AKIBAPRUNEDLANDMARKS_H_
#define AKIBAPRUNEDLANDMARKS_H_

#include "LandmarkLabel.h"
#include "../../distance/APSP.h"
#include "../../graph/Graph.h"


namespace NetworKit {
	
	typedef enum PrunedAPSPNodeOrder {
		NATURAL,
		RANDOM,
		DEGREE,
		CLOSENESS
	} PrunedAPSPNodeOrder;

/**
 * @ingroup graph
 *
 * After preprocessing, this class can compute the exact distance
 * between two nodes very quickly.
 * 
 * See https://dl.acm.org/citation.cfm?id=2465315
 */
class PrunedAPSP: public APSP {

public:

	/**
	 * @param G The graph.
	 * @param nodeOrder The order in which the nodes should be considered during preprocessing. Ordering by degree centrality gave the best results in experiments 
	 */
	PrunedAPSP(const Graph& G, PrunedAPSPNodeOrder nodeOrder = PrunedAPSPNodeOrder::DEGREE);

	/**
	* Starts the preprocessing. This must be called before calling any of the other methods.
	*/
	virtual void run();

	/**
	 * Computes the distance between two nodes in the graph.
	 * 
	 * @param u Source node
	 * @param v Target node
	 * @return distance between the nodes 
	 */
	virtual double getDistance(node u, node v);
	
	/**
	 * Computes a path between two nodes.
	 * 
	 * @param u Source node
	 * @param v Target node
	 * @return The path from u to v as a vector of nodes
	 */
	virtual std::vector<node> getPath(node u, node v);

	/**
	 * Returns a nxn dimensional vector with distances between all node pairs
	 * @return A two-dimensional vector including all distances. Result[u][v] contains the distance from node u to node v.
	 */
	virtual std::vector<std::vector<edgeweight>> getDistances();

protected:
	/**
	 * The Query. Cuts two labels and returns the distance
	 * @param sourceLabel The label of the source node
	 * @param targetLabel The label of the target node
	 * @return The distance between the two nodes
	 */
	virtual double query(LandmarkLabel sourceLabel, LandmarkLabel targetLabel);

	/**
	 * Computes the hub node on the shortest path between two nodes.
	 * The hub node is a node whose distance is known from both the
	 * source and the target node of the query and where the
	 * distance is minimal for all possible hub nodes.
	 * 
	 * @param sourceLabel The label of the source node
	 * @param targetLabel The label of the target node
	 * @return The indices of the hubs in the source and target label
	 */
	virtual std::pair<size_t, size_t> getHubIndizes(LandmarkLabel sourceLabel, LandmarkLabel targetLabel);

	/**
	 * Sorts the nodes in the graph based on the supplied strategy
	 */
	virtual void orderNodes();
	
	std::vector<LandmarkLabel> forwardLabels;
	std::vector<LandmarkLabel> backwardLabels;
	std::vector<LandmarkLabel> *sourceLabels;
	std::vector<LandmarkLabel> *targetLabels;
    std::vector<node> nodeToRank;
    std::vector<node> rankToNode;
	
	PrunedAPSPNodeOrder nodeOrder;
};

} /* namespace NetworKit */
#endif /* AKIBAPRUNEDLANDMARKS_H_ */
