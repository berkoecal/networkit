/*
 * PrunedSearch.h
 */

#ifndef PRUNEDSEARCH_H_
#define PRUNEDSEARCH_H_

#include "../../graph/Graph.h"
#include "LandmarkLabel.h"

namespace NetworKit {

/**
* This class provides an interface to a pruned version of BFS and Dijkstras algorithm used by the PrunedAPSP algorithm
*/
class PrunedSearch {
protected:
	const Graph& G;
    const std::vector<node>& nodeToRank;
    std::vector<LandmarkLabel> &forwardLabels;
	std::vector<LandmarkLabel> &backwardLabels;
	std::vector<LandmarkLabel> *sourceLabels;
	std::vector<LandmarkLabel> *targetLabels;

public:
	/**
	* @param G The graph
	* @param noteToRank A vector maping every node id to it's ordered rank used in the PrunedAPSP algorithm
	* @param forwardlabel The labels used for forward search
	* @ param backwardLabels The labels used for backward search
	*/
    PrunedSearch(const Graph &G, const std::vector<node> &nodeToRank, std::vector<LandmarkLabel> &forwardLabels, std::vector<LandmarkLabel> &backwardLabels);

    /** Default destructor */
    virtual ~PrunedSearch() = default;

    /**
    * Runs the pruned search
    * @param u The node to start the search from
    * @param searchForward Specifies wether the search should use the edges in the graph in the direction they point or the reverse direction
    */
    virtual void run(node u, bool searchForward) = 0;

    count numOfConsideredPairs();
protected:
	/**
	* Perfom a query tuned for pruned search
	* @param source The label of the source node
	* @param target the label of the target node
	* @return the distance from source to target
	*/
	double query(LandmarkLabel sourceLabel, LandmarkLabel targetLabel);

	/**
	* Perfom a query tuned for pruned search
	* @param source The rank of the source node
	* @param target the rank of the target node
	* @return the distance from source to target
	*/
	double query(node source, node target);

	/**
	* Cache the distances in the specified label to speed up the query while performing the pruned search
	* @param sourceNodeLabel The label to cache distances from
	*/
	void cacheDistances(LandmarkLabel sourceNodeLabel);

	/**
	* Clean up the cached data
	*/
	void cleanupCache();

	/**
	* Sets up the labels for the pruned search
	* @param searchForward Wether the pruned search uses the edges in the graph in the direction they point or the reverse direction
	*/
	void setLabels(bool searchForward);
	std::vector<edgeweight> distances;
	std::vector<node> parents;
	std::vector<edgeweight> cachedDistances;
	std::vector<node> cachedNodes;
	count consideredPairs;

};

} /* namespace NetworKit */
#endif /* PRUNEDSEARCH_H_ */
