/*
 * PrunedDijkstra.h
 */

#ifndef PRUNEDDIJKSTRA_H_
#define PRUNEDDIJKSTRA_H_

#include "../../graph/Graph.h"
#include "LandmarkLabel.h"
#include "PrunedSearch.h" 

namespace NetworKit {

/**
* This class provides an implementation of the PrunedSearch interface for Dijkstras algorithm
*/
class PrunedDijkstra: public NetworKit::PrunedSearch {
public:
	/**
	* @param G The graph
	* @param noteToRank A vector maping every node id to it's ordered rank used in the PrunedAPSP algorithm
	* @param forwardlabel The labels used for forward search
	* @ param backwardLabels The labels used for backward search
	*/
	PrunedDijkstra(const Graph &G, const std::vector<node> &nodeToRank, std::vector<LandmarkLabel> &forwardLabels, std::vector<LandmarkLabel> &backwardLabels);

	/**
    * Runs the pruned search
    * @param u The node to start the search from
    * @param searchForward Specifies wether the search should use the edges in the graph in the direction they point or the reverse direction
    */
	virtual void run(node u, bool searchForward);
};

} /* namespace NetworKit */
#endif /* PRUNEDDIJKSTRA_H_ */
