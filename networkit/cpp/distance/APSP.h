/*
 * APSP.h
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#ifndef APSP_H_
#define APSP_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"

namespace NetworKit {

/**
 * @ingroup distance
 * Class for all-pair shortest path algorithm.
 */
class APSP: public Algorithm {

public:

	/**
	 * Creates the APSP class for @a G.
	 *
	 * @param G The graph.
	 */
	APSP(const Graph& G);

	virtual ~APSP() = default;

	/** Computes the shortest paths from each node to all other nodes. */
	void run() override;

	/**
	* @return string representation of algorithm and parameters.
	*/
	virtual std::string toString() const override;

	/**
	 * Returns a vector of weighted distances between node pairs.
 	 *
 	 * @return The shortest-path distances from each node to any other node in the graph.
	 */
	std::vector<std::vector<edgeweight> > getDistances() const {if (!hasRun) throw std::runtime_error("Call run method first"); return distances;}


	/**
	 * Returns the distance from u to v or infinity if u and v are not connected.
	 *
	 */
	edgeweight getDistance(node u, node v) const {if (!hasRun) throw std::runtime_error("Call run method first"); return distances[u][v];}
	
	/**
	 * Returns eccentricity of a node u (= max_{v node of G} d(u,v)).
	 * It is calculated during the BFS iteration.
	 *
	 */
	
	edgeweight getEccentricity(node u) const {if (!hasRun) throw std::runtime_error("Call run method first"); return eccentricity[u];}
	
	/**
	* @return True if algorithm can run multi-threaded.
	*/
	virtual bool isParallel() const override { return true; }


protected:

	const Graph& G;
	std::vector<std::vector<edgeweight> > distances;
	std::vector<edgeweight> eccentricity;
};

} /* namespace NetworKit */

#endif /* APSP_H_ */
