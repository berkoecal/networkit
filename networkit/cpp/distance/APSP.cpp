/*
 * APSP.cpp
 *
 *  Created on: 07.07.2015
 *      Author: Arie Slobbe
 */

#include "APSP.h"
#include "../auxiliary/Log.h"
#include "Dijkstra.h"
#include "BFS.h"

namespace NetworKit {

APSP::APSP(const Graph& G) : Algorithm(), G(G) {}

void APSP::run() {
	std::vector<edgeweight> distanceVector(G.upperNodeIdBound(), 0.0);
	distances.resize(G.upperNodeIdBound(), distanceVector);
	eccentricity.resize(G.upperNodeIdBound(),0.0);
	if (G.isWeighted()) {
		G.parallelForNodes([&](node u){
			Dijkstra dijk(G, u);
			dijk.run();
			distances[u] = dijk.getDistances();
		});
	} else {
		G.parallelForNodes([&](node u){
			BFS bfs(G, u, true, true);
			bfs.run();
			distances[u] = bfs.getDistances();
			
			auto sorted_nodes = bfs.getNodesSortedByDistance();
			auto const& eccentricity_node = sorted_nodes.back();
			eccentricity[u] = distances[u][eccentricity_node];
		});
	}
	hasRun = true;
}

std::string NetworKit::APSP::toString() const {
	return "All-Pairs Shortest Path Algorithm";
}

} /* namespace NetworKit */
