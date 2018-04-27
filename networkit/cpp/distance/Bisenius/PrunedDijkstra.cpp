#include "PrunedDijkstra.h"
#include "../../auxiliary/PrioQueue.h"

namespace NetworKit {

PrunedDijkstra::PrunedDijkstra(const Graph &G, const std::vector<node> &nodeToRank, std::vector<LandmarkLabel> &forwardLabels, std::vector<LandmarkLabel> &backwardLabels) : PrunedSearch(G, nodeToRank, forwardLabels, backwardLabels) {}

void PrunedDijkstra::run(node source, bool searchForward) {
	std::vector<node> visited_nodes;
	std::vector<LandmarkLabel> &labels = searchForward ? forwardLabels : backwardLabels;
	setLabels(searchForward);
	
	node rank = nodeToRank[source];
	
	cacheDistances((*targetLabels)[rank]);
	
	distances[source] = 0;
	visited_nodes.push_back(source);
	Aux::PrioQueue<edgeweight, node> pq(distances.size());
	pq.insert(0, source);
	auto relax([&](node u, node v, edgeweight weight) {
		if (distances[v] > distances[u] + weight) {
			distances[v] = distances[u] + weight;
			parents[v] = nodeToRank[u];
			visited_nodes.push_back(v);
			pq.decreaseKey(distances[v], v);
		}
	});
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	
	while (pq.size() > 0) {
		auto minimum = pq.extractMin();
		
		node currentNode = minimum.second;
		node currentRank = nodeToRank[currentNode];
		if(query(currentRank, rank) <= distances[currentNode]) {
			continue;
		}
		labels[currentRank].nodes.push_back(rank);
		labels[currentRank].distances.push_back(distances[currentNode]);
		labels[currentRank].parents.push_back(parents[currentNode]);
		if(searchForward)
			G.forEdgesOf(currentNode, relax);
		else
			G.forInEdgesOf(currentNode, relax);
	}


	for(auto v: visited_nodes) {
		distances[v] = infDist;
	}
	cleanupCache();
}

} /* namespace NetworKit */
