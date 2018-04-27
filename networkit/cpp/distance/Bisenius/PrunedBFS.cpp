#include "PrunedBFS.h"

namespace NetworKit {

PrunedBFS::PrunedBFS(const Graph &G, const std::vector<node> &nodeToRank, std::vector<LandmarkLabel> &forwardLabels, std::vector<LandmarkLabel> &backwardLabels) : PrunedSearch(G, nodeToRank, forwardLabels, backwardLabels) {
	count z = G.upperNodeIdBound();
	visited.resize(z);
	parents.resize(z);
}

void PrunedBFS::run(node v, bool searchForward) {
	node rank = nodeToRank[v];
	count c = 0;
	std::vector<node> visited_nodes;
	std::vector<LandmarkLabel> &labels = searchForward ? forwardLabels : backwardLabels;
	setLabels(searchForward);
	
	cacheDistances((*targetLabels)[rank]);
	
	std::queue<node> q;
	q.push(v);
	visited_nodes.push_back(v);
	visited[v] = true;
	distances[v] = 0;
	parents[v] = G.upperNodeIdBound() + 1;

	while (! q.empty()) {
		node u = q.front();
		q.pop();
		node neighborRank = nodeToRank[u];
		c++;
		if(query(neighborRank, rank) <= distances[u]) {
			continue;
		}
		
		labels[neighborRank].nodes.push_back(rank);
		labels[neighborRank].distances.push_back(distances[u]);
		labels[neighborRank].parents.push_back(parents[u]);
		
		
		auto visit ([&](node w) {
				if (!visited[w]) {
					q.push(w);
					visited_nodes.push_back(w);
					visited[w] = true;
					distances[w] = distances[u] + 1;
					parents[w] = nodeToRank[u];
				}
			});
		if(searchForward) {
			G.forNeighborsOf(u, visit);
		} else {
			G.forInNeighborsOf(u, visit);
		}
	}

	for(auto u : visited_nodes) {
		visited[u] = false;
		parents[u] = G.upperNodeIdBound() + 1;
	}
	
	cleanupCache();
	consideredPairs = c;
}

} /* namespace NetworKit */
