#include "PrunedSearch.h"

#define CACHED_QUERY

namespace NetworKit {

PrunedSearch::PrunedSearch(const Graph &G, const std::vector<node> &nodeToRank, std::vector<LandmarkLabel> &forwardLabels, std::vector<LandmarkLabel> &backwardLabels) : G(G), nodeToRank(nodeToRank), forwardLabels(forwardLabels), backwardLabels(backwardLabels) {
	count z = G.upperNodeIdBound();
	edgeweight infDist = std::numeric_limits<edgeweight>::max();
	distances.resize(z, infDist);
	parents.resize(z);
	cachedDistances.resize(z, infDist);
}

double PrunedSearch::query(LandmarkLabel sourceLabel, LandmarkLabel targetLabel) {
	
#ifdef CACHED_QUERY
	
	double infDist = std::numeric_limits<edgeweight>::max();
	
	double distance = infDist;
	
	for (size_t indexTarget = 0; indexTarget < sourceLabel.nodes.size(); indexTarget++) {
		
		double cachedDist = cachedDistances[sourceLabel.nodes[indexTarget]];
		double sourceDist = sourceLabel.distances[indexTarget];
		
		if (cachedDist == infDist || sourceDist == infDist) {
			continue;
		}
		
		double temp_distance = cachedDist + sourceDist;
		
		if(temp_distance < distance) {
			distance = temp_distance;
		}
	}
	
	return distance;
#else
	size_t indexSource = 0;
	size_t indexTarget = 0;
	double distance = std::numeric_limits<double>::max();
	while(indexSource < sourceLabel.nodes.size() && indexTarget < targetLabel.nodes.size()) {
		auto comp = sourceLabel.nodes[indexSource] - targetLabel.nodes[indexTarget];
		if(comp < 0) {
			++indexSource;
		} else if(comp > 0) {
			++indexTarget;
		} else {
			auto temp_distance = sourceLabel.distances[indexSource] + targetLabel.distances[indexTarget];
			if(temp_distance < distance) {
				distance = temp_distance;
			}
			++indexTarget;
			++indexSource;
		}
	}
	return distance;
#endif
}

double PrunedSearch::query(node sourceRank, node targetRank) {
	return query((*sourceLabels)[sourceRank], (*targetLabels)[targetRank]);
}

void PrunedSearch::cacheDistances(LandmarkLabel sourceNodeLabel) {
#ifdef CACHED_QUERY
	for (count i = 0; i < sourceNodeLabel.nodes.size(); i++) {
		node nodeId = sourceNodeLabel.nodes[i];
		cachedNodes.push_back(nodeId);
		cachedDistances[nodeId] = sourceNodeLabel.distances[i];
	}
#endif
}

void PrunedSearch::cleanupCache() {
#ifdef CACHED_QUERY
	for (node nodeId : cachedNodes) {
		cachedDistances[nodeId] = std::numeric_limits<edgeweight>::max();
	}
	cachedNodes.clear();
#endif
}

void PrunedSearch::setLabels(bool searchForward) {
	sourceLabels = searchForward ? &forwardLabels : &backwardLabels;
	targetLabels = searchForward ? &backwardLabels : &forwardLabels;
}

count PrunedSearch::numOfConsideredPairs() {
	count temp = consideredPairs;
	consideredPairs = 0;
	return temp;
}



} /* namespace NetworKit */
