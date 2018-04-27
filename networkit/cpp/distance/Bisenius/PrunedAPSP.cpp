#include <memory>
#include "PrunedAPSP.h"
#include "PrunedSearch.h"
#include "PrunedBFS.h"
#include "PrunedDijkstra.h"
#include "../../centrality/ApproxCloseness.h"
#include "../../auxiliary/Timer.h"
namespace NetworKit {

PrunedAPSP::PrunedAPSP(const Graph& G, PrunedAPSPNodeOrder nodeOrder) : APSP(G), nodeOrder(nodeOrder) {}

void PrunedAPSP::run() {
	auto maxNodeId = G.upperNodeIdBound();
	forwardLabels =  std::vector<LandmarkLabel>(maxNodeId);
	sourceLabels = &forwardLabels;
	targetLabels = &forwardLabels;
	
	// we need two sets of labels for directed graphs
	if(G.isDirected()) {
		backwardLabels =  std::vector<LandmarkLabel>(maxNodeId);
		targetLabels = &backwardLabels;
	}
	
    nodeToRank = std::vector<node>(maxNodeId);
    rankToNode = std::vector<node>(maxNodeId);
	
	orderNodes();

	// initialize pruned search algorithm: BFS for unweighted graphs, Dijkstra
	// for weighted ones
    std::unique_ptr<PrunedSearch> prunedSearch;
    if(G.isWeighted()) {
		prunedSearch.reset( new PrunedDijkstra(G, nodeToRank, *sourceLabels, *targetLabels));
	} else {
    	prunedSearch.reset( new PrunedBFS(G, nodeToRank, *sourceLabels, *targetLabels));
	}
	
    count sumCounter=0;
	// run pruned search from each node
	for (node v : rankToNode) {
		prunedSearch->run(v, true);
		sumCounter+=prunedSearch->numOfConsideredPairs();
		if (G.isDirected()) {
			prunedSearch->run(v, false);
		}
	}

	INFO("Number of visited pairs in Preprocessing of Pruned: ", sumCounter, " with percentage of ", (sumCounter/static_cast<double>((G.numberOfNodes()*G.numberOfNodes()/2)))*100,"%");
}

double PrunedAPSP::getDistance(node u, node v) {
	return query((*targetLabels)[nodeToRank[u]], (*sourceLabels)[nodeToRank[v]]);
}

double PrunedAPSP::query(LandmarkLabel sourceLabel, LandmarkLabel targetLabel) {
	size_t indexSource = 0;
	size_t indexTarget = 0;

	double distance = std::numeric_limits<double>::max();

	while(indexSource < sourceLabel.nodes.size() && indexTarget < targetLabel.nodes.size()) {
		int64_t comp = sourceLabel.nodes[indexSource] - targetLabel.nodes[indexTarget];
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
}

std::pair<size_t, size_t> PrunedAPSP::getHubIndizes(LandmarkLabel sourceLabel, LandmarkLabel targetLabel) {
	size_t indexSource = 0;
	size_t indexTarget = 0;
	
	std::pair<size_t, size_t> indizes(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
	
	double distance = std::numeric_limits<double>::max();

	/* Since the nodes in the labels are always ordered by node ID, we can
	 * use a rather simple algorithm to determine all entries which
	 * are present in both supplied labels. We keep a counter for
	 * both the source and the target label. In each iteration, we compare
	 * the node IDs which both counters point to. If the node ID in the
	 * source label is lower, we increment the counter for the source label.
	 * If it is larger, we increment in the counter for the target label.
	 * If they are equal, we have found a common node and compute the
	 * resulting distance and check if it is a new global minimum.
	 */
	while(indexSource < sourceLabel.nodes.size() && indexTarget < targetLabel.nodes.size()) {
		int64_t comp = sourceLabel.nodes[indexSource] - targetLabel.nodes[indexTarget];
		if(comp < 0) {
			++indexSource;
		} else if(comp > 0) {
			++indexTarget;
		} else {
			auto temp_distance = sourceLabel.distances[indexSource] + targetLabel.distances[indexTarget];
			if(temp_distance < distance) {
				distance = temp_distance;
				indizes.first = indexSource;
				indizes.second = indexTarget;
			}
			++indexTarget;
			++indexSource;
		}
	}
	
	return indizes;
}

std::vector<node> PrunedAPSP::getPath(node u, node v) {
	
	LandmarkLabel uLabels = (*targetLabels)[nodeToRank[u]];
	LandmarkLabel vLabels = (*sourceLabels)[nodeToRank[v]];
	
	size_t infIndex = std::numeric_limits<size_t>::max();
	
	std::pair<size_t, size_t> hubIndizes = getHubIndizes(uLabels, vLabels);
	
	// If there is no path between u and v, return an empty vector
	if (hubIndizes.first == infIndex || hubIndizes.second == infIndex) {
		return std::vector<node>();
	}
	
	node hubRank = uLabels.nodes[hubIndizes.first];
	
	std::vector<node> path;
	
	// Walk to the hub node from u
	node currentRank = nodeToRank[u];
	size_t currentIndex = hubIndizes.first;
	LandmarkLabel currentLabels = (*targetLabels)[currentRank];
	
	path.push_back(rankToNode[currentRank]);

	while (currentRank != hubRank && currentRank < G.upperNodeIdBound()) {
		// new rank is rank of parent
		currentRank = currentLabels.parents[currentIndex];
		// new index is index of hub node in parent label
		currentLabels = (*targetLabels)[currentRank];
		currentIndex = std::distance(currentLabels.nodes.begin(), std::lower_bound(currentLabels.nodes.begin(), currentLabels.nodes.end(), hubRank));
		path.push_back(rankToNode[currentRank]);
	}
	
	currentRank = nodeToRank[v];
	currentIndex = hubIndizes.second;
	currentLabels = (*sourceLabels)[currentRank];
	
	std::vector<node> reversePath;
	reversePath.push_back(rankToNode[currentRank]);
	while (currentRank != hubRank && currentRank < G.upperNodeIdBound()) {
		// new rank is rank of parent
		currentRank = currentLabels.parents[currentIndex];
		
		// new index is index of hub node in parent label
		currentLabels = (*sourceLabels)[currentRank];
		currentIndex = std::distance(currentLabels.nodes.begin(), std::lower_bound(currentLabels.nodes.begin(), currentLabels.nodes.end(), hubRank));
		reversePath.push_back(rankToNode[currentRank]);
	}
	
	//INFO("Reverse path to hub: ", reversePath);
	
	// remove last element from reverse path to avoid having the hub element twice
	reversePath.pop_back();
	
	std::reverse(reversePath.begin(), reversePath.end());
	
	path.insert(path.end(), reversePath.begin(), reversePath.end());
	
	return path;
	
}

std::vector<std::vector<edgeweight>> PrunedAPSP::getDistances() {
	if (!hasRun) {
		throw new std::runtime_error("Call run() first before accessing distances.");
	}
	
	bool isDirected = G.isDirected();
	
	G.forNodePairs([&](node u, node v) {
		distances[u][v] = getDistance(u, v);
		
		if (isDirected) {
			distances[v][u] = getDistance(v, u);
		} else {
			distances[v][u] = distances[u][v];
		}
	});
	
	return distances;
	
}

void PrunedAPSP::orderNodes() {
	std::vector<node> nodes = G.nodes();
	
	switch(nodeOrder) {
		case PrunedAPSPNodeOrder::DEGREE:
			std::sort(nodes.begin(), nodes.end(), [&](node a, node b) {
				return G.degree(a) > G.degree(b);
			});
			break;
		case PrunedAPSPNodeOrder::RANDOM:
			std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
			break;
		case PrunedAPSPNodeOrder::CLOSENESS:
			{	
				count nSamples = std::log(G.numberOfNodes()) / 0.01;
				ApproxCloseness closeness(G, nSamples);
				closeness.run();
				std::vector<std::pair<node, double>> ranking = closeness.ranking();
				nodes.clear();
				for (std::pair<node, double> pair : ranking) {
					nodes.push_back(pair.first);
				}
			}
			break;
		case NATURAL:  
		default:
			break;
	}
	
	for (size_t i = 0; i < nodes.size(); i++) {
		node nodeId = nodes[i];
		nodeToRank[nodeId] = i;
		rankToNode[i] = nodeId;
	}
}

} /* namespace NetworKit */
