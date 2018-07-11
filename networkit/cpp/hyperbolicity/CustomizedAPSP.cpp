/*
 * CustomizedBFS.cpp
 *
 *  Created on: 21.02.2018
 *      Author: Berk Öcal
 */
#define _BSD_SOURCE

#include <sys/time.h>
#include "CustomizedAPSP.h"
#include "../centrality/DegreeCentrality.h"
#include "../centrality/TopCloseness.h"
#include "../centrality/GroupDegree.h"

namespace NetworKit{

CustomizedAPSP::CustomizedAPSP(const Graph& G, CentralNodeMethod method): G(G),
		relevantPairs(G.upperNodeIdBound(), true),
		distances(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
		ecc(G.upperNodeIdBound()),
		isPivot(G.upperNodeIdBound()),
		minPivotDist(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
		maxPivotDist(G.upperNodeIdBound(), 0),
		nearestPivot(G.upperNodeIdBound()),
		farestPivot(G.upperNodeIdBound()),
		farnessToPivots(G.upperNodeIdBound()),
		minFarness(std::numeric_limits<edgeweight>::max()),
		pivots(G.upperNodeIdBound()),
		numberOfLandmarks(G.numberOfNodes()),
		maxComparisons(G.numberOfNodes()),
		method(method)
		{};

CustomizedAPSP::CustomizedAPSP(const Graph& G, const count numberOfLandmarks, count maxComparisons, CentralNodeMethod method): 	G(G),
		relevantPairs(G.upperNodeIdBound(), true),
		distances(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
		ecc(G.upperNodeIdBound()),
		isPivot(G.upperNodeIdBound()),
		minPivotDist(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
		maxPivotDist(G.upperNodeIdBound(), 0),
		nearestPivot(G.upperNodeIdBound()),
		farestPivot(G.upperNodeIdBound()),
		farnessToPivots(G.upperNodeIdBound()),
		minFarness(std::numeric_limits<edgeweight>::max()),
		pivots(numberOfLandmarks),
		numberOfLandmarks(numberOfLandmarks),
		maxComparisons(maxComparisons),
		method(method)
		{};

void CustomizedAPSP::run(){
	computeCentralNodes();

	//#pragma omp parallel for
	for(auto i = 0; i < pivots.size(); ++i){
		const node source = pivots[i];
		count bound = G.upperNodeIdBound();
		edgeweight _eccentricity = 0;

		std::queue<node> q;
		distances.set(source,source,0);

		std::vector<bool> visited(bound, false);
		q.push(source);
		visited[source] = true;

		while (!q.empty()) {
			node u = q.front();
			q.pop();
			bool relevant = true;
			G.forNeighborsOf(u, [&](node v) {
				edgeweight d = distances.element(source,u)+1;
				if (!visited[v]) {
					q.push(v);
					visited[v] = true;
					distances.set(source,v, d);
					if(d > maxPivotDist[v]){
						farestPivot[v] = source;
					}
					maxPivotDist[v] = std::max(maxPivotDist[v], d);
					if (d < minPivotDist[v]) {
						minPivotDist[v] = d;
						nearestPivot[v] = source;
					}
					//update eccentricity
					if(_eccentricity < d){
						_eccentricity = d;
					}
					farnessToPivots[v] += d;
					//(source, u) not far apart
					relevant = false;
				}else{
					// v was visited but is in the (i+1)-th hierarchy level ----> (s,u) can not be far apart
					if(distances.element(source, v) == distances.element(source, u)+1){
						//(source, u) not far apart
						relevant = false;
					}
				}
			});
			if(!relevant){
				relevantPairs.set(source,u, false);
			}
		}

		ecc[source] = _eccentricity;
	}

	hasRun = true;
}

edgeweight CustomizedAPSP::distUpperBoundByComp(node u, node v) {

	if(u==v){
		return 0.0;
	}

	if(G.hasEdge(u,v)){
		return 1.0;
	}

	if(distances.element(u,v) != std::numeric_limits<edgeweight>::max()){
		return distances.element(u,v);
	}

	//d(u,v) <= min_{j central node} d(u,j) + d(j,v)
	edgeweight smallestSeen = std::numeric_limits<edgeweight>::max();
	for(auto i=0; i< maxComparisons; ++i){
		edgeweight upperBound = distances.element(u,pivots[i]) + distances.element(pivots[i],v);
		if(upperBound < smallestSeen){
			smallestSeen = upperBound;
		}
	}

	return smallestSeen;
}

edgeweight CustomizedAPSP::distLowerBoundByComp(node u, node v) {

	if(u==v){
		return 0.0;
	}

	if(G.hasEdge(u,v)){
		return 1.0;
	}

	if(distances.element(u,v) != std::numeric_limits<edgeweight>::max()){
		return distances.element(u,v);
	}

	//d(u,v) <= min_{j central node} d(u,j) + d(j,v)
	edgeweight largestSeen = 0;
	for(auto i=0; i< maxComparisons; ++i){
		edgeweight lowerBound = std::abs(distances.element(u,pivots[i]) - distances.element(pivots[i],v));
		if(lowerBound > largestSeen){
			largestSeen = lowerBound;
		}
	}

	return largestSeen;
}

edgeweight CustomizedAPSP::distUpperBound(node u, node v) {

	if (u == v) {
		return 0;
	}

	return std::min(maxPivotDist[u] + minPivotDist[v],
			minPivotDist[u] + maxPivotDist[v]);
}

edgeweight CustomizedAPSP::distUpperBoundAdvanced(node u, node v){
	if (u == v) {
		return 0;
	}

	return std::min(distances.element(u, nearestPivot[u]) + distances.element(nearestPivot[u], v),
			distances.element(u, nearestPivot[v]) + distances.element(nearestPivot[v], v));
}

edgeweight CustomizedAPSP::distLowerBoundAdvanced(node u, node v){
	if(u==v){
		return 0;
	}

	return std::max(std::abs(distances.element(u, farestPivot[u]) - distances.element(farestPivot[u], v)),
			std::abs(distances.element(u, farestPivot[v]) - distances.element(farestPivot[v], v)));
}

void CustomizedAPSP::setBoundedDistances(){
	for(node u=0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			distances.set(u,v, distUpperBoundByComp(u,v));
		}
	}
}

SymMatrix<edgeweight, node> const& CustomizedAPSP::getExactDistances(){
	if(pivots.size() < G.numberOfNodes()){
		throw std::runtime_error("Exact distances can only be provided if all nodes of the G are picked as central nodes!");
	}
	return distances;
}


void CustomizedAPSP::computeCentralNodes(){
	auto get_wall_time = []()->double{
		struct timeval time;
		if(gettimeofday(&time, NULL)){
			return 0;
		}

		return (double) time.tv_sec + (double) time.tv_usec * .000001;
	};

	switch(method) {
		case CentralNodeMethod::NATURAL:
			{
					pivots = G.nodes();
					pivots.resize(numberOfLandmarks);
					for(auto i=0; i<numberOfLandmarks; ++i){
						isPivot[pivots[i]]=true;
					}
			}
			break;
		case CentralNodeMethod::DEGREE:
			{
				double startDeg = get_wall_time();
				DegreeCentrality degCen(G);
				degCen.run();
				auto nodesByDegree = degCen.ranking();
				for(auto i=0; i<numberOfLandmarks; ++i){
					pivots.push_back(nodesByDegree[i].first);
					isPivot[nodesByDegree[i].first] = true;
				}
				double endDeg = get_wall_time();
				if(G.numberOfNodes() >= 1000){
				INFO("DEGREE CENTRALITIY Running Time within CustomizedAPSP: ", endDeg-startDeg);
				}
			}
			break;
		case CentralNodeMethod::TOPCLOSENESS:
			{
				double startTOP = get_wall_time();
				TopCloseness topK(G, numberOfLandmarks);
				topK.run();
				pivots = topK.topkNodesList();
				for(auto i=0; i<numberOfLandmarks; ++i){
					isPivot[pivots[i]]=true;
				}
				double endTOP = get_wall_time();
				if(G.numberOfNodes() >= 1000){
				INFO("TOPCLOSENESS Running Time within CustomizedAPSP: ", endTOP-startTOP);
				}
			}
			break;
		case CentralNodeMethod::RANDOM:
			{
				auto rng = std::default_random_engine{};
				pivots = G.nodes();
				std::shuffle(pivots.begin(), pivots.end(), rng /*Aux::Random::getURNG()*/);
				pivots.resize(numberOfLandmarks);
				break;
			}
			break;
		case CentralNodeMethod::GROUPDEGREE:
		{
			double startGDEG = get_wall_time();
			GroupDegree groupDeg(G, numberOfLandmarks);
			groupDeg.run();
			std::vector<node> nodesByGroupDegree = groupDeg.groupMaxDegree();
			for(count i=0; i<numberOfLandmarks; ++i){
				//INFO(nodesByGroupDegree[i]);
				pivots[i] = nodesByGroupDegree[i];
				isPivot[nodesByGroupDegree[i]] = true;
			}
			double endGDEG = get_wall_time();
			if(G.numberOfNodes() >= 1000){
			INFO("GROUPDEGREE CENTRALITIY Running Time within CustomizedAPSP: ", endGDEG-startGDEG);
			}
		}
	}
}

}

