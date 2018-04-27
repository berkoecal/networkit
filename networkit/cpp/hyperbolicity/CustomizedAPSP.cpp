/*
 * CustomizedBFS.cpp
 *
 *  Created on: 21.02.2018
 *      Author: Berk Öcal
 */

#include "CustomizedAPSP.h"
#include "../centrality/DegreeCentrality.h"
#include "../centrality/TopCloseness.h"

namespace NetworKit{

CustomizedAPSP::CustomizedAPSP(const Graph& G, CentralNodeMethod method): graph(G),
		relevantPairs(G.upperNodeIdBound(), true),
		distances(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
		eccentricity(G.upperNodeIdBound()),
		numberOfLandmarks(G.numberOfNodes()),
		method(method)
		{};

CustomizedAPSP::CustomizedAPSP(const Graph& G, const count numberOfLandmarks, CentralNodeMethod method): 	graph(G),
		relevantPairs(G.upperNodeIdBound(), true),
		distances(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
		eccentricity(G.upperNodeIdBound()),
		numberOfLandmarks(numberOfLandmarks),
		method(method)
		{};

void CustomizedAPSP::run(){
	computeCentralNodes();

	for(auto i = 0; i < nodesToConsider.size(); ++i){
		const node source = nodesToConsider[i];
		count bound = graph.upperNodeIdBound();
		edgeweight _eccentricity = 0.0;

		std::queue<node> q;
		distances.set(source,source,0.0);

		std::vector<bool> visited(bound, false);
		q.push(source);
		visited[source] = true;

		while (!q.empty()) {
			node u = q.front();
			q.pop();
			bool relevant = true;
			graph.forNeighborsOf(u, [&](node v) {
				if (!visited[v]) {
					q.push(v);
					visited[v] = true;
					edgeweight newDist = distances.element(source,u)+1;
					distances.set(source,v, newDist);
					//update eccentricity
					if(_eccentricity < newDist){
						_eccentricity = newDist;
					}
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

		eccentricity[source] = _eccentricity;
	}

	//vorerst mit dieser if Abfrage
//	if(nodesToConsider.size() < graph.numberOfNodes()){
//		setBoundedDistances();
//	}
}

edgeweight CustomizedAPSP::computeBoundedDistance(node u, node v){
	if(u==v){
		return 0.0;
	}

//	if(graph.hasEdge(u,v)){
//		return 1.0;
//	}

	if(distances.element(u,v) != std::numeric_limits<edgeweight>::max()){
		return distances.element(u,v);
	}

	//d(u,v) <= min_{j central node} d(u,j) + d(j,v)
	edgeweight smallestSeen = std::numeric_limits<edgeweight>::max();
	for(auto i=0; i< nodesToConsider.size(); ++i){
		edgeweight upperBound = distances.element(u,nodesToConsider[i]) + distances.element(nodesToConsider[i],v);
		if(upperBound < smallestSeen){
			smallestSeen = upperBound;
		}
	}

	return smallestSeen;
}

void CustomizedAPSP::setBoundedDistances(){
	for(node u=0; u < graph.numberOfNodes(); ++u){
		for(node v=u; v < graph.numberOfNodes(); ++v){
			distances.set(u,v, computeBoundedDistance(u,v));
		}
	}
}

SymMatrix<edgeweight, node> const& CustomizedAPSP::getExactDistances() const{
	if(nodesToConsider.size() < graph.numberOfNodes()){
		throw std::runtime_error("Exact distances can only be provided if all nodes of the graph are picked as central nodes!");
	}
	return distances;
}


void CustomizedAPSP::computeCentralNodes(){
	switch(method) {
		case CentralNodeMethod::NATURAL:
			{
					nodesToConsider = graph.nodes();
					nodesToConsider.resize(numberOfLandmarks);
			}
			break;
		case CentralNodeMethod::DEGREE:
			{
				DegreeCentrality degCen(graph);
				degCen.run();
				auto nodesByDegree = degCen.ranking();
				for(auto i=0; i<numberOfLandmarks; ++i){
					nodesToConsider.push_back(nodesByDegree[i].first);
				}
			}
			break;
		case CentralNodeMethod::TOPCLOSENESS:
			{
				TopCloseness topK(graph, numberOfLandmarks);
				topK.run();
				nodesToConsider = topK.topkNodesList();
			}
			break;
		case CentralNodeMethod::RANDOM:
			{
				auto rng = std::default_random_engine{};
				nodesToConsider = graph.nodes();
				std::shuffle(nodesToConsider.begin(), nodesToConsider.end(), rng /*Aux::Random::getURNG()*/);
				nodesToConsider.resize(numberOfLandmarks);
				break;
			}
			break;
	}
}

SymMatrix<bool,node> CustomizedAPSP::getRelevantPairs() const{return relevantPairs;}
//std::set<NodeTupleWithDist, CustomCompare> const& CustomizedAPSP::getRelPairs() const {return relPairs;}
SymMatrix<edgeweight, node> CustomizedAPSP::getDistances() const{return distances;}
edgeweight CustomizedAPSP::getDistance(node u, node v) {return computeBoundedDistance(u,v);}
edgeweight CustomizedAPSP::getEccentricity(node v) const{return eccentricity[v];}

}

