/*
 * CustomizedBFS.cpp
 *
 *  Created on: 21.02.2018
 *      Author: berkoec
 */

#include "CustomizedBFS.h"

namespace NetworKit{

CustomizedBFS::CustomizedBFS(const Graph& G): graph(G),
		relevantPairs(G.upperNodeIdBound(), true),
		boundedDistances(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
		eccentricity(G.upperNodeIdBound()),
		nodesToConsider(graph.nodes()){};

CustomizedBFS::CustomizedBFS(const Graph& G, const std::vector<node>& nodes): 	graph(G),
		relevantPairs(G.upperNodeIdBound(), true),
		boundedDistances(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
		eccentricity(G.upperNodeIdBound()),
		nodesToConsider(nodes){};

void CustomizedBFS::run(){
	//#pragma omp parallel for
	for(auto iter = nodesToConsider.begin(); iter != nodesToConsider.end(); ++iter){
		const node source = *iter;
		count bound = graph.upperNodeIdBound();
		edgeweight _eccentricity = 0.0;

		std::queue<node> q;
		boundedDistances.set(source,source,0.0);

		std::vector<bool> visited(bound, false);
		q.push(source);
		visited[source] = true;

		while (!q.empty()) {
			node u = q.front();
			q.pop();
			graph.forNeighborsOf(u, [&](node v) {
				if (!visited[v]) {
					q.push(v);
					visited[v] = true;
					edgeweight newDist = boundedDistances.element(source,u)+1;
					boundedDistances.set(source,v, newDist);
					//update eccentricity
					if(_eccentricity < newDist){
						_eccentricity = newDist;
					}
					//(source, u) not far apart
					relevantPairs.set(source,u, false);
				}else{
					// v was visited but is in the (i+1)-th hierarchy level ----> (s,u) can not be far apart
					if(boundedDistances.element(source, v) == boundedDistances.element(source, u)+1){
						//(source, u) not far apart
						relevantPairs.set(source,u, false);
					}
				}
			});
		}

		eccentricity[source] = _eccentricity;
	}

	setBoundedDistances();

	//TODO Parallelize!
	//TODO Abortion Criterion
//	graph.parallelForNodes([&](const node source){
//		count bound = graph.upperNodeIdBound();
//		edgeweight _eccentricity = 0;
//
//		std::queue<node> q;
//		boundedDistances.set(source,source,0);
//
//		std::vector<bool> visited(bound, false);
//		q.push(source);
//		visited[source] = true;
//
//		while (!q.empty()) {
//			node u = q.front();
//			q.pop();
//
//			graph.forNeighborsOf(u, [&](node v) {
//				if (!visited[v]) {
//					q.push(v);
//					visited[v] = true;
//					edgeweight newDist = boundedDistances.element(source,u)+1;
//					boundedDistances.set(source,v, newDist);
//					//update eccentricity
//					if(_eccentricity < newDist){
//						_eccentricity = newDist;
//					}
//					//(source, u) not far apart
//					relevantPairs.set(source,u, false);
//				}else{
//					// v was visited but is in the (i+1)-th hierarchy level ----> (s,u) can not be far apart
//					if(boundedDistances.element(source, v) == boundedDistances.element(source, u)+1){
//						//(source, u) not far apart
//						relevantPairs.set(source,u, false);
//					}
//				}
//			});
//		}
//
//		eccentricity[source] = _eccentricity;
//	});
}

SymMatrix<edgeweight, node> CustomizedBFS::getExactDistances(){
	if(nodesToConsider.size() < graph.numberOfNodes()){
		throw std::runtime_error("Exact distances can only be provided if all nodes of the graph are picked as central nodes!");
	}
	return boundedDistances;
}

edgeweight CustomizedBFS::computeBoundedDistance(node u, node v){
	if(u==v){
		return 0.0;
	}
	if(boundedDistances.element(u,v) != std::numeric_limits<edgeweight>::max()){
		return boundedDistances.element(u,v);
	}
	//d(u,v) <= min_{j central node} d(u,j) + d(j,v)
	edgeweight smallestSeen = std::numeric_limits<edgeweight>::max();
	for(auto iter=nodesToConsider.begin(); iter != nodesToConsider.end(); ++iter){
		const node j = *iter;
		edgeweight upperBound = boundedDistances.element(u,j) + boundedDistances.element(j,v);
		if(upperBound < smallestSeen){
			smallestSeen = upperBound;
		}
	}
	return smallestSeen;
}

void CustomizedBFS::setBoundedDistances(){
	for(node u = 0; u < graph.numberOfNodes(); ++u){
		for(node v=u; v < graph.numberOfNodes(); ++v){
			boundedDistances.set(u,v, computeBoundedDistance(u,v));
		}
	}
}

SymMatrix<bool,node> CustomizedBFS::getRelevantPairs(){return relevantPairs;}
SymMatrix<edgeweight, node> CustomizedBFS::getBoundedDistances(){return boundedDistances;}
edgeweight CustomizedBFS::getBoundedDistance(node u, node v){return boundedDistances.element(u,v);}
edgeweight CustomizedBFS::getEccentricity(node v){return eccentricity[v];}

}

