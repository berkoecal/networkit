/*
 * CustomizedBFS.cpp
 *
 *  Created on: 21.02.2018
 *      Author: berkoec
 */

#include "CustomizedBFS.h"

namespace NetworKit{

CustomizedBFS::CustomizedBFS(const Graph& G): 	graph(G),
												farApart(G.numberOfNodes(), true),
												distances(G.numberOfNodes(), std::numeric_limits<edgeweight>::max()),
												eccentricity(G.numberOfNodes()){};

void CustomizedBFS::run(){
	//TODO Parallelize!
	//TODO Abortion Criterion
	graph.parallelForNodes([&](const node source){
		count bound = graph.upperNodeIdBound();
		edgeweight _eccentricity = 0;

		std::queue<node> q;
		distances.set(source,source,0);

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
					edgeweight newDist = distances.element(source,u)+1;
					distances.set(source,v, newDist);
					//update eccentricity
					if(_eccentricity < newDist){
						_eccentricity = newDist;
					}
					//(source, u) not far apart
					farApart.set(source,u, false);
				}else{
					// v was visited but is in the (i+1)-th hierarchy level ----> (s,u) can not be far apart
					if(distances.element(source, v) == distances.element(source, u)+1){
						//(source, u) not far apart
						farApart.set(source,u, false);
					}
				}
			});
		}

		eccentricity[source] = _eccentricity;
	});
}

SymMatrix<bool,node> CustomizedBFS::getFarApartPairs(){return farApart;}

SymMatrix<edgeweight, node> CustomizedBFS::getDistances(){return distances;}

edgeweight CustomizedBFS::getDistance(node u, node v){ return distances.element(u,v);}

edgeweight CustomizedBFS::getEccentricity(node v){return eccentricity[v];}

}

