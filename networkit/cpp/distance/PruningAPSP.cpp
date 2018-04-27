/*
 * PruningAPSP.cpp
 *       Created on: 01.12.2017
 *       Author: Berk Öcal
 */

#include "PruningAPSP.h"

#include <cassert>
#include <mutex>

namespace NetworKit{

std::mutex L_mutex;
  
PruningAPSP::PruningAPSP(const Graph& G):
      Algorithm(), 
      G(G), 
      L(G.upperNodeIdBound(), std::vector<NodeWithDistance>()),
      distances(G.upperNodeIdBound())
{
  
}

void PruningAPSP::run()
{				
	//TODO order vertices by degree	
	assert(not (G.isDirected() or G.isWeighted()));
	
	//preprocessing (i.e. labeling)
	for(node v=0; v < G.upperNodeIdBound(); ++v){
	    prunedBFS(v);
	}
	

	hasRun = true;
}

void PruningAPSP::runParallel()
{
    assert(not (G.isDirected() or G.isWeighted()));
    
    G.parallelForNodes([&](node v){
	prunedBFS(v);
    });
    
    hasRun = true;
}

void PruningAPSP::prunedBFS(node v)
{
	edgeweight const infDist = std::numeric_limits<edgeweight>::max();
	std::queue<node> q;
	std::vector<edgeweight> P(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max());
	
	P[v] = 0;
	q.push(v);
	std::vector<edgeweight> T = computeDistanceVector(v);
	T[v] = 0;
	

	while(not q.empty()){
	    node const u = q.front();
	    q.pop();
	   
	    //pruning step
	    edgeweight const queryValue = mergingComputation(T,u);
	    if(queryValue <= P[u]){
	      continue;
	    }	

	    L_mutex.lock();  
	    L[u].emplace_back(v, P[u]);
	    L_mutex.unlock();
	    
	    std::vector<node> neighbors = G.neighbors(u);
	    for(auto const& neighbor: neighbors){
	      if(P[neighbor] == infDist){
		  assert(neighbor < G.upperNodeIdBound());
		  P[neighbor] = P[u]+1;
		  q.push(neighbor);
	      }
	    } 
	}
			
}

std::vector<edgeweight> PruningAPSP::computeDistanceVector(node v) const{
	std::vector<edgeweight> T(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max());
	L_mutex.lock();
	for(auto const& nodedistpair : L[v]){
	  node w = nodedistpair.getNode();
	  T[w] = nodedistpair.getDistance();
	}
	L_mutex.unlock();
	return T;
} 

edgeweight PruningAPSP::mergingComputation(std::vector<edgeweight> const& T, node u) const{
	edgeweight minimum = std::numeric_limits<edgeweight>::max();
	L_mutex.lock();
	for(auto const& nodedistpair : L[u]){
	  node w = nodedistpair.getNode();
	  if(nodedistpair.getDistance() + T[w]<= minimum){
	    minimum = nodedistpair.getDistance() + T[w];
	  }
	}
	L_mutex.unlock();
	return minimum;
}

edgeweight PruningAPSP::mergingComputation(node v, node u) const{
	std::vector<edgeweight> const T = computeDistanceVector(v);	
	return mergingComputation(T, u); 
}
}
