/*
 * Hyperbolicity.cpp
 *
 *  Created on: 05.11.2017
 *      Author: Berk Öcal
 */
#include<set>
#include<iostream>
#include<algorithm>
#include "Hyperbolicity.h"
#include "../distance/APSP.h"
#include "../components/ConnectedComponents.h"
#include "../auxiliary/Log.h"
#include <valgrind/callgrind.h>
#include <queue>

namespace NetworKit{

class NodeTupleWithDist
{
public:

	NodeTupleWithDist(node index_a, node index_b, edgeweight dist_): indexA(index_a), indexB(index_b), dist(dist_)
	{
	  //catch errors beforehand
	  //assert(index_a < index_b);
	  //assert
	}
	node index_a() const{return indexA;}
	node index_b() const{return indexB;}
	edgeweight distance() const {return dist;}

private:
	node indexA;
	node indexB;
	edgeweight dist;
};

//class SpecialBFS : public SSSP{
//
//public:
//	SpecialBFS(const Graph& G, node source, bool storePaths=true, bool storeNodesSortedByDistance=false, node target = none): SSSP(G, source, storePaths, storeNodesSortedByDistance, target){
//	}
//
//	virtual void run(){
//		edgeweight infDist = std::numeric_limits<edgeweight>::max();
//			count z = G.upperNodeIdBound();
//			distances.clear();
//			distances.resize(z, infDist);
//
//			std::vector<bool> visited;
//			visited.resize(z, false);
//
//			if (storePaths) {
//				previous.clear();
//				previous.resize(z);
//				npaths.clear();
//				npaths.resize(z, 0);
//				npaths[source] = 1;
//			}
//
//			if (storeNodesSortedByDistance) {
//				std::vector<node> empty;
//				std::swap(nodesSortedByDistance, empty);
//			}
//
//			std::queue<node> q;
//			q.push(source);
//			visited[source] = true;
//			distances[source] = 0;
//			bool breakWhenFound = (target != none);
//			while (! q.empty()) {
//				node u = q.front();
//				q.pop();
//
//				if (storeNodesSortedByDistance) {
//					nodesSortedByDistance.push_back(u);
//				}
//				if (breakWhenFound && target == u) {
//					break;
//				}
//
//				// TRACE("current node in BFS: " , u);
//		//		TRACE(distances);
//
//				bool is_leaf_of_BFS = true;			//check if u is leaf of our BFS-tree, i.e. neighbors are either visited or not existent.
//
//				// insert untouched neighbors into queue
//				G.forNeighborsOf(u, [&](node v) {
//					// TRACE("scanning neighbor ", v);
//					if (!visited[v]) {
//						is_leaf_of_BFS = false;
//						q.push(v);
//						visited[v] = true;
//						distances[v] = distances[u] + 1;
//
//						if (storePaths) {
//							previous[v] = {u};
//							npaths[v] = npaths[u];
//						}
//					} else if (storePaths && (distances[v] == distances[u] + 1)) {
//						previous[v].push_back(u); 	// additional predecessor
//						npaths[v] += npaths[u]; 	// all the shortest paths to u are also shortest paths to v now
//					}
//				});
//
//				if(is_leaf_of_BFS){
//					farApartPairs.emplace(source, u, distances[u]);
//				}
//			}
//
//	}
//
//	bool comparator(NodeTupleWithDist const& a, NodeTupleWithDist const& b){
//		return a.distance() >= b.distance();
//	}
//
//	std::set<NodeTupleWithDist> getFarApartPairs(){ return farApartPairs;}
//
//private:
//	std::set<NodeTupleWithDist> farApartPairs;
//};
//
//class SpecialAPSP: public Algorithm {
//
//public:
//	SpecialAPSP(const Graph& G): Algorithm(), G(G) {}
//
//	void run(){
//		std::vector<edgeweight> distanceVector(G.upperNodeIdBound(), 0.0);
//		distances.resize(G.upperNodeIdBound(), distanceVector);
//		eccentricity.resize(G.upperNodeIdBound(),0.0);
//
//		if (!G.isWeighted()) {
//			G.parallelForNodes([&](node u){
//				SpecialBFS bfs(G, u, true, true);
//				bfs.run();
//				distances[u] = bfs.getDistances();
//
//				auto sorted_nodes = bfs.getNodesSortedByDistance();
//				auto const& eccentricity_node = sorted_nodes.back();
//				eccentricity[u] = distances[u][eccentricity_node];
//
//				allFarApartPairs.insert(bfs.getFarApartPairs().begin(), bfs.getFarApartPairs().end());
//			});
//		}else{
//			throw std::runtime_error("Graph is not unweighted!");
//		}
//		hasRun = true;
//	}
//
//	std::vector<std::vector<edgeweight> > getDistances() const {if (!hasRun) throw std::runtime_error("Call run method first"); return distances;}
//
//	edgeweight getDistance(node u, node v) const {if (!hasRun) throw std::runtime_error("Call run method first"); return distances[u][v];}
//
//	edgeweight getEccentricity(node u) const {if (!hasRun) throw std::runtime_error("Call run method first"); return eccentricity[u];}
//
//	virtual bool isParallel() const override { return true; }
//
//protected:
//	const Graph& G;
//	std::vector<std::vector<edgeweight> > distances;
//	std::vector<edgeweight> eccentricity;
//	std::set<NodeTupleWithDist> allFarApartPairs; //not sorted by distance
//};


Hyperbolicity::Hyperbolicity(const Graph& G): Algorithm(), graph(G) {}

void Hyperbolicity::run(){
  	ConnectedComponents cc(graph);
  	cc.run();	
       /* 
	* We (temporarily) assume the graph to be connected 
	*/     
	if(cc.numberOfComponents() > 1){throw std::runtime_error("Graph is not connected");}

	std::clock_t start;
    double duration;
    start = std::clock();

    hyperbolicity_value = HYP();

    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    INFO("Gesamtzeit von HYP:", duration);
    INFO("Hyperbolicity mit HYP: ", hyperbolicity_value);
	hasRun = true;
}

double Hyperbolicity::HYP(){

	/*
	* Generate set (ordered by distances) consisting of
	* pairs of nodes (called NodeTupleWithDist) and
	* compute pairwise distances.
	*/	
	auto comp = [](NodeTupleWithDist const& a, NodeTupleWithDist const& b)-> bool {return a.distance() >= b.distance();};
	auto ordered_tuples = std::set<NodeTupleWithDist, decltype(comp)>(comp);	
	APSP allPairShortestPaths(graph);
	allPairShortestPaths.run();
	auto distances = allPairShortestPaths.getDistances();
	/*
	* Create NodeTupleWithDist objects only in appropriate order of nodes:
	* allowed: (2,5) not allowed: (5,2)
	*/
	

	for(node u=0; u < graph.numberOfNodes(); ++u){
	  for(node v=u+1; v < graph.numberOfNodes(); ++v){
		  bool farApart = true;
		  for(node w=0; w < graph.numberOfNodes(); ++w){
			  if(w==u or w==v){
				  continue;
			  }
			  if(distances[u][w] == distances[v][w]+distances[v][u]
			  or distances[v][w] == distances[w][u]+distances[u][v]){
				  farApart=false;
				  break;
			  }
		  }
		  if(farApart){
			 ordered_tuples.emplace(u,v, distances[u][v]);
		  }
	  }
	}
	
	/*
	 * Compute central node according to closeness centrality.
	 * Since we already have computed pairwise distances there is no need to apply the closeness class.
	 */
	node central_node = centralNode(distances);
	
	edgeweight h_diff = 0;
	std::vector<bool> seen(graph.numberOfNodes(), false);
	std::set<node> was_seen;
	std::vector<std::vector<node>> mate(graph.numberOfNodes());

	std::clock_t start;
    double duration;
    start = std::clock();

	for(auto it_1 = ordered_tuples.begin(); it_1 != ordered_tuples.end(); ++it_1){
	    node x = it_1->index_a();
	    node y = it_1->index_b();

	    //TODO parallelisieren
	    //Lemma5
	    for(auto const v: was_seen){	
		  if(is_valuable(distances,x,y,v, h_diff/2, allPairShortestPaths.getEccentricity(v), central_node)){
		    for(const node w: mate[v]){
			if(is_acceptable(distances, x,y,w, h_diff/2, allPairShortestPaths.getEccentricity(w))){
			    //calculation of maximum and updating hyperbolicity
			    auto const sum1 = it_1->distance() + distances[v][w];
			    auto const sum2 = distances[x][v] + distances[y][w];
			    auto const sum3 = distances[x][w] + distances[y][v];
			    
			    edgeweight first = std::max(sum1,sum2);
			    edgeweight second = std::max(std::min(sum1,sum2), sum3);
			    
			    h_diff = std::max(h_diff, first-second);	    
			    if(it_1->distance() <= h_diff){
			    		goto Ende;
			    		return h_diff/2;
			    }
			}
		    }
		  }		
	    }
	    
	    mate[x].push_back(y);
	    mate[y].push_back(x);
	    was_seen.insert(x);
	    was_seen.insert(y);
	    
	}
	Ende:
	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	INFO("Laufzeit vom Hauptblock:", duration);
	
	return h_diff/2;
}

bool Hyperbolicity::is_acceptable(std::vector<std::vector<edgeweight>> const& distances, 
				  node const& x, node const& y, node const& v, 
				  edgeweight const& current_lower_bound, 
				  edgeweight const& eccentricity)
{
	edgeweight const& distXV = distances[x][v];
	edgeweight const& distYV = distances[y][v];
	edgeweight const& distXY = distances[x][y];
	
	//Corollary 7
	if(distXV <= current_lower_bound or distYV <= current_lower_bound){
	    return false;
	}
	//Lemma8
	if(2*eccentricity - distXV - distYV < 4*current_lower_bound + 2 - distXY){
	    return false;
	}
	//Lemma 9
	if(eccentricity + distXY - 3*current_lower_bound - 3/2 < std::max(distXY, distYV)){
	    return false;
	}
	
	return true;
}

bool Hyperbolicity::is_valuable(std::vector<std::vector<edgeweight>> const& distances, 
				node const& x, node const& y, node const& v, 
				edgeweight const& current_lower_bound, 
				edgeweight const& eccentricity,
				node const& central_node
 			      )
{
	if(is_acceptable(distances, x, y, v, current_lower_bound, eccentricity)){
	    edgeweight const& distXV = distances[x][v];
	    edgeweight const& distYV = distances[y][v];
	    edgeweight const& distXY = distances[x][y];
	    edgeweight const& distVC = distances[central_node][v];
	    if((distXY- distXV - distYV + distVC)/2 > current_lower_bound){
		return true;
	    }
	}
	
	return false;
}

node Hyperbolicity::centralNode(std::vector<std::vector<edgeweight>> const& distances){
	edgeweight minimum_value = std::numeric_limits<edgeweight>::max();
	//dummy initialization
	node most_central_node=0;
	for(node u = 0; u < graph.upperNodeIdBound(); ++u){
	    //dummy initialization
	    edgeweight sum = 0;
	    for(const edgeweight dist: distances[u]){
		sum+= dist;
	    }
	    if(minimum_value > sum){
	      minimum_value = sum;
	      most_central_node = u;
	    }
	}
	
	return most_central_node;
}



double Hyperbolicity::naiveAlgorithm(){
	
	//TODO consider all connected components
	//TODO Far-Apart Pairs
  
	/*
	 * Generate set (ordered by distances) consisting of
	 * pairs of nodes (called NodeTupleWithDist) and
	 * compute pairwise distances.
	 */	
	auto comp = [](NodeTupleWithDist const& a, NodeTupleWithDist const& b)-> bool {return a.distance() >= b.distance();};
	auto ordered_tuples = std::set<NodeTupleWithDist, decltype(comp)>(comp);	
	APSP allPairShortestPaths(graph);
	allPairShortestPaths.run();
	
	/*
	* Create NodeTupleWithDist objects only in appropiate order of nodes,
	*/
	
	for(node i=0; i < graph.numberOfNodes(); ++i){
	  for(node j=i+1; j < graph.numberOfNodes(); ++j){
	    ordered_tuples.emplace(i,j, allPairShortestPaths.getDistance(i,j));
	  }
	}
	
	edgeweight h_diff = 0;	
	
		
	for(auto it_1 = ordered_tuples.begin(); it_1 != ordered_tuples.end(); ++it_1){
	  for(auto it_2 = ordered_tuples.begin(); it_2 != it_1; ++it_2){
	    //Determine largest two sums and compare to current value of h_diff which 2*currently_largest_hyperbolicity
	    auto const s1 = it_1->distance() + it_2->distance(); // d(a,b) + d(c,d)
	    auto const s2 = allPairShortestPaths.getDistance(it_1->index_a(), it_2->index_a()) +  // d(a,c) + d(b,d)
	    allPairShortestPaths.getDistance(it_1->index_b(), it_2->index_b());
	    auto const s3 = allPairShortestPaths.getDistance(it_1->index_a(),it_2->index_b()) +   // d(a,d) + d(b,c)
	    allPairShortestPaths.getDistance(it_1->index_b(),it_2->index_a());   
	    
	    edgeweight first = std::max(s1,s2);
	    edgeweight second = std::max(std::min(s1,s2), s3);
	    
	    h_diff = std::max(h_diff, first-second);	    
	    if(it_1->distance() <= h_diff){
	      return h_diff/2;
	    }	  
	  }
	}	
	return h_diff/2;
}
}
