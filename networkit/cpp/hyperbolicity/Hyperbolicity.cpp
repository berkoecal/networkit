/*
 * Hyperbolicity.cpp
 *
 *  Created on: 05.11.2017
 *      Author: Berk Öcal
 */
#include<set>
#include<iostream>
#include<algorithm>
#include <queue>
//#include <mutex>
#include "Hyperbolicity.h"
#include "../distance/APSP.h"
#include "../components/ConnectedComponents.h"
#include "../auxiliary/Log.h"
#include <valgrind/callgrind.h>
#include "CustomizedBFS.h"

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

	CustomizedBFS cBFS(graph);
	cBFS.run();
	SymMatrix<bool, node>  const& farApartPairs = cBFS.getFarApartPairs();
	SymMatrix<edgeweight, node> const& distances = cBFS.getDistances();

	/*
	* Create NodeTupleWithDist objects only in appropriate order of nodes:
	* allowed: (2,5) not allowed: (5,2)
	*/
	
	//Add all far-apart pairs to the set
	for(node u=0; u < graph.numberOfNodes(); ++u){
		for(node v = u; v < graph.numberOfNodes(); ++v){
			if(farApartPairs.element(u,v)){
				ordered_tuples.emplace(u,v, distances.element(u,v));
			}
		}
	}

	INFO("Number of Far Apart Pairs: ", ordered_tuples.size(), " with percantage of ", static_cast<double>(ordered_tuples.size()/static_cast<double>((graph.numberOfNodes()*graph.numberOfNodes()/2))));

	/*
	 * Compute central node according to closeness centrality.
	 * Since we already have computed pairwise distances there is no need to apply the closeness class.
	 */

	std::clock_t start_c;
    double duration_c;
    start_c = std::clock();
	node central_node = centralNode(distances);
	
    duration_c = (std::clock() - start_c) / (double) CLOCKS_PER_SEC;

    INFO("Time of Central Node Computation: ", duration_c);

	edgeweight h_diff = 0;
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
    		if(is_valuable(distances,x,y,v, h_diff/2, cBFS.getEccentricity(v), central_node)){
    			for(const node w: mate[v]){
    				if(is_acceptable(distances, x,y,w, h_diff/2, cBFS.getEccentricity(w))){

    					//calculation of maximum and updating hyperbolicity: max(h_diff, S1 - max(S2,S3))
    					h_diff = std::max(h_diff, it_1->distance() + distances.element(v,w) -
    							std::max(distances.element(x,v) + distances.element(y,w),
    									distances.element(x,w) + distances.element(y,v)));

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

bool Hyperbolicity::is_acceptable(SymMatrix<edgeweight, node> const& distances,
				  node const& x, node const& y, node const& v, 
				  edgeweight const& current_lower_bound, 
				  edgeweight const& eccentricity)
{
	edgeweight const& distXV = distances.element(x,v);
	edgeweight const& distYV = distances.element(y,v);
	edgeweight const& distXY = distances.element(x,y);
	
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

bool Hyperbolicity::is_valuable(SymMatrix<edgeweight, node> const& distances,
				node const& x, node const& y, node const& v, 
				edgeweight const& current_lower_bound, 
				edgeweight const& eccentricity,
				node const& central_node
 			      )
{
	if(is_acceptable(distances, x, y, v, current_lower_bound, eccentricity)){
	    edgeweight const& distXV = distances.element(x,v);
	    edgeweight const& distYV = distances.element(y,v);
	    edgeweight const& distXY = distances.element(x,y);
	    edgeweight const& distVC = distances.element(central_node,v);
	    if((distXY- distXV - distYV + distVC)/2 > current_lower_bound){
		return true;
	    }
	}
	
	return false;
}

//TODO parallelize with parallelSumForNodes

node Hyperbolicity::centralNode(SymMatrix<edgeweight, node> const& distances){
	edgeweight minimum_value = std::numeric_limits<edgeweight>::max();
	//dummy initialization
	node most_central_node=0;
	for(node u = 0; u < graph.numberOfNodes(); ++u){
	    //dummy initialization
	    edgeweight sum = 0;
	    for(node v = u; v < graph.numberOfNodes(); ++v){
			sum+= distances.element(u,v);
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