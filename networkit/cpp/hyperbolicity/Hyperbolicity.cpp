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
#include "SymmetricMatrix.h"
#include "../centrality/TopCloseness.h"
#include "../centrality/DegreeCentrality.h"
#include "CustomizedAPSP.h"
#include "../distance/Bisenius/PrunedAPSP.h"

namespace NetworKit{

class NodeTupleWithDist
{
public:
	NodeTupleWithDist(node index_a, node index_b, /*count degree_a, count degree_b,*/ edgeweight dist_) :
		indexA(index_a),
		indexB(index_b),
		/*degreeA(degree_a),
		degreeB(degree_b),*/
		dist(dist_){}

	node index_a() const{return indexA;}
	node index_b() const{return indexB;}
	/*count degree_a() const{return degreeA;}
	count degree_b() const{return degreeB;}*/
	edgeweight distance() const {return dist;}

private:
	node indexA;
	node indexB;
	/*count degreeA;
	count degreeB;*/
	edgeweight dist;
};


Hyperbolicity::Hyperbolicity(const Graph& G): Algorithm(), G(G) {}

void Hyperbolicity::run(){
	if(G.numberOfNodes() < 4){
		hyperbolicity_value = 0;
		hasRun=true;
	}else{
		ConnectedComponents cc(G);
		cc.run();
		/*
		 * We (temporarily) assume the graph to be connected
		 */
		if(cc.numberOfComponents() > 1){throw std::runtime_error("Graph is not connected");}
		HYP();
		INFO("---------------------------------------");
		hasRun = true;
	}
}


void Hyperbolicity::HYP_AKIBA(){
	std::clock_t startALL;
	double durationALL;
	startALL = std::clock();

	auto lex_comp = [](NodeTupleWithDist const& a, NodeTupleWithDist const& b)-> bool {
		edgeweight res = a.distance() - b.distance();
		//If same distance then consider (increasing) lexicographic order
		if(res == 0){
	        int diff = a.index_a() - b.index_a();
	        if (diff == 0)
	        {
	            diff =  a.index_b() - b.index_b();
	        }
	        // (x,y) is larger than (v,w) if x < v or (x=v and y < w), and smaller otherwise
	        return diff < 0;
		}
		return res > 0;
	};

	auto ordered_tuples = std::set<NodeTupleWithDist, decltype(lex_comp)>(lex_comp);

	count n = G.numberOfNodes();
	count numOfCentralNodes = 20;
	auto k = std::min(n, numOfCentralNodes);

	std::clock_t start;
	double duration;
	start = std::clock();

	//Upper bounds
	CustomizedAPSP upperBounds(G, k, CentralNodeMethod::TOPCLOSENESS);
	upperBounds.run();
	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	INFO("Time of CustomizedAPSP", duration);

	count numOfRelPairs=0;
	SymMatrix<bool, node> const& relPairs = upperBounds.getRelevantPairs();

	std::clock_t startTuples;
	double durationTuples;
	startTuples = std::clock();

	for(node u=0; u < G.numberOfNodes(); ++u){
		for(node v = u; v < G.numberOfNodes(); ++v){
			if(relPairs.element(u,v)){
				numOfRelPairs++;
				ordered_tuples.emplace(u,v, upperBounds.getDistance(u,v));
			}
		}
	}

	durationTuples = (std::clock() - startTuples) / (double) CLOCKS_PER_SEC;
	INFO("Duration Tuple Construction ", durationTuples);
	INFO("Number of RelPairs: ", numOfRelPairs, " with percentage of ", (numOfRelPairs/static_cast<double>((G.numberOfNodes()*G.numberOfNodes()/2)))*100,"%");

	/*	TODO
	 *  1) 	Approximiere diam(G) schnell (double sweep oder 20(?) central nodes)!	Setze zunächst U=approxDiam(G), L=approxDiam(G)/2
	 *  2) 	Führe PrunedAPSP aus, d.h es wird zuerst das Prepocessing durchgeführt. Lösche im Preprocessing bereits alle NICHT far-apart Paare
	 *  	Bemerke: Der Preprocessing Schritt besucht nur sehr wenige Paare (ca. 0%-5% aller Paare des Graphen) -> es werden möglicherweise nur sehr
	 *  			 wenige eliminiert.
	 *  3)  Betrachte anschließend nur sortierte Paare P= {(x_0, y_0), (x_1,y_1), ... , (x_k, y_k) | (wobei alle distanzen ≥ L)}
	 *  4)	Prüfe für jedes Paar (x,y) in P während der Schleife, ob (x,y) far apart mit einem "guten" counternode (= ein zentraler Knoten(?))
	 *  5)  Abbruchkriterium: h_diff ≥ d(x_k, y_k) ≥ L
	 *  6) 	Setzt das Abbruchkriterium nicht an, so setze U=L und L= L - L/2 und finde alle Paare (x,y) mit U ≥ d(x,y) ≥ L und füge diese Liste (genannt P') P hinzu
	 *  	d.h P=P+P'
	 *  7) 	Iteriere obige Schritte 3) - 6)
	 */


	//Akiba
	std::clock_t start_preproc;
	double duration_preproc;
	start_preproc = std::clock();
	PrunedAPSP pruned(G);
	pruned.run();
	duration_preproc = (std::clock()-start_preproc)/(double) CLOCKS_PER_SEC;
	INFO("Time of Preprocessing Pruned: ", duration_preproc);

	//INFO("Number of RelPairs: ", ordered_tuples.size(), " with percentage of ", (ordered_tuples.size()/static_cast<double>((G.numberOfNodes()*G.numberOfNodes()/2)))*100,"%");

	// central node via TopCloseness
//	TopCloseness topCentral(G, 20);
//	topCentral.run();
//	auto list = topCentral.topkNodesList();
//	node c = list[0];
//
//	auto isNotFarApart([&](node u, node v, const double distUV) -> bool{
//		//for(auto i=0; i< list.size(); ++i){
//		//for(node w: G.neighbors(u)){
//			node w = 0;
//			auto const distUW = pruned.getDistance(u,w);
//			auto const distVW = pruned.getDistance(v,w);
//			if((distUV + distUW == distVW) or
//					distUV + distVW == distUW){
//				return true;
//				//break;
//			}
//		//}
//		return false;
//	});
//
//		//Upper bounds
//		count counter = 0;
//		count counterLarge=0;
//		for(node u=0; u < G.numberOfNodes(); ++u){
//			for(node v = u+1; v < G.numberOfNodes(); ++v){
//				const double distUV = pruned.getDistance(u,v);
//				if(distUV >= 50){
//					counterLarge++;
//					if(!isNotFarApart(u,v, distUV)){
//						counter++;
//					}
//				}
//			}
//		}
//		count n = G.numberOfNodes();
//		INFO("Number of Large Distance Pairs: ", counterLarge, " with ",
//				(counter/static_cast<double>(n*(n-1)/2))*100);
//		INFO("Number of Far Apart Pairs according to new method: ", counter, " with ",
//				(counter/static_cast<double>(n*(n-1)/2))*100);

	//------------------------------------------MAIN PART ----------------------------------------------------
	INFO("Hyperbolicity:", 0.0);
	durationALL = (std::clock() - startALL) / (double) CLOCKS_PER_SEC;
	INFO("Time of complete computation: ", durationALL);
	hyperbolicity_value = 0.0;
}

void Hyperbolicity::HYP(){
	std::clock_t startALL;
	double durationALL;
	startALL = std::clock();

	/*******************************************************************************/
	/* Preprocessing                                                             */
	/*******************************************************************************/

	auto comp = [](NodeTupleWithDist const& a, NodeTupleWithDist const& b)-> bool {
		edgeweight res = a.distance() - b.distance();
		//If same distance then consider (increasing) lexicographic order
		if(res == 0){
	        int diff = a.index_a() - b.index_a();
	        if (diff == 0)
	        {
	            diff =  a.index_b() - b.index_b();
	        }
	        // (x,y) is larger than (v,w) if x < v or (x=v and y < w), and smaller otherwise
	        return diff < 0;
		}
		return res > 0;
	};

	auto ordered_tuples = std::set<NodeTupleWithDist, decltype(comp)>(comp);

	std::clock_t startX;
	double durationX;
	startX = std::clock();
	CustomizedAPSP cAPSP(G);
	cAPSP.run();
	durationX = (std::clock() - startX) / (double) CLOCKS_PER_SEC;
	INFO("Time of CustomizedAPSP: ", durationX);

	SymMatrix<bool, node>  const& relevantPairs = cAPSP.getRelevantPairs();
	//std::set<NodeTupleWithDist, CustomCompare> const& ordered_tuples = cAPSP.getRelPairs();
	SymMatrix<edgeweight, node> const& distances = cAPSP.getExactDistances();

	std::clock_t startTuples;
	double durationTuples;
	startTuples = std::clock();
	//Diesen Schritt auslassen und direkt in CustomizedAPSP die Tupel erzeugen
	//Add all far-apart pairs to the set
	for(node u=0; u < G.numberOfNodes(); ++u){
		for(node v = u; v < G.numberOfNodes(); ++v){
			if(relevantPairs.element(u,v)){
				ordered_tuples.emplace(u,v, distances.element(u,v));
			}
		}
	}
	durationTuples = (std::clock() - startTuples) / (double) CLOCKS_PER_SEC;
	INFO("Time of Tuple Construction: ", durationTuples);

	INFO("Diameter: ", ordered_tuples.begin()->distance());
	INFO("Number of Far Apart Pairs: ", ordered_tuples.size(),
			" with percantage of ", static_cast<double>(ordered_tuples.size()/static_cast<double>((G.numberOfNodes()*G.numberOfNodes()/2)))*100,"%");

	std::clock_t start_c;
	double duration_c;
	start_c = std::clock();

	// central node via TopCloseness
	TopCloseness topCentral(G, 1);
	topCentral.run();
	node central_node = topCentral.topkNodesList()[0];

	duration_c = (std::clock() - start_c) / (double) CLOCKS_PER_SEC;
	INFO("Time of Central Node Computation: ", duration_c);



//------------------------------------------MAIN PART ----------------------------------------------------

	double h_diff = 0.0;
	std::set<node> was_seen;
	std::vector<std::vector<node>> mate(G.upperNodeIdBound());

	std::clock_t start;
	double duration;
	start = std::clock();


	size_t iteration = 0;
	node best_x=0;
	node best_y=0;
	node best_v=0;
	node best_w=0;
	auto degree_x=0;
	auto degree_y=0;
	auto degree_v=0;
	auto degree_w=0;
	count max_degree_until_abortion=0;

	size_t position_of_best=0;

	for(auto it_1 = ordered_tuples.begin(); it_1 != ordered_tuples.end(); ++it_1){
		if(it_1->distance() <= h_diff){
			goto Ende;
		}

		//INFO(it_1->distance(),"  (", it_1->index_a(), ",", it_1->index_b(), ")");

		iteration++;
		node x = it_1->index_a();
		node y = it_1->index_b();

		if(G.degree(x) > max_degree_until_abortion){
			max_degree_until_abortion = G.degree(x);
		}
		if(G.degree(y) > max_degree_until_abortion){
			max_degree_until_abortion = G.degree(y);
		}

		//Lemma5
		for(auto const v: was_seen){
			if(is_valuable(distances,x,y,v, h_diff, cAPSP.getEccentricity(v), central_node) ){
				for(const node w: mate[v]){ //what about quadtrupels (x,y,v,w) where e.g (x,v,y,w) was already considered before
					if(is_acceptable(distances, x,y,w, h_diff, cAPSP.getEccentricity(w))){
						double S1 = it_1->distance() + distances.element(v,w);
						double S2 = distances.element(x,v) + distances.element(y,w);
						double S3 = distances.element(x,w) + distances.element(y,v);
						double largest_diff = std::max(std::max(S1, S2), S3) - std::max(std::min(S1,S2),S3);
						if(h_diff < largest_diff){
							best_x = x;
							best_y = y;
							best_v = v;
							best_w = w;
							degree_x = G.degree(x);
							degree_y = G.degree(y);
							degree_v = G.degree(v);
							degree_w = G.degree(w);
							position_of_best = iteration;
							h_diff = largest_diff;
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
	INFO("Max degree until abortion: ", max_degree_until_abortion);
	INFO("Laufzeit vom Hauptblock:", duration);
	INFO("Number of pairs considered (Abortion Position): ", iteration);
	INFO("First quadtrupel attaining maximum:", "(",best_x,",",best_y,",",best_v, ",",best_w, ") with iteration position of (x,y) at ", position_of_best, "/", ordered_tuples.size(),
			" and distances d(x,y)= ", distances.element(best_x,best_y), ", d(v,w)= ", distances.element(best_v,best_w),
			", d(x,v)= ",distances.element(best_x,best_v), ", d(y,w)= ", distances.element(best_y, best_w),
			", d(x,w)= ",distances.element(best_x,best_w), ", d(y,v)= ", distances.element(best_y,best_v));
	INFO("Degrees: ", degree_x, " ", degree_y, " ", degree_v, " ", degree_w);
	INFO("h_diff: ", h_diff);

	INFO("Hyperbolicity: ", h_diff/2);
	durationALL = (std::clock() - startALL) / (double) CLOCKS_PER_SEC;
	INFO("Time of complete computation: ", durationALL);
	hyperbolicity_value = h_diff/2;
}

bool Hyperbolicity::is_acceptable(SymMatrix<edgeweight, node> const& distances,
				  node const& x, node const& y, node const& v, 
				  edgeweight const& current_max_diff, 
				  edgeweight const& eccentricity)
{
	edgeweight const& distXV = distances.element(x,v);
	edgeweight const& distYV = distances.element(y,v);
	edgeweight const& distXY = distances.element(x,y);
	
	//Corollary 7
	if(distXV <= current_max_diff/2. or distYV <= current_max_diff/2.){
	    return false;
	}
	//Lemma8
	if(2*eccentricity - distXV - distYV < 2*current_max_diff + 2 - distXY){
	    return false;
	}
	//Lemma 9
	if(eccentricity + distXY - 3*(current_max_diff/2.) - 3/2. < std::max(distXV, distYV)){
	    return false;
	}
	
	//INFO("v=",v, " is acceptable!");

	return true;
}

bool Hyperbolicity::is_valuable(SymMatrix<edgeweight, node> const& distances,
		node const& x, node const& y, node const& v,
		edgeweight const& current_max_diff,
		edgeweight const& eccentricity,
		node const& central_node
)
{
	if(is_acceptable(distances, x, y, v, current_max_diff, eccentricity)){
		edgeweight const& distXV = distances.element(x,v);
		edgeweight const& distYV = distances.element(y,v);
		edgeweight const& distXY = distances.element(x,y);
		edgeweight const& distVC = distances.element(central_node,v);

		if((distXY- distXV - distYV)/2. + distVC > current_max_diff/2.){
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
	for(node u = 0; u < G.numberOfNodes(); ++u){
	    //dummy initialization
	    edgeweight sum = 0;
	    for(node v = u; v < G.numberOfNodes(); ++v){
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
	APSP allPairShortestPaths(G);
	allPairShortestPaths.run();
	
	/*
	* Create NodeTupleWithDist objects only in appropiate order of nodes,
	*/
	
	for(node i=0; i < G.numberOfNodes(); ++i){
	  for(node j=i+1; j < G.numberOfNodes(); ++j){
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
