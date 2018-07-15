/*
 * Hyperbolicity.cpp
 *
 *  Created on: 05.11.2017
 *      Author: Berk Öcal
 */

#define _BSD_SOURCE

#include<sys/time.h>
#include<set>
#include<iostream>
#include<algorithm>
#include <queue>
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
#include "../algebraic/CSRMatrix.h"
#include "../centrality/GroupDegree.h"
#include "../distance/PrunedLabeling.h"


namespace NetworKit{

class NodeTupleWithDist
{
public:
	NodeTupleWithDist(node index_a, node index_b, edgeweight dist_/*, edgeweight lowerBound_ = 0.0*/) :
		indexA(std::min(index_a, index_b)),
		indexB(std::max(index_a, index_b)),
		dist(dist_)/*,
		lowerB(lowerBound_)*/
{}

	node index_a() const{return indexA;}
	node index_b() const{return indexB;}
	edgeweight distance() const {return dist;}
	//edgeweight lowerBound() const {return lowerB;}

private:
	node indexA;
	node indexB;
	edgeweight dist;
	//edgeweight lowerB;
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
		HYP_AKIBA();
		//INFO("---------------------------------------");
		hasRun = true;
	}
}


void Hyperbolicity::HYP_AKIBA(){

	auto get_wall_time = []()->double{
		struct timeval time;
		if(gettimeofday(&time, NULL)){
			return 0;
		}

		return (double) time.tv_sec + (double) time.tv_usec * .000001;
	};


	/*******************************************************************************/
	/* PREPROCESSING                                                             */
	/*******************************************************************************/
	double startALL = get_wall_time();
	edgeweight inf = std::numeric_limits<edgeweight>::max();
	count n = G.numberOfNodes();

	/*
	 * This part corresponds to the selection of landmarks
	 * and the construction of upper bounds.
	 * @param k: number of landmarks
	 * @param maxComp: limit for U', i.e. it decided how many comparisons U' is allowed to make.
	 * @param CentralNodeMethod: can be picked as GROUPDEGREE, TOPCLOSENESS, RANDOM and DEGREE.
	 */
	count numOfCentralNodes = ceil(n*0.1);
	auto k = std::min(n, numOfCentralNodes);
	auto maxComp = std::min(numOfCentralNodes, (count) 20);

	INFO("n=", n, " m=", G.numberOfEdges());

	double start = get_wall_time();
	CustomizedAPSP customAPSP(G, k, maxComp, CentralNodeMethod::GROUPDEGREE);
	customAPSP.run();
	double end = get_wall_time();

	SymMatrix<bool, node> const& relPairsFromCAPSP = customAPSP.getRelevantPairs();
	SymMatrix<edgeweight, node>& distances = customAPSP.getDistances(); // stores exact distances
	/*
	 * These are comparison methods to sort the list of pairs P.
	 * If distances (or bounds) are equal, one can apply further ordering strategies.
	 */

	auto lex_comp = [](NodeTupleWithDist const& a, NodeTupleWithDist const& b)-> bool {
		edgeweight res = a.distance() - b.distance();
		if(res == 0){
			int diff = a.index_a() - b.index_a();
			if (diff == 0)
			{
				diff =  a.index_b() - b.index_b();
			}
			return diff < 0;
		}
		return res > 0;
	};

//	auto lower_bound_comp = [](NodeTupleWithDist const& a, NodeTupleWithDist const& b)-> bool {
//		edgeweight res = a.distance() - b.distance();
//		if(res == 0){
//			edgeweight lower = a.lowerBound() - b.lowerBound();
//			if(lower == 0){
//				int diff = a.index_a() - b.index_a();
//				if (diff == 0)
//				{
//					diff =  a.index_b() - b.index_b();
//				}
//				return diff < 0;
//			}
//			return lower > 0;
//		}
//
//		//hier werden alle tuple mit gleicher distance und gleichem lowerBound ausgeschlossen.
//		return res > 0;
//	};


	auto degree_comp = [&](NodeTupleWithDist const& firstTuple, NodeTupleWithDist const& secTuple)-> bool{
		edgeweight res = firstTuple.distance() - secTuple.distance();
		if(res == 0){
			count degA = G.degree(firstTuple.index_a()) + G.degree(firstTuple.index_b());
			count degB = G.degree(secTuple.index_a()) + G.degree(secTuple.index_b());
			return degA >= degB;
		}
		return res > 0;
	};
	std::vector<NodeTupleWithDist> ordered_tuples;

	//Akiba
	std::clock_t start_preproc;
	double duration_preproc;
	start_preproc = std::clock();
	PrunedLabeling<> pruned(G);
	pruned.ConstructIndex();
	duration_preproc = (std::clock()-start_preproc)/(double) CLOCKS_PER_SEC;

	//TopK
	double startCentral = get_wall_time();
	TopCloseness topCentral(G, 1);
	topCentral.run();
	auto central_node = topCentral.topkNodesList()[0];
	double endCentral = get_wall_time();
	if(G.numberOfNodes() >= 1000){
	INFO("Time of CustomizedAPSP with ", k, " nodes :", end-start);
	INFO("Time of central node computation: ", endCentral-startCentral);
	INFO("Time of Preprocessing Pruned: ", duration_preproc);
	}

	//This will be not necessary once the getDistances method in topCloseness is implemented
	BFS bfs(G, central_node);
	bfs.run();
	auto distC = bfs.getDistances();

	count capspRelPairs=0;

	double startTuples = get_wall_time();
/*
 * NEW with lower bounds
 */
//	for(node u=0; u < G.numberOfNodes(); ++u){
//		for(node v = u; v < G.numberOfNodes(); ++v){
//			if(relPairsFromCAPSP.element(u,v)){
//				capspRelPairs++;
//			}
//			auto dist = distances.element(u,v);
//			if(relPairsFromCAPSP.element(u,v)){
//				//maybe weaker upper bound, but fast to compute.
//				auto eccLower = std::min(customAPSP.eccLowerBound(u), customAPSP.eccLowerBound(v));
//				if(dist == inf){
//					dist = customAPSP.distUpperBoundAdvanced(u,v);
//					if(dist >= eccLower){
//						ordered_tuples.emplace(u,v,dist, customAPSP.distLowerBoundAdvanced(u,v));
//					}
//				}else{
//					if(dist >= eccLower){
//						ordered_tuples.emplace(u,v, dist, dist);
//					}
//				}
//			}
//		}
//	}

	for(node u=0; u < G.numberOfNodes(); ++u){
		for(node v = u; v < G.numberOfNodes(); ++v){
			if(relPairsFromCAPSP.element(u,v)){
				capspRelPairs++;
			}
			edgeweight dist = customAPSP.getDistance(u,v); //FEHLER VLLT HIER
			if(relPairsFromCAPSP.element(u,v)){
				if(dist == inf){
					dist = customAPSP.distUpperBoundAdvanced(u,v);
					//dist = customAPSP.distUpperBoundByComp(u,v);
					ordered_tuples.emplace_back(u,v,dist);
				}else{
					ordered_tuples.emplace_back(u,v, dist);
				}
			}
		}
	}

	std::sort(ordered_tuples.begin(), ordered_tuples.end(), lex_comp);

	double endTuples = get_wall_time();
	if(G.numberOfNodes() >= 1000){
	INFO("Duration Tuple Construction ", endTuples-startTuples);
	INFO("Number of CAPSP RelPairs: ", capspRelPairs, " with percentage of ", (capspRelPairs/static_cast<double>((G.numberOfNodes()*(G.numberOfNodes()-1)/2)))*100,"%");
	INFO("Number of RelPairs Total (ordered_tuples size): ", ordered_tuples.size(), " with percentage of ", (ordered_tuples.size()/static_cast<double>((G.numberOfNodes()*(G.numberOfNodes()-1)/2)))*100,"%");
	}

	/****************************************************************************
	 * LAMBDA EXPRESSIONS
	 ****************************************************************************/


	auto numProperQueries = 0;
	auto numTotalQueries = 0;
	double sumQueryTimes = 0;
	/*
	 * This method realizes Option A and B. Currently it applies Option A
	 * storing distances that are already queried.
	 * We can easily count the number of total queries and the number of proper queries, i.e.
	 * not taking repetitions into account. To count the number or time of (proper) queries,
	 * uncomment the corresponding variables.
	 */
	auto exact_dist = [&](node v, node w) -> edgeweight{
		//numTotalQueries++;
		if(distances.element(v,w) != inf){
			return distances.element(v,w);
		}
		//double startQ = get_wall_time();
		double dist = pruned.QueryDistance(v,w);
		//double endQ = get_wall_time();
		//numProperQueries++;
		distances.set(v,w, dist);
		//sumQueryTimes += (endQ - startQ);
		return dist;
	};


	/*******************************************************************************/
	/* MAIN PART                                                             */
	/*******************************************************************************/


	double h_diff = 0.0;
	std::vector<bool> was_seen(G.numberOfNodes(), false);
	std::vector<node> nodesToConsider(G.numberOfNodes());
	size_t nodesToConsiderCounter = 0;
	double wasSeenTime = 0;
	double accValTime = 0;
	std::vector<std::vector<node>> mate(n);
	std::vector<edgeweight> largestDistToMates(n,0);

	/* NEW */
	std::vector<bool> acceptable(n, 0);
	std::vector<node> valuable(n);
	/* ---- */

	std::clock_t start_main;
	double duration_main;
	start_main = std::clock();

	// ------ Tracking Data ------
	size_t iteration = 0;
	node best_x=0;
	node best_y=0;
	node best_v=0;
	node best_w=0;
	count quadtruples_visited=0;
	count position_of_best=0;
	//-----------------------------

	for(auto tuple = ordered_tuples.begin(); tuple != ordered_tuples.end(); ++tuple){
		if(tuple->distance() <= h_diff){
			goto Ende;
		}

		iteration++;
		node x = tuple->index_a();
		node y = tuple->index_b();

		edgeweight distXY = exact_dist(x,y);

		mate[x].push_back(y);
		mate[y].push_back(x);

		largestDistToMates[x] = std::max(largestDistToMates[x], distXY);
		largestDistToMates[y] = std::max(largestDistToMates[y], distXY);

		if(!(was_seen[x])){
			was_seen[x] = true;
			nodesToConsider[nodesToConsiderCounter++] = x;
		}

		if(!(was_seen[y])){
			was_seen[y] = true;
			nodesToConsider[nodesToConsiderCounter++] = y;
		}

		/*	Implements the measurement of
		 *  \phi_i and \gamma_i
		 */

//		static size_t debugCounter = 0;
//		INFO("Ratio", nodesToConsiderCounter/ (double) n,
//		" Progress: ", "ordered_tuples size: ", ordered_tuples.size(),
//		" debugCounter: ",debugCounter," Ratio: ",debugCounter++/ (double) ordered_tuples.size());

		/*************************************************************************
		 * Computation of acceptable and valuable nodes
		 ************************************************************************/

		double startAccVal = get_wall_time();

		auto hdiffPlusOne = h_diff + 1;
		count numberOfValNodes = 0;
		auto condval = distXY - h_diff;
		auto condacc1 = 3*hdiffPlusOne- 2*distXY;
		auto condacc2 = 2*hdiffPlusOne - distXY;

		if(condacc1 > 0){
			// see Lemma 9: if 3delta_L + 3/2 < d(x,y) or equivalently 3*h_diff+3 < 2*d(x,y) then by Lemma 9 node v is already NOT skippable
			//Lemma5
			for(size_t i=0; i < nodesToConsiderCounter; ++i){
				node v = nodesToConsider[i];
				//Lemma 9 with exact distances
				edgeweight distXV = exact_dist(x,v);
				edgeweight distYV = exact_dist(y,v);
				if(2*(largestDistToMates[v] - std::max(distXV,distYV)) >= condacc1){
					//Corollary 7 with exact distances
					if((2*distXV >= hdiffPlusOne) and (2*distYV >= hdiffPlusOne)){
						//Lemma 8 with exact distances
						if(2*largestDistToMates[v] - distXV - distYV >= condacc2){
							//node v is acceptable!
							acceptable[v] = true;
							//check if valuable
							if(condval + 2*distC[v] > distXV + distYV){
								valuable[numberOfValNodes++] = v;
							}
						}
					}
				}
			}
		}else{
			for(size_t i=0; i < nodesToConsiderCounter; ++i){
				node v = nodesToConsider[i];
				edgeweight distXV = exact_dist(x,v);
				edgeweight distYV = exact_dist(y,v);
				//Corollary 7 with exact distances
				if((2*distXV >= hdiffPlusOne) and (2*distYV >= hdiffPlusOne)){
					//Lemma 8 with exact distances
					if(2*largestDistToMates[v] - distXV - distYV >= condacc2){
						acceptable[v] = true;
						//check if valuable
						if(condval + 2*distC[v] > distXV + distYV){
							valuable[numberOfValNodes++] = v;
						}
					}
				}
			}
		}

		/*
		 * For Option B: Version with upper bounds
		 * This version of valuable and acceptable node computation seems to be unfeasible
		 * due to less exclusions by upper bounds. Better use the computation scheme above.
		 * In the thesis both versions are presented.
		 */

//		if(condacc1 > 0){
//			// see Lemma 9: if 3delta_L + 3/2 < d(x,y) or equivalently 3*h_diff+3 < 2*d(x,y) then by Lemma 9 node v is already NOT skippable
//			//Lemma5
//			for(auto const v: was_seen){
//				auto boundXV = customAPSP.distUpperBoundAdvanced(x,v);
//				auto boundYV = customAPSP.distUpperBoundAdvanced(y,v);
//				if((2*boundXV >= hdiffPlusOne) and (2*boundYV >= hdiffPlusOne)){
//					if((2*(largestDistToMates[v] - std::max(boundXV, boundYV)) >= condacc1)
//							or (2*(largestDistToMates[v] - std::max(exact_dist(x,v), exact_dist(y,v))) >= condacc1)){
//						if((2*largestDistToMates[v] - boundXV - boundYV >= condacc2)
//								or (2*largestDistToMates[v] - exact_dist(x,v) - exact_dist(y,v) >= condacc2)){
//							//node v is acceptable!
//							acceptable[v] = true;
//							//check if valuable
//							if(condval + 2*distC[v] > boundXV + boundYV
//									or condval + 2*distC[v] > exact_dist(x,v) + exact_dist(y,v)){
//								valuable[numberOfValNodes] = v;
//								numberOfValNodes++;
//							}
//						}
//					}
//				}
//			}
//		}else{
//			for(auto const v: was_seen){
//				auto boundXV = customAPSP.distUpperBoundAdvanced(x,v);
//				auto boundYV = customAPSP.distUpperBoundAdvanced(y,v);
//				//Corollary 7 with bounds and exact distances
//				if((2*boundXV >= hdiffPlusOne) and (2*boundYV >= hdiffPlusOne)){
//					//Lemma 8 with bounds and exact distances
//					if(2*largestDistToMates[v] - boundXV - boundYV >= condacc2
//							or 2*largestDistToMates[v] - exact_dist(x,v) - exact_dist(y,v) >= condacc2){
//						acceptable[v] = true;
//						//check if valuable
//						if(condval + 2*distC[v] > boundXV + boundYV
//								or condval + 2*distC[v] > exact_dist(x,v) + exact_dist(y,v) ){
//							valuable[numberOfValNodes] = v;
//							numberOfValNodes++;
//						}
//					}
//				}
//			}
//		}

		double endAccVal = get_wall_time();
		accValTime += (endAccVal - startAccVal);

		//----------------------------------------------------------------------------------------

		for(auto i=0; i < numberOfValNodes; ++i){
			const node v = valuable[i];
			for(const node w: mate[v]){
				if(acceptable[w]){
					double S2 = exact_dist(x,v) + exact_dist(y,w);
					double S3 = exact_dist(x,w) + exact_dist(y,v);
					double largest_diff = distXY + exact_dist(v,w) - std::max(S2, S3);
					quadtruples_visited++;
					if(h_diff < largest_diff){
						best_x = x;
						best_y = y;
						best_v = v;
						best_w = w;
						position_of_best = iteration;
						h_diff = largest_diff;
					}

					if(tuple->distance() <= h_diff){
						goto Ende;
					}
				}
			}
		}


		numberOfValNodes=0;
		valuable.clear();
		valuable.resize(n);
		acceptable.clear();
		acceptable.resize(n);
	}

	Ende:
	duration_main = (std::clock() - start_main) / (double) CLOCKS_PER_SEC;
	if(G.numberOfNodes() >= 1000){
	INFO("Running time of main part:", duration_main);
	INFO("Number of proper queries: ", numProperQueries);
	INFO("Number of total queries: ", numTotalQueries);
	INFO("Total running time of all queries (Akiba): ", sumQueryTimes);
	INFO("AccVal time: ", accValTime);
	INFO("Abortion Position (iteration): ", iteration);
	INFO("First quadtrupel attaining maximum:", "(",best_x,",",best_y,",",best_v, ",",best_w, ") with iteration position of (x,y) at ", position_of_best, "/", ordered_tuples.size(),
			" and distances d(x,y)= ", exact_dist(best_x,best_y), ", d(v,w)= ", exact_dist(best_v,best_w),
			", d(x,v)= ",exact_dist(best_x,best_v), ", d(y,w)= ", exact_dist(best_y, best_w),
			", d(x,w)= ",exact_dist(best_x,best_w), ", d(y,v)= ", exact_dist(best_y,best_v));
	INFO("Tuples visited: ", quadtruples_visited);
	INFO("h_diff: ", h_diff);
	INFO("Hyperbolicity: ", h_diff/2);

	double endALL = get_wall_time();
	INFO("Time of complete computation: ", endALL-startALL);
	}
	hyperbolicity_value = h_diff/2;
}



void Hyperbolicity::HYP(){
	auto get_wall_time = []()->double{
		struct timeval time;
		if(gettimeofday(&time, NULL)){
			return 0;
		}

		return (double) time.tv_sec + (double) time.tv_usec * .000001;
	};


	double startALL = get_wall_time();

	/*******************************************************************************/
	/* Preprocessing                                                             */
	/*******************************************************************************/

	double startX = get_wall_time();
	CustomizedAPSP customAPSP(G);
	customAPSP.run();
	double endX = get_wall_time();
	INFO("Time of CustomizedAPSP: ", endX-startX);

	SymMatrix<bool, node>  const& relevantPairs = customAPSP.getRelevantPairs();
	SymMatrix<edgeweight, node> const& distances = customAPSP.getExactDistances();
	std::vector<edgeweight> const& farness = customAPSP.getFarnessVec();
	std::vector<edgeweight> const& ecc = customAPSP.getEccentricityVec();



	/******************************All types of comparisons ********************************/

	auto lex_comp = [](NodeTupleWithDist const& firstTuple, NodeTupleWithDist const& secondTuple)-> bool {
		int res = firstTuple.distance() - secondTuple.distance();
		//If same distance then consider (increasing) lexicographic order
		if(res == 0){
	        int diff = (firstTuple.index_a() - secondTuple.index_a());
	        if (diff == 0)
	        {
	            diff =  firstTuple.index_b() - secondTuple.index_b();
	        }
	        // (x,y) is larger than (v,w) if x < v or (x=v and y < w), and smaller otherwise
	        return diff < 0;
		}
		return res > 0;
	};

	auto farness_comp = [&](NodeTupleWithDist const& firstTuple, NodeTupleWithDist const& secondTuple)-> bool {
		int res = firstTuple.distance() - secondTuple.distance();
		if(res == 0){
			int diff = farness[firstTuple.index_a()] + farness[firstTuple.index_b()] - farness[secondTuple.index_a()] - farness[secondTuple.index_b()];
			return diff <= 0;
		}

		return res > 0;
	};


	auto degree_comp =[&](NodeTupleWithDist const& firstTuple, NodeTupleWithDist const& secondTuple)-> bool {
		int res = firstTuple.distance() - secondTuple.distance();
		if(res == 0){
			int diff = G.degree(firstTuple.index_a()) + G.degree(firstTuple.index_b()) - G.degree(secondTuple.index_a()) - G.degree(secondTuple.index_b());
			return diff >= 0;
		}
		return res > 0;
	};

	/************************************************************************************/


	//Note: farness_comp causes additional time overhead in the tuple construction phase
	//		since it conducts additional comparisons if distances are the same.
	//auto ordered_tuples = std::set<NodeTupleWithDist, decltype(lex_comp)>(lex_comp);
	std::vector<NodeTupleWithDist> ordered_tuples;
	//std::vector<NodeTupleWithDist> ordered_tuples;

	//Tuple construction
	double startTuples = get_wall_time();
	//Add all far-apart pairs to the set
	for(node u=0; u < G.numberOfNodes(); ++u){
		for(node v = u; v < G.numberOfNodes(); ++v){
			if(relevantPairs.element(u,v)){
				ordered_tuples.emplace_back(u,v, distances.element(u,v));
			}
		}
	}
	std::sort(ordered_tuples.begin(), ordered_tuples.end(), lex_comp);

	//std::sort(ordered_tuples.begin(), ordered_tuples.end(), farness_comp);

	double endTuples = get_wall_time();
	INFO("Time of Tuple Construction: ", endTuples-startTuples);

	INFO("Diameter: ", ordered_tuples.begin()->distance());
	INFO("Number of Far Apart Pairs: ", ordered_tuples.size(),
			" with percantage of ", static_cast<double>(ordered_tuples.size()/static_cast<double>((G.numberOfNodes()*(G.numberOfNodes()-1)/2)))*100,"%");

	// central node via TopCloseness
	double start_c = get_wall_time();
	TopCloseness topCentral(G, 1);
	topCentral.run();
	node central_node = topCentral.topkNodesList()[0];
	double end_c = get_wall_time();
	INFO("Time of Central Node Computation (TopCloseness): ", end_c - start_c);

	//central node via usual algorithm
//	double start_c = get_wall_time();
//	node central_node = centralNode(distances);
//	double end_c = get_wall_time();
//	INFO("Time of Central Node Computation (Naive): ", end_c - start_c);


//------------------------------------------MAIN PART ----------------------------------------------------

	edgeweight h_diff = 0;
	std::set<node> was_seen;
	std::vector<bool> wasSeen(G.upperNodeIdBound(), false);
	std::vector<std::vector<node>> mate(G.upperNodeIdBound());
	std::vector<edgeweight> largestDistToMates(G.upperNodeIdBound());

	double startMAIN = get_wall_time();

	size_t iteration = 0;
	node best_x=0;
	node best_y=0;
	node best_v=0;
	node best_w=0;
	auto degree_x=0;
	auto degree_y=0;
	auto degree_v=0;
	auto degree_w=0;
	count quadtruples_visited=0;
	edgeweight eccL_xy = 0;
	edgeweight eccL_vw = 0;
	edgeweight distByCompOfBestxy = 0;
	edgeweight distByCompOfBestvw = 0;
	edgeweight distEugOfBestxy = 0;
	edgeweight distEugOfBestvw = 0;

	size_t position_of_best=0;

	for(auto tuple = ordered_tuples.begin(); tuple != ordered_tuples.end(); ++tuple){
		if(tuple->distance() <= h_diff){
			goto Ende;
		}

		//INFO(it_1->distance(),"  (", it_1->index_a(), ",", it_1->index_b(), ")");

		iteration++;
		node x = tuple->index_a();
		node y = tuple->index_b();
		//Lemma5
		//TODO hier clevere Sortierung wie in hyp.c von Borassi!
		for(auto const v: was_seen){
//		for(auto rit = was_seen.rbegin(); rit != was_seen.rend(); ++rit){
//			const node v = *rit;
			if(is_valuable(distances,x,y,v, h_diff, largestDistToMates[v], central_node) ){
				for(const node w: mate[v]){
					if(is_acceptable(distances, x,y,w, h_diff, largestDistToMates[w])){
						edgeweight S1 = tuple->distance() + distances.element(v,w);
						edgeweight S2 = distances.element(x,v) + distances.element(y,w);
						edgeweight S3 = distances.element(x,w) + distances.element(y,v);
						edgeweight largest_diff = std::max(std::max(S1, S2), S3) - std::max(std::min(S1,S2),S3);
						//double largest_diff = S1 - std::max(S2,S3);
						quadtruples_visited++;
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
					if(tuple->distance() <= h_diff){
						eccL_xy = std::max(customAPSP.eccLowerBound(x), customAPSP.eccLowerBound(y));
						eccL_vw = std::max(customAPSP.eccLowerBound(v), customAPSP.eccLowerBound(w));
						distByCompOfBestxy = tuple->distance();
						distByCompOfBestvw = customAPSP.distUpperBoundByComp(v,w);
						distEugOfBestxy = customAPSP.distUpperBound(x,y);
						distEugOfBestvw = customAPSP.distUpperBound(v,w);
						goto Ende;
					}
				}
			}
		}

		mate[x].push_back(y);
		mate[y].push_back(x);
		if(largestDistToMates[x] < distances.element(x,y)){
			largestDistToMates[x] = distances.element(x,y);
		}
		if(largestDistToMates[y] < distances.element(x,y)){
			largestDistToMates[y] = distances.element(x,y);
		}
		was_seen.insert(x);
		was_seen.insert(y);
	}

	Ende:
	double endMAIN = get_wall_time();
	INFO("Laufzeit vom Hauptblock:", endMAIN-startMAIN);
	INFO("Number of pairs considered (Abortion Position): ", iteration);
	INFO("First quadtrupel attaining maximum:", "(",best_x,",",best_y,",",best_v, ",",best_w, ") with iteration position of (x,y) at ", position_of_best, "/", ordered_tuples.size(),
			" and distances d(x,y)= ", distances.element(best_x,best_y), ", d(v,w)= ", distances.element(best_v,best_w),
			", d(x,v)= ",distances.element(best_x,best_v), ", d(y,w)= ", distances.element(best_y, best_w),
			", d(x,w)= ",distances.element(best_x,best_w), ", d(y,v)= ", distances.element(best_y,best_v));
	INFO("Quadtruples visited: ", quadtruples_visited);
	INFO("Degrees: ", degree_x, " ", degree_y, " ", degree_v, " ", degree_w);
	INFO("h_diff: ", h_diff);
	INFO("eccLower(x,y) = ", eccL_xy, " distComp(x,y) = ", distByCompOfBestxy, " distEug(x,y) = ", distEugOfBestxy);
	INFO("eccLower(v,w) = ", eccL_vw, " distComp(v,w) = ", distByCompOfBestvw, " distEu(v,w) = ", distEugOfBestvw);

	INFO("Hyperbolicity: ", h_diff/2);
	double endALL = get_wall_time();
	INFO("Time of complete computation: ", endALL-startALL);
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
	//3delta_L + 3/2 >= d(x,y)
//	bool acceptable = false;
//
//    if (3 * (current_max_diff+1) >= 2 * distXY) {
//    	auto condacc1 = 3 * (current_max_diff+1) - 2 * distXY;
//    	auto condacc2 = 2 * (current_max_diff+1) - distXY;
//    	if (2 * (eccentricity - distXV) >= condacc1) {
//    		if (2 * (eccentricity - distYV) >= condacc1) {
//    			if(2*distXV >= current_max_diff+1 && 2 * distYV >= current_max_diff+1) {
//    				if (2 * eccentricity >= condacc2 + distXV + distYV ) {
//    					acceptable = true;
//    				}
//    			}
//    		}
//    	}
//    }else{
//    	if(2 * distXV >= current_max_diff+1 && 2 * distYV >= current_max_diff+1) {
//    		if (2 * eccentricity + distXY >= 2 * (current_max_diff+1) + distXV + distYV ) {
//    			acceptable = true;
//    		}
//    	}
//    }
//
//    return acceptable;
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

//        if (2 * distVC + distXY - current_max_diff > distXV + distYV) {
//            return true;
//        }
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
	    for(node v = 0; v < G.numberOfNodes(); ++v){
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
	    edgeweight s1 = it_1->distance() + it_2->distance(); // d(a,b) + d(c,d)
	    edgeweight s2 = allPairShortestPaths.getDistance(it_1->index_a(), it_2->index_a()) +  // d(a,c) + d(b,d)
	    allPairShortestPaths.getDistance(it_1->index_b(), it_2->index_b());
	    edgeweight s3 = allPairShortestPaths.getDistance(it_1->index_a(),it_2->index_b()) +   // d(a,d) + d(b,c)
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
