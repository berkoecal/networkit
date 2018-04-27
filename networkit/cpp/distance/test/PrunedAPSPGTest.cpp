/*
 * PrunedAPSPGTest.cpp
 *       Created on: 01.12.2017
 *       Author: Berk Öcal
 */



#ifndef NOGTEST

#include <cstdio>
#include <ctime>
#include "PrunedAPSPGTest.h"
#include "../Bisenius/PrunedAPSP.h"
#include "../PruningAPSP.h"
#include "../APSP.h"
#include "../../hyperbolicity/CustomizedAPSP.h"
#include "../../hyperbolicity/SymmetricMatrix.h"
#include "../../centrality/DegreeCentrality.h"
#include "../../io/METISGraphReader.h"
#include "../../io/SNAPGraphReader.h"

namespace NetworKit {
  
TEST_F(PrunedAPSPGTest, testPrunedAPSP) {
    
	//Graph of Paper Akiba et al.
	int n = 12;
	Graph G(n);
	G.addEdge(0,7);
	G.addEdge(0,9);
	G.addEdge(0,11);
	G.addEdge(0,8);
	G.addEdge(0,6);
	G.addEdge(0,5);
	G.addEdge(0,3);
	G.addEdge(1,7);
	G.addEdge(1,4);
	G.addEdge(1,9);
	G.addEdge(2,11);
	G.addEdge(2,10);
	G.addEdge(3,4);
	G.addEdge(3,5);
	G.addEdge(4,7);
	G.addEdge(5,6);
	G.addEdge(6,8);
	G.addEdge(8,11);
	G.addEdge(9,10);
	
	APSP usual(G);
	usual.run();
	DEBUG("APSP-Result for (6,11): ", usual.getDistance(5,10));
	DEBUG("APSP-Result for (2,11): ", usual.getDistance(1,10));
	DEBUG("APSP-Result for (8,2): ", usual.getDistance(7,2));
	DEBUG("APSP-Result for (4,9): ", usual.getDistance(3,8));	
	PrunedAPSP pruned(G);
	pruned.run();
	DEBUG("APSPPruned-Result for (6,11): ", pruned.getDistance(5,10));
	DEBUG("APSPPruned-Result for (2,11): ", pruned.getDistance(1,10));
	DEBUG("APSPPruned-Result for (8,2): ", pruned.getDistance(7,2));
	DEBUG("APSPPruned-Result for (4,9): ", pruned.getDistance(3,8));
	
	EXPECT_TRUE(pruned.getDistance(3,7) == usual.getDistance(3,7));
}

TEST_F(PrunedAPSPGTest, testPrunedAPSPFarApart) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input/snap/ca-GrQc.txt");

//	Graph G(8);
//	G.addEdge(0,1);
//	G.addEdge(0,2);
//	G.addEdge(1,3);
//	G.addEdge(1,4);
//	G.addEdge(2,5);
//	G.addEdge(3,6);
//	G.addEdge(3,7);

	INFO("n= ", G.numberOfNodes());
	INFO("m= ", G.numberOfEdges());

}

TEST_F(PrunedAPSPGTest, testPrunedAPSP2) {

	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/power.graph");

	INFO("n= ", G.numberOfNodes());
	INFO("m= ", G.numberOfEdges());

	std::clock_t start3;
	double duration3;
	start3 = std::clock();

	PrunedAPSP pruning(G);
	pruning.run();

	duration3 = (std::clock() - start3) / (double) CLOCKS_PER_SEC;

	INFO("TIME Pruning: ", duration3);

	std::clock_t start;
	double duration;
	start = std::clock();

	CustomizedAPSP cAPSP(G);
	cAPSP.run();

	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;

	INFO("Time CustomizedAPSP: ", duration);

	auto& distances = cAPSP.getExactDistances();


//	int errors=0;
//	double errorSum=0.0;
//	double relativeError=0.0;

//	std::clock_t start2;
//	double duration_cAPSP;
//	start2 = std::clock();
//	CustomizedAPSP cAPSP2(G, 20);
//	cAPSP2.run();
//	duration_cAPSP = (std::clock() - start2) / (double) CLOCKS_PER_SEC;
//	INFO("Duration CustomizedAPSP Degree20: ", duration_cAPSP);


//	for(node u= 0; u < G.numberOfNodes(); ++u){
//		for(node v=u; v < G.numberOfNodes(); ++v){
//			if((cAPSP.getBoundedDistance(u,v) != cAPSP2.getBoundedDistance(u,v))
//					and cAPSP2.getBoundedDistance(u,v) != std::numeric_limits<edgeweight>::max()){
//				errors++;
//				errorSum +=(cAPSP2.getBoundedDistance(u,v) - cAPSP.getBoundedDistance(u,v));
//				relativeError+= (cAPSP2.getBoundedDistance(u,v) - cAPSP.getBoundedDistance(u,v))/cAPSP.getBoundedDistance(u,v);
//			}
//		}
//	}
//
//	INFO("Errors:", errors);
//	INFO("Error Average:", errorSum/(double)errors);
//	INFO("Relative Error Average: ", relativeError/(double)errors);


	//QueryTime
//	std::clock_t start_query;
//	double duration_query;
//	start_query = std::clock();
//	for(node u=0; u < G.numberOfNodes(); ++u){
//		for(node v=u+1; v < G.numberOfNodes(); ++v){
//			pruning.getDistance(u,v);
//		}
//	}
//	duration_query = (std::clock() - start_query) / (double) CLOCKS_PER_SEC;
//	INFO("Time Pruning Query: ", duration_query);


}

TEST_F(PrunedAPSPGTest, testPrunedAPSP3) {
  
	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/astro-ph.graph");
	
	std::clock_t start2;
	double duration2;
	start2 = std::clock();
	
	APSP usual(G);
	usual.run();	
	duration2 = (std::clock() - start2) / (double) CLOCKS_PER_SEC;
	
	INFO("Time of APSP: ", duration2);
	
	std::clock_t start;
	double duration;	
	start = std::clock();
	
	PruningAPSP pruned(G);
	pruned.runParallel();
	INFO("PrunedAPSP (parallel) has run! Next compute distances...");
	for(node u=0; u < G.upperNodeIdBound(); ++u){
	  for(node v=u+1; v < G.upperNodeIdBound(); ++v){
	    EXPECT_TRUE(pruned.getDistance(u,v)==usual.getDistance(u,v));
	  }
	}
	
	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	
	INFO(" Time of PrunedAPSP(parallel): ", duration);
  
 	
}

TEST_F(PrunedAPSPGTest, testPrunedAPSP4) {
  auto Reader = METISGraphReader();
  Graph G = Reader.read("input/astro-ph.graph");
  
  std::clock_t start2;
  double duration2;
  start2 = std::clock();
  
  APSP usual(G);
  usual.run();
  
  duration2 = (std::clock() - start2) / (double) CLOCKS_PER_SEC;
  
  INFO("Time of APSP: ", duration2);
  
  
  for(node const u: G.nodes()){
    for(node const v: G.nodes()){
      EXPECT_TRUE(usual.getDistance(u,v) == usual.getDistance(u,v));
    }
  }
}
} /* namespace NetworKit */

#endif /*NOGTEST */
