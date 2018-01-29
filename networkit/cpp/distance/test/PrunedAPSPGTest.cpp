/*
 * PrunedAPSPGTest.cpp
 *       Created on: 01.12.2017
 *       Author: Berk Öcal
 */



#ifndef NOGTEST

#include <cstdio>
#include <ctime>
#include "PrunedAPSPGTest.h"
#include "../PrunedAPSP.h"
#include "../APSP.h"
#include "../../io/METISGraphReader.h"

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

TEST_F(PrunedAPSPGTest, testPrunedAPSP2) {

	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/power.graph");
	
	std::clock_t start2;
	double duration2;
	start2 = std::clock();
	
	APSP usual(G);
	usual.run();
	
	duration2 = (std::clock() - start2) / (double) CLOCKS_PER_SEC;
	
	DEBUG("Time of APSP: ", duration2);
	
	std::clock_t start;
	double duration;	
	start = std::clock();
	
	PrunedAPSP pruned(G);
	pruned.run();
	
	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	
	DEBUG(" Time of PrunedAPSP: ", duration);
	
// 	std::clock_t startParallel;
// 	double durationParallel;	
// 	startParallel = std::clock();
// 	
// 	PrunedAPSP prunedParallel(G);
// 	prunedParallel.runParallel();
// 	
// 	durationParallel = (std::clock() - startParallel) / (double) CLOCKS_PER_SEC;
// 	
// 	DEBUG(" Time of PrunedAPSP parallel: ", durationParallel);
	
	for(node const u: G.nodes()){
	  for(node const v: G.nodes()){
	    EXPECT_TRUE(usual.getDistance(u,v) == pruned.getDistance(u,v));
	  }
	}	
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
	
	PrunedAPSP pruned(G);
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