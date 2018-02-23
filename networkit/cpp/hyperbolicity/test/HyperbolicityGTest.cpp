/*
 * HyperbolicityGTest.cpp
 *       Created on: 18.11.2017
 *       Author: Berk Öcal
 */

#include <cstdio>
#include <ctime>
#include "HyperbolicityGTest.h"
#include "../Hyperbolicity.h"
#include "../../io/METISGraphReader.h"
#include "../../io/KONECTGraphReader.h"
#include "../CustomizedBFS.h"
#include "../../distance/APSP.h"

namespace NetworKit{
  TEST_F(HyperbolicityGTest, testHyperbolicity) {
 
    int n = 4;
    Graph G(n);
    
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(3, 0);
    
    Hyperbolicity hyperbolicity(G);
    hyperbolicity.run(); 
  }
  
  TEST_F(HyperbolicityGTest, testHyperbolicityMETIS) {

	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/mapping/grid-8x8x8-dist-arch.graph");

	std::clock_t start;
	double duration;
	start = std::clock();

    Hyperbolicity hyperbolicity(G);
    hyperbolicity.run();

    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityKONECT) {

	auto Reader = KONECTGraphReader();
	Graph G = Reader.read("input/foodweb-baydry.konect");

	std::clock_t start;
	double duration;
	start = std::clock();

    Hyperbolicity hyperbolicity(G);
    hyperbolicity.run();

    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDistCheck) {

	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/power.graph");

	APSP _apsp(G);
	_apsp.run();

	CustomizedBFS cBFS(G);
	cBFS.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			EXPECT_TRUE(_apsp.getDistance(u,v) == cBFS.getDistance(u,v));
		}
	}

  }

  TEST_F(HyperbolicityGTest, testHyperbolicity2) {
    
	//4x4 grid

    int n = 16;
    Graph G(n);
    
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(6, 7);
    G.addEdge(8, 9);
    G.addEdge(9, 10);
    G.addEdge(10, 11);
    G.addEdge(12, 13);
    G.addEdge(13, 14);
    G.addEdge(14, 15);
    G.addEdge(0, 4);
    G.addEdge(1, 5);
    G.addEdge(2, 6);
    G.addEdge(3, 7);
    G.addEdge(4, 8);
    G.addEdge(5, 9);
    G.addEdge(6, 10);
    G.addEdge(7, 11);
    G.addEdge(8, 12);
    G.addEdge(9, 13);
    G.addEdge(10, 14);
    G.addEdge(11, 15);

    Hyperbolicity hyperbolicity(G);
    hyperbolicity.run(); 
  }
  


  TEST_F(HyperbolicityGTest, testHyperbolicity3) {
		 int n = 11;
		 Graph G(n);
		 G.addEdge(0,1);
		 G.addEdge(0,2);
		 G.addEdge(1,3);
		 G.addEdge(1,4);
		 G.addEdge(2,5);
		 G.addEdge(1,3);
		 G.addEdge(4,6);
		 G.addEdge(4,7);
		 G.addEdge(5,8);
		 G.addEdge(8,9);
		 G.addEdge(9,10);

    std::clock_t start;
    double duration;
    start = std::clock();

    Hyperbolicity hyperbolicity(G);
    hyperbolicity.run();

    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;

    INFO("Duration:", duration);
  }

  TEST_F(HyperbolicityGTest, testCustomizedBFS){
	 int n = 3;
	 Graph G(n);
	 G.addEdge(0,1);
	 G.addEdge(1,2);
	 G.addEdge(2,0);

	 CustomizedBFS cBFS(G);
	 cBFS.run();

	 SymMatrix<bool, node> sym_matrix = cBFS.getFarApartPairs();

	 for(size_t i=0; i < sym_matrix.size(); ++i){
		 for(size_t j=0; j < sym_matrix.size(); ++j){
			 INFO("(",i ,",", j, "):", sym_matrix.element(i,j));
		 }
	 }
  }

  TEST_F(HyperbolicityGTest, testCustomizedBFS2){
	 int n = 3;
	 Graph G(n);
	 G.addEdge(0,1);
	 G.addEdge(0,2);

	 CustomizedBFS cBFS(G);
	 cBFS.run();

	 SymMatrix<bool, node> sym_matrix = cBFS.getFarApartPairs();

	 for(size_t i=0; i < sym_matrix.size(); ++i){
		 for(size_t j=0; j < sym_matrix.size(); ++j){
			 INFO("(",i ,",", j, "):", sym_matrix.element(i,j));
		 }
	 }
  }
}
