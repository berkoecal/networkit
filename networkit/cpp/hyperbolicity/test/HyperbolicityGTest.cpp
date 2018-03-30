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
#include "../../io/SNAPGraphReader.h"
#include "../CustomizedBFS.h"
#include "../../distance/APSP.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../centrality/TopCloseness.h"
#include "../../distance/BFS.h"
#include "../../components/ConnectedComponents.h"
#include "../../components/BiconnectedComponents.h"

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
  
  TEST_F(HyperbolicityGTest, testHyperbolicityMETIS2){
		auto Reader = METISGraphReader();
		Graph G = Reader.read("input/power.graph");

		Hyperbolicity hyperbolicity(G);
		hyperbolicity.run();

  }

  TEST_F(HyperbolicityGTest, testHyperbolicityMETIS) {

	//Graph g = ErdosRenyiGenerator(500, 0.1, false).generate();

	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/power.graph");

	BiconnectedComponents bcc(G);
	bcc.run();

	INFO("Number of biconnected comps: ", bcc.numberOfComponents());

	//compute hyperbolicity for each biconnected component and save maximum value among them
	double maximum_seen = 0.0;
	for(count i=0; i < bcc.numberOfComponents(); ++i){
		auto ith_comp = bcc.getComponents()[i];
		Graph component;

		std::vector<node> map(G.numberOfNodes());
		for(node u: ith_comp){
			node newNode = component.addNode();
			map[u]=newNode;
		}

		auto edges = G.edges();

		for(auto edge : edges){
			node first = edge.first;
			node second = edge.second;

			if(std::find(ith_comp.begin(), ith_comp.end(), first) != ith_comp.end() and
					std::find(ith_comp.begin(), ith_comp.end(), second) != ith_comp.end()){
				component.addEdge(map[first], map[second]);
			}
		}

		INFO("Anzahl der Knoten und Kanten der Komponente: ", component.numberOfNodes(), " und ", component.numberOfEdges());

		Hyperbolicity hyperbolicity(component);
		hyperbolicity.run();
		double value = hyperbolicity.getHyperbolicity();

		if(value > maximum_seen){
			maximum_seen = value;
		}

	}

	INFO("Final Hyperbolicity over all connected components: ", maximum_seen);


//  	ConnectedComponents cc(G);
//  	cc.run();
//
//  	INFO("Number of components: ", cc.numberOfComponents());
//
//    double maximum_seen = 0.0;
//    for(count i=0; i < cc.numberOfComponents(); ++i){
//    	auto component = cc.getComponents()[i];
//    	Graph temp = G;
//
//    	for(node u=0; u < G.numberOfNodes(); ++u){
//    		if(cc.componentOfNode(u) != i){
//    			temp.removeNode(u);
//    		}
//    	}
//
//    	Graph graphComponent;
//
//    	//auto edges = temp.edges();
//    	std::vector<node> map(G.numberOfNodes());
//    	for(node u: temp.nodes()){
//    		node newNode = graphComponent.addNode();
//    		map[u]=newNode;
//    	}
//
//    	auto edges = temp.edges();
//
//    	for(auto edge : edges){
//    		node first = edge.first;
//    		node second = edge.second;
//    		graphComponent.addEdge(map[first], map[second]);
//    	}
//
//    	INFO("Anzahl der Knoten und Kanten der Komponente: ", graphComponent.numberOfNodes(), " und ", graphComponent.numberOfEdges());
//
//    	Hyperbolicity hyperbolicity(graphComponent);
//        hyperbolicity.run();
//    	auto value = hyperbolicity.getHyperbolicity();
//
//    	if(value > maximum_seen){
//    		maximum_seen = value;
//    	}
//  }
//
// 		INFO("Final Hyperbolicity over all connected components: ", maximum_seen);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityRandom) {

	Graph G = ErdosRenyiGenerator(500, 0.5, false).generate();

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

  TEST_F(HyperbolicityGTest, testHyperbolicitySNAP) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input/snap/CA-GrQc.txt");

	  BiconnectedComponents bcc(G);
	  bcc.run();

	  INFO("Number of biconnected comps: ", bcc.numberOfComponents());

	  //compute hyperbolicity for each biconnected component and save maximum value among them
	  double maximum_seen = 0.0;
	  for(count i=0; i < bcc.numberOfComponents(); ++i){
		  auto ith_comp = bcc.getComponents()[i];
		  Graph component;

		  std::vector<node> map(G.numberOfNodes());
		  for(node u: ith_comp){
			  node newNode = component.addNode();
			  map[u]=newNode;
		  }

		  auto edges = G.edges();

		  for(auto edge : edges){
			  node first = edge.first;
			  node second = edge.second;

			  if(std::find(ith_comp.begin(), ith_comp.end(), first) != ith_comp.end() and
				std::find(ith_comp.begin(), ith_comp.end(), second) != ith_comp.end()){
				  component.addEdge(map[first], map[second]);
			  }
		  }

		  INFO("Anzahl der Knoten und Kanten der Komponente: ", component.numberOfNodes(), " und ", component.numberOfEdges());

		  Hyperbolicity hyperbolicity(component);
		  hyperbolicity.run();
		  double value = hyperbolicity.getHyperbolicity();

		  if(value > maximum_seen){
			  maximum_seen = value;
		  }

	  }

	  INFO("Final Hyperbolicity over all connected components: ", maximum_seen);

//	  ConnectedComponents cc(G);
//	  cc.run();
//
//	  INFO("Number of components: ", cc.numberOfComponents());
//
//	  double maximum_seen = 0.0;
//	  for(count i=0; i < cc.numberOfComponents(); ++i){
//		  auto component = cc.getComponents()[i];
//		  Graph temp = G;
//
//		  for(node u=0; u < G.numberOfNodes(); ++u){
//			  if(cc.componentOfNode(u) != i){
//				  temp.removeNode(u);
//			  }
//		  }
//
//		  Graph graphComponent;
//
//		  //auto edges = temp.edges();
//		  std::vector<node> map(G.numberOfNodes());
//		  for(node u: temp.nodes()){
//			  node newNode = graphComponent.addNode();
//			  map[u]=newNode;
//		  }
//
//		  auto edges = temp.edges();
//
//		  for(auto edge : edges){
//			  node first = edge.first;
//			  node second = edge.second;
//			  graphComponent.addEdge(map[first], map[second]);
//		  }
//
//		  INFO("Anzahl der Knoten und Kanten der Komponente: ", graphComponent.numberOfNodes(), " und ", graphComponent.numberOfEdges());
//
//		  Hyperbolicity hyperbolicity(graphComponent);
//		  hyperbolicity.run();
//		  double value = hyperbolicity.getHyperbolicity();
//
//		  if(value > maximum_seen){
//			  maximum_seen = value;
//		  }
//	  }
//
//	  INFO("Final Hyperbolicity over all connected components: ", maximum_seen);

  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDistCheck) {

	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/fe_4elt2.graph");

	std::clock_t start;
	double duration_apsp;
	start = std::clock();
	APSP _apsp(G);
	_apsp.run();
	duration_apsp = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	INFO("Duration of APSP ", duration_apsp);

	TopCloseness topk(G, 20);
	topk.run();

	std::clock_t start2;
	double duration_cBFS;
	start2 = std::clock();
	CustomizedBFS cBFS(G);
	cBFS.run();
	duration_cBFS = (std::clock() - start2) / (double) CLOCKS_PER_SEC;
	INFO("Duration CustomizedBFS: ", duration_cBFS);

	int errors=0;
	int errorSum=0;

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cBFS.getBoundedDistance(u,v)){
				errors++;
				errorSum +=(cBFS.getBoundedDistance(u,v) - _apsp.getDistance(u,v));
			}
		}
	}

	INFO("Errors:", errors);
	INFO("Error Average:", errorSum/errors);

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
	  auto Reader = METISGraphReader();
	  Graph G = Reader.read("input/wing.graph");

	  //closeness
	  std::clock_t start1;
	  double duration_closeness;
	  start1 = std::clock();

	  TopCloseness topk(G, 100);
	  topk.run();

	  duration_closeness = (std::clock() - start1) / (double) CLOCKS_PER_SEC;
	  INFO("Duration of Closeness:", duration_closeness);

	  //customizedBFS for all nodes
	  std::clock_t start2;
	  double duration_cBFS;
	  start2 = std::clock();
	  CustomizedBFS cBFS(G);
	  cBFS.run();
	  duration_cBFS = (std::clock() - start2) / (double) CLOCKS_PER_SEC;
	  INFO("Duration of CustomizedBFS:", duration_cBFS);

	  //customizedBFS for 100 nodes
	  std::clock_t start3;
	  double duration_cBFS_hundret;
	  start3 = std::clock();
	  CustomizedBFS cBFSHundret(G, topk.topkNodesList());
	  cBFSHundret.run();
	  duration_cBFS_hundret = (std::clock() - start3) / (double) CLOCKS_PER_SEC;
	  INFO("Duration of CustomizedBFSHundret:", duration_cBFS_hundret);

	  //APSP
	  std::clock_t start4;
	  double duration_apsp;
	  start4 = std::clock();
	  APSP apsp(G);
	  apsp.run();
	  duration_apsp = (std::clock() - start4) / (double) CLOCKS_PER_SEC;
	  INFO("Duration of APSP:", duration_apsp);

   }

  TEST_F(HyperbolicityGTest, testCustomizedBFS2){
	 int n = 3;
	 Graph G(n);
	 G.addEdge(0,1);
	 G.addEdge(0,2);

	 CustomizedBFS cBFS(G, G.nodes());
	 cBFS.run();

	 SymMatrix<bool, node> sym_matrix = cBFS.getRelevantPairs();

	 for(size_t i=0; i < sym_matrix.size(); ++i){
		 for(size_t j=0; j < sym_matrix.size(); ++j){
			 INFO("(",i ,",", j, "):", sym_matrix.element(i,j));
		 }
	 }
  }
}
