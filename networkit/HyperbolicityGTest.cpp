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
#include "../../distance/APSP.h"
#include "../../generators/ErdosRenyiGenerator.h"
#include "../../centrality/TopCloseness.h"
#include "../../centrality/ApproxBetweenness.h"
#include "../../centrality/DegreeCentrality.h"
#include "../../distance/BFS.h"
#include "../../components/ConnectedComponents.h"
#include "../../components/BiconnectedComponents.h"
#include "../CustomizedAPSP.h"
//#include "../../distance/Akiba/src/pruned_landmark_labeling.h"

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
	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/4elt.graph");

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
	  Graph G = Reader.read("input_borassi/032233cit-HepPh");

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

		  if(component.numberOfNodes() > 4){
			  INFO("Anzahl der Knoten und Kanten der Komponente: ", component.numberOfNodes(), " und ", component.numberOfEdges());
		  }

		  Hyperbolicity hyperbolicity(component);
		  hyperbolicity.run();
		  double value = hyperbolicity.getHyperbolicity();

		  if(value > maximum_seen){
			  maximum_seen = value;
		  }

	  }

	  INFO("Final Hyperbolicity over all connected components: ", maximum_seen);

  }

  TEST_F(HyperbolicityGTest, testHyperbolicitySNAP2) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/024982cit-HepTh");

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

		  if(component.numberOfNodes() > 4){
			  INFO("Anzahl der Knoten und Kanten der Komponente: ", component.numberOfNodes(), " und ", component.numberOfEdges());
		  }

		  Hyperbolicity hyperbolicity(component);
		  hyperbolicity.run();
		  double value = hyperbolicity.getHyperbolicity();

		  if(value > maximum_seen){
			  maximum_seen = value;
		  }

	  }

	  INFO("Final Hyperbolicity over all connected components: ", maximum_seen);

  }

  TEST_F(HyperbolicityGTest, testHyperbolicitySNAP3) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/009577ca-AstroPh");

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

		  if(component.numberOfNodes() > 4){
			  INFO("Anzahl der Knoten und Kanten der Komponente: ", component.numberOfNodes(), " und ", component.numberOfEdges());
		  }

		  Hyperbolicity hyperbolicity(component);
		  hyperbolicity.run();
		  double value = hyperbolicity.getHyperbolicity();

		  if(value > maximum_seen){
			  maximum_seen = value;
		  }

	  }

	  INFO("Final Hyperbolicity over all connected components: ", maximum_seen);

  }

  TEST_F(HyperbolicityGTest, testHyperbolicitySNAP4) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/020095p2p-Gnutella30");

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

		  if(component.numberOfNodes() > 4){
			  INFO("Anzahl der Knoten und Kanten der Komponente: ", component.numberOfNodes(), " und ", component.numberOfEdges());
		  }

		  Hyperbolicity hyperbolicity(component);
		  hyperbolicity.run();
		  double value = hyperbolicity.getHyperbolicity();

		  if(value > maximum_seen){
			  maximum_seen = value;
		  }

	  }

	  INFO("Final Hyperbolicity over all connected components: ", maximum_seen);

  }
  TEST_F(HyperbolicityGTest, testHyperbolicityDistCheck) {

	auto Reader = METISGraphReader();
	Graph G = Reader.read("input/power.graph");

//	Graph G(8);
//	G.addEdge(0,1);
//	G.addEdge(0,2);
//	G.addEdge(1,3);
//	G.addEdge(1,4);
//	G.addEdge(2,5);
//	G.addEdge(3,6);
//	G.addEdge(3,7);


	int errors=0;
	int errorsForLargeDistances;
	double errorSum=0.0;
	double errorSumForLargeDistances=0.0;
	double relativeError=0.0;
	double relativeErrorForLargeDistances=0.0;

	std::clock_t start2;
	double duration_cAPSP;
	start2 = std::clock();
	//CustomizedAPSP cAPSP(G, topk.topkNodesList());
	count k = 1000;
	count maxComp = 20;
	CustomizedAPSP cAPSP(G, k, maxComp);
	cAPSP.run();
	duration_cAPSP = (std::clock() - start2) / (double) CLOCKS_PER_SEC;
	INFO("Duration CustomizedAPSP with ",k, " nodes: ", duration_cAPSP);

	std::clock_t start;
	double duration_apsp;
	start = std::clock();

	CustomizedAPSP _apsp(G);
	_apsp.run();
	duration_apsp = (std::clock() - start) / (double) CLOCKS_PER_SEC;
	INFO("Duration of CustomizedAPSP ", duration_apsp);

	count largeDistance = 30;

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.getDistance(u,v)){
				errors++;
				errorSum +=(cAPSP.getDistance(u,v) - _apsp.getDistance(u,v));
				relativeError+= (cAPSP.getDistance(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);

				if(_apsp.getDistance(u,v) >= largeDistance){
					errorsForLargeDistances++;
					errorSumForLargeDistances +=(cAPSP.getDistance(u,v) - _apsp.getDistance(u,v));
					relativeErrorForLargeDistances+= (cAPSP.getDistance(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
				}
			}
		}
	}

	INFO("Errors:", errors);
	INFO("Error Average:", errorSum/errors);
	INFO("Relative Error Average: ", relativeError/errors);
	INFO("Errors (large Distance):", errorsForLargeDistances);
	INFO("Error Average (large Distance):", errorSumForLargeDistances/errorsForLargeDistances);
	INFO("Relative Error Average (large Distance): ", relativeErrorForLargeDistances/errorsForLargeDistances);

  }

//  TEST_F(HyperbolicityGTest, testAKIBA) {
////	  auto Reader = SNAPGraphReader();
////	  Graph G = Reader.read("input_borassi/045425soc-Slashdot0811");
//	  PrunedLandmarkLabeling<> pll;
//	  pll.ConstructIndex(samples/graph_example.tsv);
//	  cout << pll.QueryDistance(1, 4) << endl;
//  }

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


  TEST_F(HyperbolicityGTest, testCustomizedAPSP){
	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input/snap/CA-AstroPh.txt");

//	  //closeness
//	  std::clock_t start1;
//	  double duration_closeness;
//	  start1 = std::clock();
//
//
//	  duration_closeness = (std::clock() - start1) / (double) CLOCKS_PER_SEC;
//	  INFO("Duration of Closeness:", duration_closeness);

	  //CustomizedAPSP for all nodes
	  std::clock_t start2;
	  double duration_cAPSP;
	  start2 = std::clock();
	  CustomizedAPSP cAPSP(G);
	  cAPSP.run();
	  auto distances = cAPSP.getExactDistances();
	  duration_cAPSP = (std::clock() - start2) / (double) CLOCKS_PER_SEC;
	  INFO("Duration of CustomizedAPSP:", duration_cAPSP);

//	  //CustomizedAPSP for 100 nodes
//	  std::clock_t start3;
//	  double duration_cAPSP_hundret;
//	  start3 = std::clock();
//	  CustomizedAPSP cAPSPHundret(G, topk.topkNodesList());
//	  cAPSPHundret.run();
//	  duration_cAPSP_hundret = (std::clock() - start3) / (double) CLOCKS_PER_SEC;
//	  INFO("Duration of CustomizedAPSPHundret:", duration_cAPSP_hundret);

	  //APSP
	  std::clock_t start4;
	  double duration_apsp;
	  start4 = std::clock();
	  APSP apsp(G);
	  apsp.run();
	  duration_apsp = (std::clock() - start4) / (double) CLOCKS_PER_SEC;
	  INFO("Duration of APSP:", duration_apsp);

	  for(node u=0; u < G.numberOfNodes(); ++u){
		  for(node v=u; v < G.numberOfNodes(); ++v){
			  EXPECT_TRUE(distances.element(u,v) == apsp.getDistance(u,v));
		  }
	  }

   }

  TEST_F(HyperbolicityGTest, testCustomizedAPSP2){

  }
}
