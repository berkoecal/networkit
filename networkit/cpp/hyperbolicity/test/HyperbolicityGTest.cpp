/*
 * HyperbolicityGTest.cpp
 *       Created on: 18.11.2017
 *       Author: Berk Öcal
 */
#define _BSD_SOURCE

#include<sys/time.h>
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
#include "../../distance/Diameter.h"
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

  TEST_F(HyperbolicityGTest, testHyperbolicityCloseness) {
	  auto get_wall_time = []()->double{
		  struct timeval time;
		  if(gettimeofday(&time, NULL)){
			  return 0;
		  }

		  return (double) time.tv_sec + (double) time.tv_usec * .000001;
	  };
	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/032233cit-HepPh");
	  double begin = get_wall_time();
	  auto k = (G.numberOfNodes())*0.1;
	  TopCloseness topK(G, k);
	  topK.run();
	  double end = get_wall_time();
	  INFO("Topk-Closeness with ", k, " nodes took ", end-begin);
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
	  Graph G = Reader.read("input_borassi/snap/Amazon0302.txt");

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

		  if(component.numberOfNodes() > 1000){
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
	  Graph G = Reader.read("input_borassi/045425soc-Slashdot0811");

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

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe1) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/004527p2p-Gnutella08");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe2) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/004583wiki-Vote");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe3) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/004667oregon2_010331");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe4) {

  	  auto Reader = SNAPGraphReader();
  	  Graph G = Reader.read("input_borassi/008362p2p-Gnutella04");
  	  Hyperbolicity hyperbolicity(G);
  	  hyperbolicity.run();
  	  double value = hyperbolicity.getHyperbolicity();

  	  INFO("Final Hyperbolicity over all connected components: ", value);
    }
  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe5) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/008773ca-CondMat");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe6) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/009577ca-AstroPh");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe7) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/009710ASEdges10_2011.edgelist");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe8) {

  	  auto Reader = SNAPGraphReader();
  	  Graph G = Reader.read("input_borassi/010969email-Enron");
  	  Hyperbolicity hyperbolicity(G);
  	  hyperbolicity.run();
  	  double value = hyperbolicity.getHyperbolicity();

  	  INFO("Final Hyperbolicity over all connected components: ", value);
    }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe9) {

  	  auto Reader = SNAPGraphReader();
  	  Graph G = Reader.read("input_borassi/013301p2p-Gnutella25");
  	  Hyperbolicity hyperbolicity(G);
  	  hyperbolicity.run();
  	  double value = hyperbolicity.getHyperbolicity();

  	  INFO("Final Hyperbolicity over all connected components: ", value);
    }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe10) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/015549as_20100120.caida");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe11) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/022153email-EuAll");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe12) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/020095p2p-Gnutella30");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe13) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/024982cit-HepTh");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe14) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/030935soc-Epinions1");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe15) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/032233cit-HepPh");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe16) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/045425soc-Slashdot0811");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe17) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/050219soc-sign-epinions");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe18) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/004950ca-HepPh");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityAbgabe19) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/003702ca-HepTh");
	  Hyperbolicity hyperbolicity(G);
	  hyperbolicity.run();
	  double value = hyperbolicity.getHyperbolicity();

	  INFO("Final Hyperbolicity over all connected components: ", value);
  }


/*
  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe1) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/004527p2p-Gnutella08");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;
	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe2) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/004583wiki-Vote");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe3) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/004667oregon2_010331");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe4) {

  	  auto Reader = SNAPGraphReader();
  	  Graph G = Reader.read("input_borassi/008362p2p-Gnutella04");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
    }
  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe5) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/008773ca-CondMat");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe6) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/009577ca-AstroPh");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe7) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/009710ASEdges10_2011.edgelist");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe8) {

  	  auto Reader = SNAPGraphReader();
  	  Graph G = Reader.read("input_borassi/010969email-Enron");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
    }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe9) {

  	  auto Reader = SNAPGraphReader();
  	  Graph G = Reader.read("input_borassi/013301p2p-Gnutella25");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
    }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe10) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/015549as_20100120.caida");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe11) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/022153email-EuAll");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe12) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/020095p2p-Gnutella30");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe13) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/024982cit-HepTh");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe14) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/030935soc-Epinions1");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe15) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/032233cit-HepPh");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe16) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/045425soc-Slashdot0811");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe17) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/050219soc-sign-epinions");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe18) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/004950ca-HepPh");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }

  TEST_F(HyperbolicityGTest, testHyperbolicityDiamAbgabe19) {

	  auto Reader = SNAPGraphReader();
	  Graph G = Reader.read("input_borassi/003702ca-HepTh");
	  Diameter diam(G,1);
	  diam.run();
	  auto val = diam.getDiameter();
	  auto value = val.first;

	  INFO("Diameter: ", value);
  }
*/


  TEST_F(HyperbolicityGTest, testHyperbolicityDistCheck1) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/004527p2p-Gnutella08");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, ceil(n*0.05), 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, ceil(n*0.05), 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, ceil(n*0.05), 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundAdvanced(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundAdvanced(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundAdvanced(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

  }

TEST_F(HyperbolicityGTest, testHyperbolicityDistCheck2) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/013301p2p-Gnutella25");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, ceil(n*0.05), 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, ceil(n*0.05), 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, ceil(n*0.05), 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundAdvanced(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundAdvanced(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundAdvanced(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

 }

TEST_F(HyperbolicityGTest, testHyperbolicityDistCheck3) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/010969email-Enron");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, ceil(n*0.05), 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, ceil(n*0.05), 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, ceil(n*0.05), 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundAdvanced(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundAdvanced(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundAdvanced(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

 }

TEST_F(HyperbolicityGTest, testHyperbolicityDistCheck4) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/024982cit-HepTh");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, ceil(n*0.05), 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, ceil(n*0.05), 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, ceil(n*0.05), 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundAdvanced(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundAdvanced(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundAdvanced(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

 }

TEST_F(HyperbolicityGTest, testHyperbolicityDistCheck5) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/032233cit-HepPh");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, ceil(n*0.05), 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, ceil(n*0.05), 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, ceil(n*0.05), 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundAdvanced(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundAdvanced(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundAdvanced(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundAdvanced(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

 }


TEST_F(HyperbolicityGTest, testHyperbolicityDistCheckU1) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/004527p2p-Gnutella08");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, 20, 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, 20, 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, 20, 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundByComp(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundByComp(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundByComp(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

 }

TEST_F(HyperbolicityGTest, testHyperbolicityDistCheckU2) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/013301p2p-Gnutella25");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, 20, 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, 20, 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, 20, 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundByComp(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundByComp(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundByComp(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

 }

TEST_F(HyperbolicityGTest, testHyperbolicityDistCheckU3) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/010969email-Enron");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, 20, 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, 20, 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, 20, 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundByComp(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundByComp(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundByComp(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

 }

TEST_F(HyperbolicityGTest, testHyperbolicityDistCheckU4) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/024982cit-HepTh");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, 20, 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, 20, 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, 20, 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundByComp(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundByComp(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundByComp(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

 }

TEST_F(HyperbolicityGTest, testHyperbolicityDistCheckU5) {

	auto Reader = SNAPGraphReader();
	Graph G = Reader.read("input_borassi/032233cit-HepPh");
	auto n = G.numberOfNodes();

	int errorsGD=0;
	double errorSumGD=0.0;
	double relativeErrorGD=0.0;

	int errorsT=0;
	double errorSumT=0.0;
	double relativeErrorT=0.0;

	int errorsD=0;
	double errorSumD=0.0;
	double relativeErrorD=0.0;

	CustomizedAPSP _apsp(G);
	_apsp.run();

	CustomizedAPSP cAPSP(G, 20, 20, CentralNodeMethod::GROUPDEGREE);
	cAPSP.run();

	CustomizedAPSP cAPSP2(G, 20, 20, CentralNodeMethod::TOPCLOSENESS);
	cAPSP2.run();

	CustomizedAPSP cAPSP3(G, 20, 20, CentralNodeMethod::DEGREE);
	cAPSP3.run();

	for(node u= 0; u < G.numberOfNodes(); ++u){
		for(node v=u; v < G.numberOfNodes(); ++v){
			if(_apsp.getDistance(u,v) != cAPSP.distUpperBoundByComp(u,v)){
				errorsGD++;
				errorSumGD +=(cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorGD+= (cAPSP.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP2.distUpperBoundByComp(u,v)){
				errorsT++;
				errorSumT +=(cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorT+= (cAPSP2.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}

			if(_apsp.getDistance(u,v) != cAPSP3.distUpperBoundByComp(u,v)){
				errorsD++;
				errorSumD +=(cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v));
				relativeErrorD+= (cAPSP3.distUpperBoundByComp(u,v) - _apsp.getDistance(u,v))/_apsp.getDistance(u,v);
			}
		}
	}

	INFO("Relative Error Average (Group Degree 10%): ", relativeErrorGD/(double)errorsGD);
	INFO("Relative Error Average (TOPCLOSENESS 10%): ", relativeErrorT/(double) errorsT);
	INFO("Relative Error Average (Degree 10%): ", relativeErrorD/ (double) errorsD);

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
