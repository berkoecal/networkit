/*
 * Hyperbolicity.h
 *
 *  Created on: 05.11.2017
 *      Author: Berk Öcal
 */

#ifndef HYPERBOLICITY_H
#define HYPERBOLICITY_H

#include "../graph/Graph.h"
#include "../base/Algorithm.h"

namespace NetworKit{
  
/**
 * @ingroup hyperbolicity
 * Class for computing the hyperbolicity of a graph
 */
  

class Hyperbolicity : public Algorithm
{
public:
	/**
	 * Creates the Hyperbolicity class for @a G
	 * 
	 * @param G The graph
	 */
	
	Hyperbolicity(const Graph& G);
	
	virtual ~Hyperbolicity() = default;
	/**
	 * Computes the hyperbolicity for the graph
	 */ 
	void run() override;
	
	int getHyperbolicity() const {if (!hasRun) throw std::runtime_error("Call run method first"); return hyperbolicity_value;}


protected:

	const Graph& G;
	int hyperbolicity_value;
	std::vector<std::vector<edgeweight> > distances;
};

} /* namespace NetworKit */

#endif // HYPERBOLICITY_H
