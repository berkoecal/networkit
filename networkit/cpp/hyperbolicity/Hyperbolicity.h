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
#include "SymmetricMatrix.h"

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
	double getHyperbolicity() const {if (!hasRun) throw std::runtime_error("Call run method first"); return hyperbolicity_value;}

protected:
	const Graph& graph;
	double hyperbolicity_value=0;
	double naiveAlgorithm();
	double HYP(); 				//original name of corresponding paper BCCM15
	node centralNode(SymMatrix<edgeweight, node> const& distances);
	bool is_acceptable(SymMatrix<edgeweight, node> const& distances, node const& x, node const& y, node const& v, edgeweight const& current_lower_bound, edgeweight const& eccentricity);
	bool is_valuable(SymMatrix<edgeweight, node> const& distances,node const& x, node const& y, node const& v, edgeweight const& current_lower_bound, edgeweight const& eccentricity, node const& central_node);
};
 
} /* namespace NetworKit */

#endif // HYPERBOLICITY_H
