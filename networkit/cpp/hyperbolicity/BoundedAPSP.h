/*
 * BoundedAPSP.h
 *
 *  Created on: 04.05.2018
 *      Author: berkoec
 */

#ifndef NETWORKIT_CPP_HYPERBOLICITY_BOUNDEDAPSP_H_
#define NETWORKIT_CPP_HYPERBOLICITY_BOUNDEDAPSP_H_

#include "SymmetricMatrix.h"
#include "../graph/Graph.h"
#include "../auxiliary/PrioQueue.h"
#include <queue>

namespace NetworKit {

//enum class CentralNodeMethod{
//		NATURAL,
//		DEGREE,
//		TOPCLOSENESS,
//		RANDOM
//};

class BoundedAPSP {
public:
	BoundedAPSP();
	virtual ~BoundedAPSP();

//private:
//	std::vector<node> nodesForBFS;
//	SymMatrix<edgeweight, node> exactDistances;
//	const count numberOfLandmarks;
//	CentralNodeMethod method;
//	edgeweight queryBound(node u, node v);
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_HYPERBOLICITY_BOUNDEDAPSP_H_ */
