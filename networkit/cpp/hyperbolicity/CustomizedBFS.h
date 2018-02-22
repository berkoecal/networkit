/*
 * CustomizedBFS.h
 *
 *  Created on: 18.02.2018
 *      Author: Berk Öcal
 */
#ifndef CUSTOMIZEDBFS_H
#define CUSTOMIZEDBFS_H

#include "SymmetricMatrix.h"
#include "../graph/Graph.h"
#include <queue>

namespace NetworKit{

class CustomizedBFS{
public:
	CustomizedBFS(const Graph& G);
	virtual void run();
	SymMatrix<bool, node> getFarApartPairs();
	SymMatrix<edgeweight, node> getDistances();
	edgeweight getDistance(node u, node v);
	edgeweight getEccentricity(node v);

private:
	const Graph& graph;
	SymMatrix<bool, node> farApart;
	SymMatrix<edgeweight, node> distances;
	std::vector<edgeweight> eccentricity;

};

}
#endif /*CUSTOMIZEDBFS_H */
