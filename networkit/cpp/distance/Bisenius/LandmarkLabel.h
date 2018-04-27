/*
 * LandmarkLabel.h
 */

#ifndef LANDMARKLABEL_H_
#define LANDMARKLABEL_H_

#include "../../graph/Graph.h"

namespace NetworKit {

/**
* This class encapsulates the landmark labels user in the PrunedAPSP algorithm.
*/
class LandmarkLabel {
	
public:
	LandmarkLabel();
	std::vector<node> nodes;
	std::vector<double> distances;
	std::vector<node> parents;
};

} /* namespace NetworKit */
#endif /* LANDMARKLABEL_H_ */
