#define _USE_MATH_DEFINES
#include "DGGS_3D.h"

#include <cmath>

DGGS_3D::DGGS_3D(bool volume, int numInitialFaces, double surfaceRadialRatio, double gridRadius)
	: volume(volume), numInitialFaces(numInitialFaces), surfaceRadialRatio(surfaceRadialRatio), gridRadius(gridRadius) {

	double meanArea = 4.0 * M_PI / numInitialFaces;
	idealDegenChildren = surfaceRadialRatio / (2.0 * sqrt(meanArea));

	degenChildren = (int)(idealDegenChildren + 0.5);
	extraLayers = degenChildren - 1;

	layerBounds.push_back(gridRadius);
	layerSL.push_back(0);
	maxSL = 0;
}


// Generate layers for the next level of subdivision
void DGGS_3D::generateNextSL() {

	int start = 0;
	for (int i = 0; i < maxSL; i++) {
		start = 2 * start + extraLayers * i + 1;
	}
	int numElements = layerBounds.size() - start;

	double maxRad = layerBounds[start];
	double midRad = maxRad / 2.0;
	double minRad = 0.0;

	layerBounds.push_back(midRad);
	layerSL.push_back(layerSL[start]);

	// Create children of degenerate layer
	if (idealDegenChildren > 0.75) {

		double prevUpper = 0.0;
		for (int i = 0; i < degenChildren; i++) {

			// Splitting point changes if we want equal volume
			double upperRad;
			if (volume) {

				int n = degenChildren - i - 1;
				// Equal cell volume formula
				upperRad = (n == 0.0) ? maxRad : pow((n * pow(prevUpper, 3) + pow(maxRad, 3)) / (1 + n), 1.0 / 3.0);
				prevUpper = upperRad;
			}
			else {

				int j = degenChildren - 1 - i;
				double uPerc = (degenChildren - j) / degenChildren;
				double lPerc = j / degenChildren;

				upperRad = uPerc * maxRad + lPerc * midRad;
			}
			layerBounds.push_back(upperRad);
			layerSL.push_back(layerSL[start]); // no surface subivision
		}
	}
	// No extra layers and top layer gets subdivided
	else {
		layerBounds.push_back(maxRad);
		layerSL.push_back(layerSL[start] + 1); // surface subivision
	}

	// Create two children layers for all non-degenerate layers
	for (int i = 1; i <= numElements - 1; i ++) {

		minRad = maxRad;
		maxRad = layerBounds[start + i];
		
		if (volume) {
			midRad = pow((pow(maxRad, 3) + pow(minRad, 3)) / 2.0, 1.0 / 3.0);
		}
		else {
			midRad = 0.5 * maxRad + 0.5 * minRad;
		}
		layerBounds.push_back(midRad);
		layerBounds.push_back(maxRad);
		layerSL.push_back(layerSL[start + i] + 1); // surface subivision
		layerSL.push_back(layerSL[start + i] + 1); // surface subivision
	}
	maxSL++;
}


// Returns the layer that spans the provided radius at the given SL
//
// radius - target radius
// SL - target subdivision level
// return - layer ID
int DGGS_3D::layer(double radius, int SL) {

	// Generate necessary number of SL
	while (maxSL < SL) {
		generateNextSL();
	}

	// Find starting index in list
	int start = 0;
	for (int i = 0; i < SL; i++) {
		start = 2 * start + extraLayers * i + 1;
	}

	int first = start - 1;
	int last = layerBounds.size() - 1;
	int mid;

	// Binary search to find layer that contains radius
	// Search for smallest element >= radius
	while (last - first > 1) {

		mid = (first + last) / 2;
		if (radius > layerBounds[mid]) {
			first = mid;
			mid++;
		}
		else {
			last = mid;
		}
	}
	return mid;
}


// Calculates information about the provided layer
// Determines depth, boundary case, degenerate layer, child of degenerate layer, and end layer
//
// layer - ID of the target layer
// return - struct with all layer info
LayerInfo DGGS_3D::layerInfo(int layer) {

	LayerInfo info;

	// Determine depth and if layer is a degen layer or end layer
	info.depth = 0;
	info.end = false;
	info.degen = (layer == 0) ? true : false;

	// Loop to determine depth and check for degen layer and end layer
	int end = 0;
	while (layer > end) {

		if (layer == end + 1) {
			info.degen = true;
		}

		end = 2 * end + extraLayers * info.depth + extraLayers + 2;
		info.depth++;

		if (layer == end) {
			info.end = true;
		}
	}

	// Determine if layer is child of a degen layer and if it is a boundary layer
	info.degenChild = (info.depth == 1) ? true : false;

	// Loop to check boundary cases and check for child of SG layer
	int inc = (1 + extraLayers) * pow(2, info.depth - 1);
	int top = end;
	int i = 0;
	while (layer <= top + 1) {

		top -= inc;

		// Condition for being child of degen layer
		if (i == info.depth - 2 && layer <= top) {
			info.degenChild = true;
		}

		// Test for the two different boundary cases
		if (layer == top) {
			info.boundary = -1;
			break;
		}
		else if (layer == top + 1) {
			info.boundary = 1;
			break;
		}
		inc /= 2;
		i++;
	}
	return info;
}


// Returns the layer parent to the given layer
//
// layer - ID of the target layer
// info - info about the target layer
// return - ID of parent layer
int DGGS_3D::parentLayer(int layer, const LayerInfo & info) {

	if (layer == 0) {
		return -1;
	}

	if (info.degenChild) {
		int parent = 0;
		for (int i = 0; i < info.depth - 1; i++) {
			parent = 2 * parent + extraLayers * i + 1;
		}
		return parent;
	}
	else {
		return (layer - extraLayers - (info.depth - 1) * extraLayers - 1) / 2;
	}
}


// Returns the layers children to the given layer
//
// layer - ID of the target layer
// info - info about the target layer
// return - vector of all children layer IDs (2 children for normal layers, potentially more for degen layers)
std::vector<int> DGGS_3D::childrenLayers(int layer, const LayerInfo& info) {

	std::vector<int> children;

	if (info.degen) {
		for (int i = 1; i <= extraLayers + 2; i++) {
			children.push_back(layer * 2 + (info.depth) * extraLayers + i);
		}
	}
	else {
		children.push_back(layer * 2 + (info.depth) * extraLayers + extraLayers + 1);
		children.push_back(layer * 2 + (info.depth) * extraLayers + extraLayers + 2);
	}

	return children;
}



// Returns the layers neighbour (above and below) the given layer
//
// first - layer above (larger radius)
// second - layer below (smaller radius)
// return - pair containing the ID of neighbour layers. -1 if given neighbour does not exist
std::pair<int, int> DGGS_3D::neighbourLayers(int layer, const LayerInfo & info) {

	int above = (info.end) ? -1 : layer + 1;
	int below = (info.degen) ? -1 : layer - 1;

	return std::pair<int, int>(above, below);
}