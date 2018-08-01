#define _USE_MATH_DEFINES
#include "DGGS_3D.h"

#include <cmath>

DGGS_3D::DGGS_3D(bool volume, int numInitialFaces, double surfaceRadialRatio, double gridRadius)
	: volume(volume), numInitialFaces(numInitialFaces), surfaceRadialRatio(surfaceRadialRatio), gridRadius(gridRadius) {

	double faceArea = 4.0 * M_PI / numInitialFaces;
	idealNumDegenLayers = surfaceRadialRatio / (2.0 * sqrt(faceArea));

	numDegenLayers = (int)(idealNumDegenLayers + 0.5);
	numExtraLayers = numDegenLayers - 1;

	layerBounds.push_back(gridRadius);
	maxSL = 0;
}


// Generate the bounds of layers for one more level of subdivision
void DGGS_3D::generateLayerBounds() {

	int start = 0;
	for (int i = 0; i < maxSL; i++) {
		start = 2 * start + numExtraLayers * i + 1;
	}
	int numElements = layerBounds.size() - start;

	double maxRad = layerBounds[start];
	double midRad = maxRad / 2.0;
	double minRad = 0.0;

	layerBounds.push_back(midRad);

	// degenerate case here
	if (numDegenLayers > 0.75) {

		// Round number of layers to nearest integer
		double prevUpper = 0.0;
		double numLayersR = (int)(numDegenLayers + 0.5);

		for (int i = 0; i < numLayersR; i++) {

			// Splitting point changes if we want equal volume
			double upperRad;
			if (volume) {

				int n = numLayersR - i - 1;
				// Equal cell volume formula
				upperRad = (n == 0.0) ? maxRad : pow((n * pow(prevUpper, 3) + pow(maxRad, 3)) / (1 + n), 1.0 / 3.0);
				prevUpper = upperRad;
			}
			else {

				int j = numLayersR - 1 - i;
				double uPerc = (numLayersR - j) / numLayersR;
				double lPerc = j / numLayersR;

				upperRad = uPerc * maxRad + lPerc * midRad;
			}

			layerBounds.push_back(upperRad);
		}
	}
	// Upper region is too short and fat: surface subdivide upper region
	else {
		layerBounds.push_back(maxRad);
	}

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
	}

	maxSL++;
}


// Returns the layer that spans the provided radius at the given SL
//
// radius - target radius
// SL - target subdivision level
// return - layer ID
int DGGS_3D::layer(double radius, int SL) {

	while (maxSL < SL) {
		generateLayerBounds();
	}

	int base = 0;
	for (int i = 0; i < SL; i++) {
		base = 2 * base + numExtraLayers * i + 1;
	}

	int first = base - 1;
	int last = layerBounds.size() - 1;
	int mid;

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

	// Determine depth and if layer is an SG layer or end layer
	info.depth = 0;
	info.end = false;
	info.degen = (layer == 0) ? true : false;

	// Loop to determine depth and check for SG layer and end layer
	int end = 0;
	while (layer > end) {

		if (layer == end + 1) {
			info.degen = true;
		}

		end = 2 * end + numExtraLayers * info.depth + numExtraLayers + 2;
		info.depth++;

		if (layer == end) {
			info.end = true;
		}
	}

	// Determine if layer is child of an SG layer and if it is a boundary layer
	info.boundary = 0;
	info.degenChild = (info.depth == 1) ? true : false;

	// Loop to check boundary cases and check for child of SG layer
	int inc = (1 + numExtraLayers) * pow(2, info.depth - 1);
	int top = end;
	int i = 0;
	while (layer <= top + 1) {

		top -= inc;

		// Condition for being child of SG layer
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
			parent = 2 * parent + numExtraLayers * i + 1;
		}
		return parent;
	}
	else {
		return (layer - numExtraLayers - (info.depth - 1) * numExtraLayers - 1) / 2;
	}
}


// Returns the layers children to the given layer
//
// layer - ID of the target layer
// info - info about the target layer
// return - vector of all children layers (2 children for normal layers and 
std::vector<int> DGGS_3D::childrenLayers(int layer, const LayerInfo& info) {

	std::vector<int> children;

	if (info.degen) {
		for (int i = 1; i <= numExtraLayers + 2; i++) {
			children.push_back(layer * 2 + (info.depth) * numExtraLayers + i);
		}
	}
	else {
		children.push_back(layer * 2 + (info.depth) * numExtraLayers + numExtraLayers + 1);
		children.push_back(layer * 2 + (info.depth) * numExtraLayers + numExtraLayers + 2);
	}

	return children;
}