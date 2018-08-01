#pragma once

#include <vector>

struct LayerInfo {
	int depth;
	int boundary;
	bool degen;
	bool degenChild;
	bool end;
};

class DGGS_3D {

public:
	DGGS_3D(bool volume, int numInitialFaces, double surfaceRadialRatio, double gridRadius);

	int layer(double radius, int SL);

	LayerInfo layerInfo(int layer);
	int parentLayer(int layer, const LayerInfo& info);
	std::vector<int> childrenLayers(int layer, const LayerInfo& info);

//private:
	bool volume;

	int numInitialFaces;
	double surfaceRadialRatio;
	double gridRadius;

	double idealNumDegenLayers;
	int numDegenLayers;
	int numExtraLayers;

	int maxSL;
	std::vector<double> layerBounds;

	void generateLayerBounds();
};

