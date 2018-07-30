#pragma once

#include "Renderable.h"
#include "SphCoord.h"

#include <glm/glm.hpp>

#include <array>
#include <string>
#include <unordered_map>

struct LayerInfo {
	int depth;
	int boundary;
	bool SG;
	bool SG_Child;
	bool end;
};

enum class CT {
	SG,
	NG
};

struct Tri {
	Tri() = default;
	Tri(int n);

	glm::dvec3 v0;
	glm::dvec3 v1;
	glm::dvec3 v2;

	std::array<Tri, 4> fourToOne() const;
	std::array<Tri, 9> nineToOne() const;
};

struct TriCell {
	Tri tri;
	double minRad;
	double maxRad;
	CT cellType;
	std::string code;

	std::vector<TriCell> subdivide(bool volume) const;
	void fillRenderable(Renderable& r, const glm::vec3& colour, bool curved) const;
};

class TestGrid {

public:
	TestGrid();

	void subdivide(bool volume);

	static LayerInfo layerInfo(int layer, int n);
	static int parentLayer(int layer, int n, const LayerInfo& info);
	static std::vector<int> childrenLayers(int layer, int n, const LayerInfo& info);

//private:
	std::unordered_map<std::string, TriCell> map;
	int curDepth;
};

