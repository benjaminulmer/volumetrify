#pragma once

#include "Renderable.h"
#include "SphCoord.h"

#include <glm/glm.hpp>

#include <string>
#include <unordered_map>

enum class CT {
	SG,
	F_NG,
	NG
};


class Cell {

public:
	SphCoord v0;
	SphCoord v1;
	SphCoord v2;
	double maxRad;
	double minRad;

	void fillRenderable(Renderable& r, const glm::vec3& colour, bool curved);

	CT cellType;
};

class TestGrid {

public:
	TestGrid();

	void subdivide(bool volume);

//private:
	std::unordered_map<std::string, Cell> map;
	int curDepth;
};

