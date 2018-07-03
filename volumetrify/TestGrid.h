#pragma once

#include "Renderable.h"
#include "SphCoord.h"

#include <glm/glm.hpp>

#include <array>
#include <string>
#include <unordered_map>

enum class CT {
	SG,
	NG
};

class SphPolygon {

};

class Cell {

public:
	SphCoord v0;
	SphCoord v1;
	SphCoord v2;
	double maxRad;
	double minRad;
	CT cellType;

	void fillRenderable(Renderable& r, const glm::vec3& colour, bool curved);

	void fourToOne(std::array<Cell, 4>& out);
	void nineToOne(std::array<Cell, 9>& out);
};

class TestGrid {

public:
	TestGrid();

	void subdivide(bool volume);

//private:
	std::unordered_map<std::string, Cell> map;
	int curDepth;
};

