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

struct Tri {
	Tri() = default;
	Tri(int n);

	glm::dvec3 v0;
	glm::dvec3 v1;
	glm::dvec3 v2;

	void fourToOne(std::array<Tri, 4>& out);
	void nineToOne(std::array<Tri, 9>& out);
};

struct TriCell {
	Tri tri;
	double minRad;
	double maxRad;
	CT cellType;

	void fillRenderable(Renderable& r, const glm::vec3& colour, bool curved);
};

class TestGrid {

public:
	TestGrid();

	void subdivide(bool volume);
	void subdivide2(bool volume);
	void subdivide3(bool volume);

//private:
	std::unordered_map<std::string, TriCell> map;
	int curDepth;
};

