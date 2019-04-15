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

	int sl = 0;

	std::array<Tri, 2> twoToOne() const;
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

	double volume() const;
};

class TestGrid {

public:
	TestGrid(int n);

	void subdivide(bool volume);

//private:
	std::unordered_map<std::string, TriCell> map;
	int curDepth;
};

