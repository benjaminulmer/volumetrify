#define _USE_MATH_DEFINES
#include "TestGrid.h"

#include "Geometry.h"

void Cell::fillRenderable(Renderable & r, const glm::vec3 & colour, bool curved) {

	r.lineColour = colour;
	glm::dvec3 o(0.0);

	Geometry::createArcR(v0.toCartesian(maxRad), v1.toCartesian(maxRad), o, r);
	Geometry::createArcR(v0.toCartesian(maxRad), v2.toCartesian(maxRad), o, r);
	Geometry::createArcR(v1.toCartesian(maxRad), v2.toCartesian(maxRad), o, r);

	Geometry::createArcR(v0.toCartesian(minRad), v1.toCartesian(minRad), o, r);
	Geometry::createArcR(v0.toCartesian(minRad), v2.toCartesian(minRad), o, r);
	Geometry::createArcR(v1.toCartesian(minRad), v2.toCartesian(minRad), o, r);

	if (!curved) {
		r.lineColour = glm::vec3(1.f, 0.f, 0.f);
		Geometry::createLineR(v0.toCartesian(maxRad), v1.toCartesian(maxRad), r);
		Geometry::createLineR(v0.toCartesian(maxRad), v2.toCartesian(maxRad), r);
		Geometry::createLineR(v1.toCartesian(maxRad), v2.toCartesian(maxRad), r);

		Geometry::createLineR(v0.toCartesian(minRad), v1.toCartesian(minRad), r);
		Geometry::createLineR(v0.toCartesian(minRad), v2.toCartesian(minRad), r);
		Geometry::createLineR(v1.toCartesian(minRad), v2.toCartesian(minRad), r);
		r.lineColour = colour;
	}

	Geometry::createLineR(v0.toCartesian(maxRad), v0.toCartesian(minRad), r);
	Geometry::createLineR(v1.toCartesian(maxRad), v1.toCartesian(minRad), r);
	Geometry::createLineR(v2.toCartesian(maxRad), v2.toCartesian(minRad), r);
}

void Cell::fourToOne(std::array<Cell, 4>& out) {

	SphCoord mid01 = SphCoord((v0.toCartesian(1.0) + v1.toCartesian(1.0)) / 2.0);
	SphCoord mid02 = SphCoord((v0.toCartesian(1.0) + v2.toCartesian(1.0)) / 2.0);
	SphCoord mid12 = SphCoord((v1.toCartesian(1.0) + v2.toCartesian(1.0)) / 2.0);

	out[0] = Cell();
	out[0].v0 = v0;
	out[0].v1 = mid01;
	out[0].v2 = mid02;

	out[1] = Cell();
	out[1].v0 = v1;
	out[1].v1 = mid01;
	out[1].v2 = mid12;

	out[2] = Cell();
	out[2].v0 = v2;
	out[2].v1 = mid02;
	out[2].v2 = mid12;

	out[3] = Cell();
	out[3].v0 = mid01;
	out[3].v1 = mid02;
	out[3].v2 = mid12;
}

TestGrid::TestGrid() {

	char num = 'a';

	// Top
	for (int i = 0; i < 5; i++) {

		std::string code1 = std::string(1, num++);
		Cell c1;
		c1.v0 = SphCoord(M_PI_2, 0.0);
		c1.v1 = SphCoord(atan(0.5), i * (2.0 * M_PI / 5.0));
		c1.v2 = SphCoord(atan(0.5), (i + 1) * (2.0 * M_PI / 5.0));
		c1.maxRad = 6371000.0 * (4.0 / 3.0);
		c1.minRad = 0.0;
		c1.cellType = CT::SG;

		std::string code2 = std::string(1, num++);
		Cell c2;
		c2.v0 = SphCoord(-atan(0.5), i * (2.0 * M_PI / 5.0) + (M_PI / 5.0));
		c2.v1 = SphCoord(atan(0.5), i * (2.0 * M_PI / 5.0));
		c2.v2 = SphCoord(atan(0.5), (i + 1) * (2.0 * M_PI / 5.0));
		c2.maxRad = 6371000.0 * (4.0 / 3.0);
		c2.minRad = 0.0;
		c2.cellType = CT::SG;

		map[code1] = c1;
		map[code2] = c2;
	}

	// Bottom
	for (int i = 0; i < 5; i++) {

		std::string code1 = std::string(1, num++);
		Cell c1;
		c1.v0 = SphCoord(-M_PI_2, 0.0);
		c1.v1 = SphCoord(-atan(0.5), i * (2.0 * M_PI / 5.0) + (M_PI / 5.0));
		c1.v2 = SphCoord(-atan(0.5), (i + 1) * (2.0 * M_PI / 5.0) + (M_PI / 5.0));
		c1.maxRad = 6371000.0 * (4.0 / 3.0);
		c1.minRad = 0.0;
		c1.cellType = CT::SG;

		std::string code2 = std::string(1, num++);
		Cell c2;
		c2.v0 = SphCoord(atan(0.5), (i + 1) * (2.0 * M_PI / 5.0));
		c2.v1 = SphCoord(-atan(0.5), i * (2.0 * M_PI / 5.0) + (M_PI / 5.0));
		c2.v2 = SphCoord(-atan(0.5), (i + 1) * (2.0 * M_PI / 5.0) + (M_PI / 5.0));
		c2.maxRad = 6371000.0 * (4.0 / 3.0);
		c2.minRad = 0.0;
		c2.cellType = CT::SG;

		map[code1] = c1;
		map[code2] = c2;
	}
	curDepth = 1;
}


void TestGrid::subdivide(bool volume) {

	for (auto p : map) {
		if (p.first.length() != curDepth) continue;

		Cell& c = p.second;

		double midRad;
		if (c.cellType == CT::SG) {
			midRad = 0.5 * c.maxRad + 0.5 * c.minRad;
		}
		else {

			if (volume) {
				double n = (c.cellType == CT::SG) ? 4.0 : 1.0;
				double base = (n / (n + 1)) * ((1 / n) * pow(c.maxRad, 3) + pow(c.minRad, 3));
				midRad = pow(base, 1.0 / 3.0);
			}
			else {
				midRad = 0.5 * c.maxRad + 0.5 * c.minRad;
			}
		}

		std::array<Cell, 4> children;
		c.fourToOne(children);

		Cell top1 = children[0];
		top1.maxRad = c.maxRad;
		top1.minRad = midRad;
		top1.cellType = CT::NG;

		Cell top2 = children[1];
		top2.maxRad = c.maxRad;
		top2.minRad = midRad;
		top2.cellType = CT::NG;

		Cell top3 = children[2];
		top3.maxRad = c.maxRad;
		top3.minRad = midRad;
		top3.cellType = CT::NG;

		Cell top4 = children[3];
		top4.maxRad = c.maxRad;
		top4.minRad = midRad;
		top4.cellType = CT::NG;

		map[p.first + "0"] = top1;
		map[p.first + "1"] = top2;
		map[p.first + "2"] = top3;
		map[p.first + "3"] = top4;

		if (c.cellType == CT::NG) {

			Cell bot1 = children[0];
			bot1.maxRad = midRad;
			bot1.minRad = c.minRad;
			bot1.cellType = CT::NG;

			Cell bot2 = children[1];
			bot2.maxRad = midRad;
			bot2.minRad = c.minRad;
			bot2.cellType = CT::NG;

			Cell bot3 = children[2];
			bot3.maxRad = midRad;
			bot3.minRad = c.minRad;
			bot3.cellType = CT::NG;

			Cell bot4 = children[3];
			bot4.maxRad = midRad;
			bot4.minRad = c.minRad;
			bot4.cellType = CT::NG;

			map[p.first + "4"] = bot1;
			map[p.first + "5"] = bot2;
			map[p.first + "6"] = bot3;
			map[p.first + "7"] = bot4;
		}
		else {
			Cell bottom = c;
			bottom.maxRad = midRad;
			bottom.cellType = CT::SG;

			map[p.first + "4"] = bottom;
		}
	}
	curDepth++;
}
