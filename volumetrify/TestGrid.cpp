#define _USE_MATH_DEFINES
#include "TestGrid.h"

#include "Geometry.h"

#define RAD 6371000.0 * (4.0 / 3.0)
#include <iostream>
Tri::Tri(int n) {

	double target = 4 * M_PI / n;
	double actual = 2 * M_PI;

	v0 = SphCoord(0.0, 0.0, false).toCartesian(1.0);
	v1 = SphCoord(0.0, 120.0, false).toCartesian(1.0);
	v2 = SphCoord(0.0, -120.0, false).toCartesian(1.0);
	glm::dvec3 p = SphCoord(90.0, 0.0, false).toCartesian(1.0);

	if (n <= 2) return;

	double upper = 1.0;
	double lower = 0.0;

	glm::dvec3 tV0, tV1, tV2;
	while (abs(target - actual) > 0.0001) {

		double t = 0.5 * upper + 0.5 * lower;
		tV0 = Geometry::geomSlerp(v0, p, t);
		tV1 = Geometry::geomSlerp(v1, p, t);
		tV2 = Geometry::geomSlerp(v2, p, t);

		std::vector<SphCoord> poly;
		poly.push_back(SphCoord(tV0));
		poly.push_back(SphCoord(tV1));
		poly.push_back(SphCoord(tV2));
		actual = SphCoord::areaPolygon(poly, 1.0);

		if (actual > target) {
			lower = t;
		}
		else {
			upper = t;
		}
	}
	v0 = tV0 * RAD;
	v1 = tV1 * RAD;
	v2 = tV2 * RAD;
}


void TriCell::fillRenderable(Renderable& r, const glm::vec3 & colour, bool curved) {

	r.lineColour = colour;
	glm::dvec3 o(0.0);

	if (curved) {
		Geometry::createArcR(glm::normalize(tri.v0) * RAD * maxRad, glm::normalize(tri.v1) * RAD * maxRad, o, r);
		Geometry::createArcR(glm::normalize(tri.v0) * RAD * maxRad, glm::normalize(tri.v2) * RAD * maxRad, o, r);
		Geometry::createArcR(glm::normalize(tri.v1) * RAD * maxRad, glm::normalize(tri.v2) * RAD * maxRad, o, r);

		Geometry::createArcR(glm::normalize(tri.v0) * RAD * minRad, glm::normalize(tri.v1) * RAD * minRad, o, r);
		Geometry::createArcR(glm::normalize(tri.v0) * RAD * minRad, glm::normalize(tri.v2) * RAD * minRad, o, r);
		Geometry::createArcR(glm::normalize(tri.v1) * RAD * minRad, glm::normalize(tri.v2) * RAD * minRad, o, r);

		Geometry::createLineR(glm::normalize(tri.v0) * RAD * maxRad, glm::normalize(tri.v0) * RAD * minRad, r);
		Geometry::createLineR(glm::normalize(tri.v1) * RAD * maxRad, glm::normalize(tri.v1) * RAD * minRad, r);
		Geometry::createLineR(glm::normalize(tri.v2) * RAD * maxRad, glm::normalize(tri.v2) * RAD * minRad, r);
	}
	else {
		Geometry::createLineR(tri.v0 * maxRad, tri.v1 * maxRad, r);
		Geometry::createLineR(tri.v0 * maxRad, tri.v2 * maxRad, r);
		Geometry::createLineR(tri.v1 * maxRad, tri.v2 * maxRad, r);

		Geometry::createLineR(tri.v0 * minRad, tri.v1 * minRad, r);
		Geometry::createLineR(tri.v0 * minRad, tri.v2 * minRad, r);
		Geometry::createLineR(tri.v1 * minRad, tri.v2 * minRad, r);

		Geometry::createLineR(tri.v0 * maxRad, tri.v0 * minRad, r);
		Geometry::createLineR(tri.v1 * maxRad, tri.v1 * minRad, r);
		Geometry::createLineR(tri.v2 * maxRad, tri.v2 * minRad, r);
	}
}


void Tri::fourToOne(std::array<Tri, 4>& out) {

	glm::dvec3 mid01 = 0.5 * v0 + 0.5 * v1;
	glm::dvec3 mid02 = 0.5 * v0 + 0.5 * v2;
	glm::dvec3 mid12 = 0.5 * v1 + 0.5 * v2;

	out[0] = Tri();
	out[0].v0 = v0;
	out[0].v1 = mid01;
	out[0].v2 = mid02;

	out[1] = Tri();
	out[1].v0 = v1;
	out[1].v1 = mid01;
	out[1].v2 = mid12;

	out[2] = Tri();
	out[2].v0 = v2;
	out[2].v1 = mid02;
	out[2].v2 = mid12;

	out[3] = Tri();
	out[3].v0 = mid01;
	out[3].v1 = mid02;
	out[3].v2 = mid12;
}


TestGrid::TestGrid() {

	std::string code = "a";
	Tri t(2000);
	TriCell cell;
	cell.tri = t;
	cell.maxRad = 1.0;
	cell.minRad = 0.0;
	cell.cellType = CT::SG;

	map[code] = cell;
	curDepth = 1;

	//char num = 'a';

	//// Top
	//for (int i = 0; i < 5; i++) {

	//	std::string code1 = std::string(1, num++);
	//	TriCell c1;
	//	c1.tri.v0 = SphCoord(M_PI_2, 0.0).toCartesian(RAD);
	//	c1.tri.v1 = SphCoord(atan(0.5), i * (2.0 * M_PI / 5.0)).toCartesian(RAD);
	//	c1.tri.v2 = SphCoord(atan(0.5), (i + 1) * (2.0 * M_PI / 5.0)).toCartesian(RAD);

	//	c1.maxRad = 1.0;
	//	c1.minRad = 0.0;
	//	c1.cellType = CT::SG;

	//	std::string code2 = std::string(1, num++);
	//	TriCell c2;
	//	c2.tri.v0 = SphCoord(-atan(0.5), i * (2.0 * M_PI / 5.0) + (M_PI / 5.0)).toCartesian(RAD);
	//	c2.tri.v1 = SphCoord(atan(0.5), i * (2.0 * M_PI / 5.0)).toCartesian(RAD);
	//	c2.tri.v2 = SphCoord(atan(0.5), (i + 1) * (2.0 * M_PI / 5.0)).toCartesian(RAD);

	//	c2.maxRad = 1.0;
	//	c2.minRad = 0.0;
	//	c2.cellType = CT::SG;

	//	map[code1] = c1;
	//	map[code2] = c2;
	//}

	//// Bottom
	//for (int i = 0; i < 5; i++) {

	//	std::string code1 = std::string(1, num++);
	//	TriCell c1;
	//	c1.tri.v0 = SphCoord(-M_PI_2, 0.0).toCartesian(RAD);
	//	c1.tri.v1 = SphCoord(-atan(0.5), i * (2.0 * M_PI / 5.0) + (M_PI / 5.0)).toCartesian(RAD);
	//	c1.tri.v2 = SphCoord(-atan(0.5), (i + 1) * (2.0 * M_PI / 5.0) + (M_PI / 5.0)).toCartesian(RAD);

	//	c1.maxRad = 1.0;
	//	c1.minRad = 0.0;
	//	c1.cellType = CT::SG;

	//	std::string code2 = std::string(1, num++);
	//	TriCell c2;
	//	c2.tri.v0 = SphCoord(atan(0.5), (i + 1) * (2.0 * M_PI / 5.0)).toCartesian(RAD);
	//	c2.tri.v1 = SphCoord(-atan(0.5), i * (2.0 * M_PI / 5.0) + (M_PI / 5.0)).toCartesian(RAD);
	//	c2.tri.v2 = SphCoord(-atan(0.5), (i + 1) * (2.0 * M_PI / 5.0) + (M_PI / 5.0)).toCartesian(RAD);

	//	c2.maxRad = 1.0;
	//	c2.minRad = 0.0;
	//	c2.cellType = CT::SG;

	//	map[code1] = c1;
	//	map[code2] = c2;
	//}
	//curDepth = 1;
}


void TestGrid::subdivide(bool volume) {

	for (auto p : map) {
		if (p.first.length() != curDepth) continue;

		TriCell& c = p.second;

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

		std::array<Tri, 4> children;
		c.tri.fourToOne(children);

		TriCell top1;
		top1.tri = children[0];
		top1.maxRad = c.maxRad;
		top1.minRad = midRad;
		top1.cellType = CT::NG;

		TriCell top2;
		top2.tri = children[1];
		top2.maxRad = c.maxRad;
		top2.minRad = midRad;
		top2.cellType = CT::NG;

		TriCell top3;
		top3.tri = children[2];
		top3.maxRad = c.maxRad;
		top3.minRad = midRad;
		top3.cellType = CT::NG;

		TriCell top4;
		top4.tri = children[3];
		top4.maxRad = c.maxRad;
		top4.minRad = midRad;
		top4.cellType = CT::NG;

		map[p.first + "0"] = top1;
		map[p.first + "1"] = top2;
		map[p.first + "2"] = top3;
		map[p.first + "3"] = top4;

		if (c.cellType == CT::NG) {

			TriCell bot1;
			bot1.tri = children[0];
			bot1.maxRad = midRad;
			bot1.minRad = c.minRad;
			bot1.cellType = CT::NG;

			TriCell bot2;
			bot2.tri = children[1];
			bot2.maxRad = midRad;
			bot2.minRad = c.minRad;
			bot2.cellType = CT::NG;

			TriCell bot3;
			bot3.tri = children[2];
			bot3.maxRad = midRad;
			bot3.minRad = c.minRad;
			bot3.cellType = CT::NG;

			TriCell bot4;
			bot4.tri = children[3];
			bot4.maxRad = midRad;
			bot4.minRad = c.minRad;
			bot4.cellType = CT::NG;

			map[p.first + "4"] = bot1;
			map[p.first + "5"] = bot2;
			map[p.first + "6"] = bot3;
			map[p.first + "7"] = bot4;
		}
		else {
			TriCell bottom = c;
			bottom.maxRad = midRad;
			bottom.cellType = CT::SG;

			map[p.first + "4"] = bottom;
		}
	}
	curDepth++;
}


void TestGrid::subdivide2(bool volume) {

	for (auto p : map) {

		char num = 'a';

		if (p.first.length() != curDepth) continue;

		TriCell& c = p.second;

		std::array<Tri, 4> children;
		c.tri.fourToOne(children);

		double midRad = 0.5 * c.maxRad + 0.5 * c.minRad;

		double numLayers = (c.cellType == CT::SG) ? 13.0 : 2.0;
		double lowerRad = (c.cellType == CT::SG) ? midRad : c.minRad;

		for (double i = 0.0; i < numLayers; i++) {

			double u1 = (numLayers - i) / numLayers;
			double l1 = i / numLayers;

			double u2 = (numLayers - i - 1) / numLayers;
			double l2 = (i + 1) / numLayers;

			double maxRad = u1 * c.maxRad + l1 * lowerRad;
			double minRad = u2 * c.maxRad + l2 * lowerRad;

			TriCell top1;
			top1.tri = children[0];
			top1.maxRad = maxRad;
			top1.minRad = minRad;
			top1.cellType = CT::NG;

			TriCell top2;
			top2.tri = children[1];
			top2.maxRad = maxRad;
			top2.minRad = minRad;
			top2.cellType = CT::NG;

			TriCell top3;
			top3.tri = children[2];
			top3.maxRad = maxRad;
			top3.minRad = minRad;
			top3.cellType = CT::NG;

			TriCell top4;
			top4.tri = children[3];
			top4.maxRad = maxRad;
			top4.minRad = minRad;
			top4.cellType = CT::NG;

			map[p.first + std::string(1, num++)] = top1;
			map[p.first + std::string(1, num++)] = top2;
			map[p.first + std::string(1, num++)] = top3;
			map[p.first + std::string(1, num++)] = top4;
		}

		if (c.cellType == CT::SG) {
			TriCell bottom = c;
			bottom.maxRad = midRad;
			bottom.cellType = CT::SG;

			map[p.first + "4"] = bottom;
		}
	}
	curDepth++;
}


void TestGrid::subdivide3(bool volume) {

	for (auto p : map) {

		char num = 'a';

		if (p.first.length() != curDepth) continue;

		TriCell& c = p.second;

		std::array<Tri, 4> children;
		c.tri.fourToOne(children);

		double midRad = 0.5 * c.maxRad + 0.5 * c.minRad;

		double numLayers = (c.cellType == CT::SG) ? 6.0 : 2.0;
		double lowerRad = (c.cellType == CT::SG) ? midRad : c.minRad;

		for (double i = 0.0; i < numLayers; i++) {

			double u1 = (numLayers - i) / numLayers;
			double l1 = i / numLayers;

			double u2 = (numLayers - i - 1) / numLayers;
			double l2 = (i + 1) / numLayers;

			double maxRad = u1 * c.maxRad + l1 * lowerRad;
			double minRad = u2 * c.maxRad + l2 * lowerRad;

			if (c.cellType == CT::NG) {

				TriCell top1;
				top1.tri = children[0];
				top1.maxRad = maxRad;
				top1.minRad = minRad;
				top1.cellType = CT::NG;

				TriCell top2;
				top2.tri = children[1];
				top2.maxRad = maxRad;
				top2.minRad = minRad;
				top2.cellType = CT::NG;

				TriCell top3;
				top3.tri = children[2];
				top3.maxRad = maxRad;
				top3.minRad = minRad;
				top3.cellType = CT::NG;

				TriCell top4;
				top4.tri = children[3];
				top4.maxRad = maxRad;
				top4.minRad = minRad;
				top4.cellType = CT::NG;

				map[p.first + std::string(1, num++)] = top1;
				map[p.first + std::string(1, num++)] = top2;
				map[p.first + std::string(1, num++)] = top3;
				map[p.first + std::string(1, num++)] = top4;
			}
			else {
				TriCell top1 = c;
				top1.maxRad = maxRad;
				top1.minRad = minRad;
				top1.cellType = CT::NG;

				map[p.first + std::string(1, num++)] = top1;
			}
		}

		if (c.cellType == CT::SG) {
			TriCell bottom = c;
			bottom.maxRad = midRad;
			bottom.cellType = CT::SG;

			map[p.first + "4"] = bottom;
		}
	}
	curDepth++;
}