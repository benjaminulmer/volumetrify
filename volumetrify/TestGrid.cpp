#define _USE_MATH_DEFINES
#include "TestGrid.h"

#include "Geometry.h"

#define RAD 6371000.0 * (4.0 / 3.0)


// Creates a triangle with spherical area 1/n * area of the whole sphere
Tri::Tri(int n) {

	double target = 4 * M_PI / n;
	double actual = 2 * M_PI;

	// Equilateral triangle around the equator
	v0 = SphCoord(0.0, 0.0, false).toCartesian(1.0);
	v1 = SphCoord(0.0, 120.0, false).toCartesian(1.0);
	v2 = SphCoord(0.0, -120.0, false).toCartesian(1.0);
	glm::dvec3 p = SphCoord(90.0, 0.0, false).toCartesian(1.0);

	if (n <= 2) return;

	double upper = 1.0;
	double lower = 0.0;

	// Binary search vertex positions until desired area is achieved
	glm::dvec3 tV0, tV1, tV2;
	while (abs(target - actual) > 0.0001) {

		double t = 0.5 * upper + 0.5 * lower;
		tV0 = Geometry::geomSlerp(v0, p, t);
		tV1 = Geometry::geomSlerp(v1, p, t);
		tV2 = Geometry::geomSlerp(v2, p, t);

		// Calculate area of new triangle
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
	// Scale final result to sphere of desired size
	v0 = tV0 * RAD;
	v1 = tV1 * RAD;
	v2 = tV2 * RAD;
}


// Four to one surface subdivision for a triangle
std::array<Tri, 4> Tri::fourToOne() const {

	std::array<Tri, 4> toReturn;

	// Edge midpoints
	glm::dvec3 mid01 = 0.5 * v0 + 0.5 * v1;
	glm::dvec3 mid02 = 0.5 * v0 + 0.5 * v2;
	glm::dvec3 mid12 = 0.5 * v1 + 0.5 * v2;

	toReturn[0] = Tri();
	toReturn[0].v0 = v0;
	toReturn[0].v1 = mid01;
	toReturn[0].v2 = mid02;
	toReturn[0].sl = sl + 1;

	toReturn[1] = Tri();
	toReturn[1].v0 = v1;
	toReturn[1].v1 = mid01;
	toReturn[1].v2 = mid12;
	toReturn[1].sl = sl + 1;

	toReturn[2] = Tri();
	toReturn[2].v0 = v2;
	toReturn[2].v1 = mid02;
	toReturn[2].v2 = mid12;
	toReturn[2].sl = sl + 1;

	toReturn[3] = Tri();
	toReturn[3].v0 = mid01;
	toReturn[3].v1 = mid02;
	toReturn[3].v2 = mid12;
	toReturn[3].sl = sl + 1;

	return toReturn;
}


// Subdivide TriCell into its children. Properly adapts to size of triangle
// Currently assumes 1-4 surface subdivision scheme.
std::vector<TriCell> TriCell::subdivide(bool volume) const {

	double ratio = 1.0;

	std::vector<TriCell> toReturn;
	std::array<Tri, 4> children = tri.fourToOne();
	char num = 'a';

	double midRad = 0.5 * maxRad + 0.5 * minRad;

	if (cellType == CT::NG) {

		// Split all surface children once along radial direction
		for (const Tri& t : children) {

			// Splitting point changes if we want equal volume
			double radSplit;
			if (volume) {
				radSplit = pow((pow(maxRad, 3) + pow(minRad, 3)) / 2.0, 1.0 / 3.0);
			}
			else {
				radSplit = midRad;
			}

			TriCell top;
			top.tri = t;
			top.maxRad = maxRad;
			top.minRad = radSplit;
			top.cellType = CT::NG;
			top.code = code + std::string(1, num++);
			toReturn.push_back(top);

			TriCell bottom;
			bottom.tri = t;
			bottom.maxRad = radSplit;
			bottom.minRad = minRad;
			bottom.cellType = CT::NG;
			bottom.code = code + std::string(1, num++);
			toReturn.push_back(bottom);
		}
	}
	// In SG case treat cell as upper and lower region made by splitting along radial midpoint
	else {

		// Caculate spherical area of triangle
		std::vector<SphCoord> points;
		points.push_back(tri.v0);
		points.push_back(tri.v1);
		points.push_back(tri.v2);
		double sa = SphCoord::areaPolygon(points, 1.0);

		// Find number of layers to split upper region into to satisfy:
		// sqrt(sa) == depth of cell
		double numLayers = ratio / (2.0 * sqrt(sa));

		// Upper region is too long and skinny: split upper region into one or more layers
		if (numLayers > 0.75) {

			// Round number of layers to nearest integer
			double prevLower = maxRad;
			double numLayersR = (int)(numLayers + 0.5);
			for (int i = 0; i < numLayersR; i++) {

				// Splitting point changes if we want equal volume
				double upperRad, lowerRad;
				if (volume) {

					int n = numLayersR - i - 1;
					upperRad = prevLower;
					// Equal cell volume formula
					lowerRad = (n == 0.0) ? midRad : pow((n * pow(prevLower, 3) + pow(midRad, 3)) / (1 + n), 1.0 / 3.0);
					prevLower = lowerRad;
				}
				else {

					double u1 = (numLayersR - i) / numLayersR;
					double l1 = i / numLayersR;

					double u2 = (numLayersR - i - 1) / numLayersR;
					double l2 = (i + 1) / numLayersR;

					upperRad = u1 * maxRad + l1 * midRad;
					lowerRad = u2 * maxRad + l2 * midRad;
				}

				TriCell cell = *this;
				cell.maxRad = upperRad;
				cell.minRad = lowerRad;
				cell.cellType = CT::NG;
				cell.code = code + std::string(1, num++);
				toReturn.push_back(cell);
			}
		}
		// Upper region is too short and fat: surface subdivide upper region
		else {
			for (const Tri& t : children) {

				TriCell top;
				top.tri = t;
				top.maxRad = maxRad;
				top.minRad = midRad;
				top.cellType = CT::NG;
				top.code = code + std::string(1, num++);
				toReturn.push_back(top);
			}
		}

		// Bottom region is not modified
		TriCell bottom = *this;
		bottom.maxRad = midRad;
		bottom.cellType = CT::SG;
		bottom.code = code + std::string(1, num++);
		toReturn.push_back(bottom);
	}
	return toReturn;
}


// Create renderable for given TriCell and put it into the provided renderable object
void TriCell::fillRenderable(Renderable& r, const glm::vec3 & colour, bool curved) const {

	r.lineColour = colour;
	glm::dvec3 o(0.0);

	if (curved) {
		r.drawMode = GL_TRIANGLES;
		r.verts.push_back(tri.v0 * maxRad);
		r.verts.push_back(tri.v1 * maxRad);
		r.verts.push_back(tri.v2 * maxRad);

		r.verts.push_back(tri.v0 * maxRad);
		r.verts.push_back(tri.v1 * maxRad);
		r.verts.push_back(tri.v0 * minRad);

		r.verts.push_back(tri.v0 * minRad);
		r.verts.push_back(tri.v1 * minRad);
		r.verts.push_back(tri.v1 * maxRad);

		r.verts.push_back(tri.v1 * maxRad);
		r.verts.push_back(tri.v2 * maxRad);
		r.verts.push_back(tri.v1 * minRad);

		r.verts.push_back(tri.v1 * minRad);
		r.verts.push_back(tri.v2 * minRad);
		r.verts.push_back(tri.v2 * maxRad);

		r.verts.push_back(tri.v2 * maxRad);
		r.verts.push_back(tri.v0 * maxRad);
		r.verts.push_back(tri.v2 * minRad);

		r.verts.push_back(tri.v2 * minRad);
		r.verts.push_back(tri.v0 * minRad);
		r.verts.push_back(tri.v0 * maxRad);

		for (int i = 0; i < 21; i++) {
			r.colours.push_back(colour);
		}


		//Geometry::createArcR(glm::normalize(tri.v0) * RAD * maxRad, glm::normalize(tri.v1) * RAD * maxRad, o, r);
		//Geometry::createArcR(glm::normalize(tri.v0) * RAD * maxRad, glm::normalize(tri.v2) * RAD * maxRad, o, r);
		//Geometry::createArcR(glm::normalize(tri.v1) * RAD * maxRad, glm::normalize(tri.v2) * RAD * maxRad, o, r);

		//Geometry::createArcR(glm::normalize(tri.v0) * RAD * minRad, glm::normalize(tri.v1) * RAD * minRad, o, r);
		//Geometry::createArcR(glm::normalize(tri.v0) * RAD * minRad, glm::normalize(tri.v2) * RAD * minRad, o, r);
		//Geometry::createArcR(glm::normalize(tri.v1) * RAD * minRad, glm::normalize(tri.v2) * RAD * minRad, o, r);

		//Geometry::createLineR(glm::normalize(tri.v0) * RAD * maxRad, glm::normalize(tri.v0) * RAD * minRad, r);
		//Geometry::createLineR(glm::normalize(tri.v1) * RAD * maxRad, glm::normalize(tri.v1) * RAD * minRad, r);
		//Geometry::createLineR(glm::normalize(tri.v2) * RAD * maxRad, glm::normalize(tri.v2) * RAD * minRad, r);
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


double TriCell::volume() const {

	std::vector<SphCoord> points;
	points.push_back(SphCoord(tri.v0));
	points.push_back(SphCoord(tri.v1));
	points.push_back(SphCoord(tri.v2));

	//return SphCoord::volumeCell(points, maxRad, minRad);

	double target = 4 * M_PI / 20;
	double thisV = target / pow(4, tri.sl);
	return thisV * (pow(maxRad, 3) - pow(minRad, 3)) / 3.0;
}


// Dummy constructor for testing
TestGrid::TestGrid(int n) {

	Tri t(n);
	TriCell cell;
	cell.tri = t;
	cell.maxRad = 1.0;
	cell.minRad = 0.0;
	cell.cellType = CT::SG;
	cell.code = "a";

	map[cell .code] = cell;
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
	//curdepth = 1;
}


// Subdivide all cells in the grid once
void TestGrid::subdivide(bool volume) {

	for (const auto& p : map) {

		if (p.first.length() != curDepth) continue;

		std::vector<TriCell> children = p.second.subdivide(volume);
		for (const TriCell& c : children) {
			map[c.code] = c;
		}
	}
	curDepth++;
}


enum {
	NONE = 0,
	ABOVE = 1,
	BELOW = -1
};


