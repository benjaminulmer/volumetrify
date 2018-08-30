#define _USE_MATH_DEFINES
#include "Program.h"

#include <GL/glew.h>
#include <glm/gtx/intersect.hpp>
#include <glm/gtx/transform.hpp>
#include <SDL2/SDL_opengl.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>

#include "Constants.h"
#include "ContentReadWrite.h"
#include "InputHandler.h"
#include "SphCoord.h"
#include "Geometry.h"
#define _N 200
Program::Program() : g(_N), v(_N) {

	window = nullptr;
	renderEngine = nullptr;
	camera = nullptr;

	longRot = 0;
	latRot = 0;

	curved = false;
	volume = false;
	fullSphere = false;
	sl = 0;

	width = height = 800;
}


// Called to start the program. Conducts set up then enters the main loop
void Program::start() {	

	//DGGS_3D dg(true, 18, 1.0, 100.0);
	//std::cout.precision(17);
	//std::cout << dg.layer(50.0, 25) << std::endl;

	//for (int i = 0; i < dg.layerSurfaceSL.size(); i++) {
	//	std::cout << i << ": " << dg.layerSurfaceSL[i] << std::endl;
	//}

	//std::cout << "Lyr\tDp\tBnd\tSG\tSGC\tEnd" << std::endl << "--------------------------------------------" << std::endl;
	//int n = 2;
	//for (int j = 0; j <= 36  ; j++) {

	//	LayerInfo i = dg.layerInfo(j);
	//	std::cout << j << ":\t"<< i.depth << "\t" << i.boundary << "\t" << i.degen << "\t" << i.degenChild << "\t" << i.end << std::endl;
	//}
	//system("pause");


	setupWindow();
	GLenum err = glewInit();
	if (glewInit() != GLEW_OK) {
		std::cerr << glewGetErrorString(err) << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	camera = new Camera();
	renderEngine = new RenderEngine(window);
	InputHandler::setUp(camera, renderEngine, this);

	// Assign buffers
	RenderEngine::assignBuffers(grid, false);
	RenderEngine::assignBuffers(ref, false);

	// Set starting radius
	scale = 1.0;
	radius = RADIUS_EARTH_M * 4.0 / 3.0;

	// Objects to draw initially
	objects.push_back(&grid);

	// Draw stuff
	grid.drawMode = GL_LINES;

	// Subdivide grid and create renderable
	for (int i = 0; i < maxSL; i++) {
		g.subdivide(false);
		v.subdivide(true);
	}
	updateRenderable();

	Tri start(_N);
	ref.verts.push_back(start.v0);
	ref.verts.push_back(start.v1);
	ref.verts.push_back(start.v2);

	ref.verts.push_back(start.v0);
	ref.verts.push_back(start.v1);
	ref.verts.push_back(glm::vec3(0.f));

	ref.verts.push_back(start.v1);
	ref.verts.push_back(start.v2);
	ref.verts.push_back(glm::vec3(0.f));

	ref.verts.push_back(start.v2);
	ref.verts.push_back(start.v0);
	ref.verts.push_back(glm::vec3(0.f));

	for (int i = 0; i < 12; i++) {
		ref.colours.push_back(glm::vec3(0.8));
	}
	ref.doubleToFloats();
	RenderEngine::setBufferData(ref, false);
	ref.drawMode = GL_TRIANGLES;
	ref.model = glm::translate(glm::vec3(0.f, 0.5f * 6371000.0 * (4.0 / 3.0), 0.f)) * glm::scale(glm::vec3(0.99f, 0.999f, 0.99f)) * glm::translate(glm::vec3(0.f, -0.5f * 6371000.0 * (4.0 / 3.0), 0.f));

	objects.push_back(&ref);
	mainLoop();
}


void Program::updateRenderable() {

	grid.clear();

	const auto& map = (volume) ? v.map : g.map;

	//double max = -1.0;
	//double min = 999999.0;
	//for (auto p : map) {
	//	if (p.first.length() != sl + 1 || (p.first[0] != 'a' && !fullSphere)) continue;

	//	double v = p.second.volume();
	//	max = (v > max) ? v : max;
	//	min = (v < min) ? v : min;
	//}

	//for (auto p : map) {
	//	if (p.first.length() != sl + 1 || (p.first[0] != 'a' && !fullSphere)) continue;

	//	double v = p.second.volume();
	//	glm::vec3 colour;

	//	// If all volumes the same return grey
	//	if (max == min) {
	//		colour = glm::vec3(0.5f, 0.5f, 0.5f);
	//	}

	//	// Set colour
	//	if ((v / max) > 0.999f && (v / max) < 1.001f) {
	//		colour = glm::vec3(1.f, 0.f, 0.f);
	//	}
	//	else if ((v / min) > 0.999f && (v / min) < 1.001f) {
	//		colour = glm::vec3(0.f, 0.f, 1.f);
	//	}
	//	else {
	//		float norm = (v - min) / (max - min);
	//		colour = glm::vec3(norm, norm, norm);
	//	}

	//	p.second.fillRenderable(grid, colour, curved);
	//}


	double mean = 0.0;
	int n = 0;
	for (auto p : map) {
		if (p.first.length() != sl + 1 || (p.first[0] != 'a' && !fullSphere)) continue;
		mean += p.second.volume();
		n++;
	}
	mean /= n;

	double sd = 0.0;
	for (auto p : map) {
		if (p.first.length() != sl + 1 || (p.first[0] != 'a' && !fullSphere)) continue;
		double v = p.second.volume();
		sd += (v - mean) * (v - mean);
	}
	sd = sqrt(sd / n);

	for (auto p : map) {
		if (p.first.length() != sl + 1 || (p.first[0] != 'a' && !fullSphere)) continue;
		double v = p.second.volume();
		glm::vec3 colour;

		double z = (v - mean) / sd;
		float sat = 0.f;
		float light = 0.5f;
		int hue = 0;

		if (z > 0.0) {
			hue = 0;
			sat = 0.3 * z;
			float c = (1.f - abs(2.f * light - 1.f)) * sat;
			int h = hue / 60;
			float x = c * (1.f - abs(h % 2 - 1));
			float m = light - 0.5f * c;
			colour = glm::vec3(c + m, x + m, m);
		}
		else {
			hue = 250;
			sat = 0.3 * -z;
			float c = (1.f - abs(2.f * light - 1.f)) * sat;
			int h = hue / 60;
			float x = c * (1.f - abs(h % 2 - 1));
			float m = light - 0.5f * c;
			colour = glm::vec3(x + m, m, c + m);
		}
		p.second.fillRenderable(grid, glm::vec3(0.f), curved);
	}
	std::cout << mean << " : " << sd << std::endl;
	grid.doubleToFloats();
	RenderEngine::setBufferData(grid, false);
}


// Creates SDL window for the program and sets callbacks for input
void Program::setupWindow() {

	if (SDL_Init(SDL_INIT_VIDEO) != 0){
		std::cerr << "SDL_Init Error: " << SDL_GetError() << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 16);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

	window = SDL_CreateWindow("sudivision framework", 10, 30, width, height, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
	if (window == nullptr) {
		//TODO: cleanup methods upon exit
		std::cerr << "SDL_CreateWindow Error: " << SDL_GetError() << std::endl;
		SDL_Quit();
		exit(EXIT_FAILURE);
	}

	SDL_GLContext context = SDL_GL_CreateContext(window);
	if (context == NULL) {
		std::cout << "OpenGL context could not be created! SDL Error: " << SDL_GetError() << std::endl;
	}
	SDL_GL_SetSwapInterval(1); // Vsync on
}


// Main loop
void Program::mainLoop() {

	while (true) {

		// Process all SDL events
		SDL_Event e;
		while (SDL_PollEvent(&e)) {
			InputHandler::pollEvent(e);
		}

		glm::dmat4 worldModel(1.f);
		double s = scale * (1.0 / RADIUS_EARTH_M) * RADIUS_EARTH_VIEW;
		worldModel = glm::scale(worldModel, glm::dvec3(s, s, s));
		worldModel = glm::rotate(worldModel, latRot, glm::dvec3(-1.0, 0.0, 0.0));
		worldModel = glm::rotate(worldModel, longRot, glm::dvec3(0.0, 1.0, 0.0));

		// Find min and max distance from camera to cell renderable - used for fading effect
		glm::vec3 cameraPos = camera->getPosition();
		float max = -1.f;
		float min = 999999999.f;

		for (const glm::vec3& vert : grid.verts) {
			float l = glm::length((float)s * vert - cameraPos);
			max = std::max(max, l);
			min = std::min(min, l);
		}
		//std::cout << max << " : " << min << std::endl;

		renderEngine->render(objects, (glm::dmat4)camera->getLookAt() * worldModel, max, min);
		SDL_GL_SwapWindow(window);
	}
}


// Updates camera rotation
// Locations are in pixel coordinates
void Program::updateRotation(int oldX, int newX, int oldY, int newY, bool skew) {

	glm::dmat4 projView = renderEngine->getProjection() * camera->getLookAt();
	glm::dmat4 invProjView = glm::inverse(projView);

	double oldXN = (2.0 * oldX) / (width) - 1.0; 
	double oldYN = (2.0 * oldY) / (height) - 1.0;
	oldYN *= -1.0;

	double newXN = (2.0 * newX) / (width) - 1.0;
	double newYN = (2.0 * newY) / (height) - 1.0;
	newYN *= -1.0;

	glm::dvec4 worldOld(oldXN, oldYN, -1.0, 1.0);
	glm::dvec4 worldNew(newXN, newYN, -1.0, 1.0);

	worldOld = invProjView * worldOld; 
	worldOld /= worldOld.w;

	worldNew = invProjView * worldNew;
	worldNew /= worldNew.w;

	glm::dvec3 rayO = camera->getPosition();
	glm::dvec3 rayDOld = glm::normalize(glm::dvec3(worldOld) - rayO);
	glm::dvec3 rayDNew = glm::normalize(glm::dvec3(worldNew) - rayO);
	double sphereRad = RADIUS_EARTH_VIEW * scale;
	glm::dvec3 sphereO = glm::dvec3(0.0);

	glm::dvec3 iPosOld, iPosNew, iNorm;

	if (glm::intersectRaySphere(rayO, rayDOld, sphereO, sphereRad, iPosOld, iNorm) && 
			glm::intersectRaySphere(rayO, rayDNew, sphereO, sphereRad, iPosNew, iNorm)) {

		double longOld = atan2(iPosOld.x, iPosOld.z);
		double latOld = M_PI_2 - acos(iPosOld.y / sphereRad);

		double longNew = atan2(iPosNew.x, iPosNew.z);
		double latNew = M_PI_2 - acos(iPosNew.y / sphereRad);

		if (skew) {
			camera->updateLatitudeRotation(latNew - latOld);
		}
		else {
			latRot += latNew - latOld;
			longRot += longNew - longOld;
		}
	}
}


// Changes scale of model
void Program::updateScale(int inc) {

	if (inc < 0) {
		scale /= 1.4f;
	}
	else {
		scale *= 1.4f;
	}
	camera->setScale(scale);
}


void Program::toggleCurved() {
	curved = !curved;
	updateRenderable();
}


void Program::toggleVolume() {
	volume = !volume;
	updateRenderable();
}


void Program::toggleFullSphere() {
	fullSphere = !fullSphere;
	updateRenderable();
}


void Program::updateSL(int inc) {
	sl += inc;
	sl = std:: min(sl, maxSL);
	sl = std::max(sl, 0);
	updateRenderable();
}