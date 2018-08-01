#pragma once

#include <iostream>

#include <SDL2/SDL.h>
#undef main

#include "Camera.h"
#include "RenderEngine.h"
#include "TestGrid.h"
#include "DGGS_3D.h"

enum class RadialBound {
	MAX,
	MIN,
	BOTH
};

class Program {

public:
	Program();

	void start();

	void updateRenderable();

	void updateRotation(int oldX, int newX, int oldY, int newY, bool skew);
	void updateScale(int inc);

	// VOLUME
	void toggleCurved();
	void toggleVolume();
	void toggleFullSphere();
	void updateSL(int inc);

private:
	TestGrid g;
	TestGrid v;

	const int maxSL = 6;
	bool curved;
	bool volume;
	bool fullSphere;
	int sl;

	SDL_Window* window;
	int width, height;

	RenderEngine* renderEngine;
	Camera* camera;

	Renderable grid;

	std::vector<const Renderable*> objects;

	double radius;

	double scale;
	double latRot;
	double longRot;

	void setupWindow();

	void mainLoop();
};
