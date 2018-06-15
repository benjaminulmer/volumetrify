#pragma once

#include <iostream>

#include <SDL2/SDL.h>
#undef main

#include "Camera.h"
#include "RenderEngine.h"
#include "TestGrid.h"

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

	void toggleCurved();
	void toggleVolume();
	void updateSL(int inc);

private:
	TestGrid g;
	TestGrid v;

	const int maxSL = 6;
	bool curved;
	bool volume;
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
