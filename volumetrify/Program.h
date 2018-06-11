#pragma once

#include <iostream>

#include <SDL2/SDL.h>
#undef main

#include "Camera.h"
#include "RenderEngine.h"

enum class RadialBound {
	MAX,
	MIN,
	BOTH
};

class Program {

public:
	Program();

	void start();

	void updateRotation(int oldX, int newX, int oldY, int newY, bool skew);
	void updateScale(int inc);

private:
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
	void insertAirSigmets();
	void insertWind();

	void airSigRender1();
	void windRender1();

	void mainLoop();
};
