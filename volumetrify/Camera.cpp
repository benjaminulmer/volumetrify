#define _USE_MATH_DEFINES
#include "Camera.h"

#include <glm/gtx/rotate_vector.hpp>

#include <cmath>

#include "Constants.h"

Camera::Camera() : zoomScale(1.3), rotScale(0.008) {
	reset();
}
#include <iostream>
// Returns view matrix for the camera
glm::dmat4 Camera::getLookAt() const {

	std::cout << "Scale: " << curScale << "\tLong: " << longitudeRotRad << "\tLat: " << latitudeRotRad << "\tTrans: " << translation.x << ", " << translation.y << ", " << translation.z << std::endl;

	// Rotate eye along longitude
	glm::dvec3 eyeTemp = glm::rotateY(eye, -longitudeRotRad);

	// Find axis then rotate eye and up along latitude
	glm::dvec3 axis = glm::cross(eyeTemp, glm::dvec3(0.0, 1.0, 0.0));

	eyeTemp = glm::rotate(eyeTemp, latitudeRotRad, axis);
	glm::dvec3 upTemp = glm::rotate(up, latitudeRotRad, axis);

	return glm::lookAt(eyeTemp + translation, centre + translation, upTemp);
}

// Returns position of the camera
glm::dvec3 Camera::getPosition() const {

	glm::dvec3 eyeTemp = glm::rotateY(eye, -longitudeRotRad);
	eyeTemp = glm::rotate(eyeTemp, latitudeRotRad, glm::cross(eyeTemp, glm::dvec3(0.0, 1.0, 0.0)));

	return eyeTemp + translation;
}

// Returns up of the camera
glm::dvec3 Camera::getUp() const {

	// Rotate eye along longitude
	glm::dvec3 eyeTemp = glm::rotateY(eye, -longitudeRotRad);

	// Find axis then rotate eye and up along latitude
	glm::dvec3 axis = glm::cross(eyeTemp, glm::dvec3(0.0, 1.0, 0.0));

	eyeTemp = glm::rotate(eyeTemp, latitudeRotRad, axis);
	return glm::rotate(up, latitudeRotRad, axis);
}

// Returns looking direction of camera
glm::dvec3 Camera::getLookDir() const {
	return glm::normalize(centre - getPosition());
}

// Sets current model scale
void Camera::setScale(double scale) {
	curScale = scale;
}

// Rotates camera along longitudinal axis (spherical coords)
void Camera::updateLongitudeRotation(double rad) {

	// If camera is upside down reverse longitude rotations
	if (cos(latitudeRotRad) > 0) {
		longitudeRotRad += rad * M_PI / 180;
	}
	else {
		longitudeRotRad -= rad * M_PI / 180;
	}
}

// Rotates camera along latitudinal axis (spherical coords)
void Camera::updateLatitudeRotation(double rad) {
	latitudeRotRad += rad * M_PI / 180;
}

// Zooms camera in or out (+1 or -1)
void Camera::updateZoom(int sign) {

	if (sign < 0) {
		eye.z = (eye.z - RADIUS_EARTH_VIEW) / zoomScale + RADIUS_EARTH_VIEW;
	}
	else {
		eye.z = (eye.z - RADIUS_EARTH_VIEW) * zoomScale + RADIUS_EARTH_VIEW;;
	}
	if (eye.z > 4.0 * RADIUS_EARTH_VIEW) eye.z = 4.0 * RADIUS_EARTH_VIEW;
}

// Translates camera along x and y of view plane
void Camera::translate(const glm::dvec3& planeTranslation) {

	glm::dvec3 pTrans(planeTranslation);

	// Scale translation based on zoom level
	double scale = (eye.z / 500.f);
	pTrans *= scale;

	// Get rotation axis
	glm::dvec3 eyeTemp = glm::rotateY(eye, -longitudeRotRad);
	glm::dvec3 axis = glm::cross(eyeTemp, glm::dvec3(0.0, 1.0, 0.0));

	// Convert screen space translation into world space translation
	glm::dvec3 rotatedTranslation = glm::rotateY(planeTranslation, -longitudeRotRad);
	rotatedTranslation = glm::rotate(rotatedTranslation, latitudeRotRad, axis);

	translation += rotatedTranslation;
}

// Reset camera to starting position
void Camera::reset() {
	eye = glm::dvec3(0.0, 0.0, RADIUS_EARTH_VIEW + 30.0);
	up = glm::dvec3(0.0, 1.0, 0.0);
	centre = glm::dvec3(0.0, 0.0, 0.0);

	longitudeRotRad = 1.62316;
	latitudeRotRad = 0.0;
	translation = glm::dvec3(2.47295, 6.70685, 5.38065);

	longitudeRotRad = -0.593408;
	latitudeRotRad = 0.863938;
	translation = glm::dvec3(2.9946, 7.7957, 3.83579);

	longitudeRotRad = -1.05592;
	latitudeRotRad = 0.253073;
	translation = glm::dvec3(-1.17255, 5.42163, -0.794411);

	longitudeRotRad = -1.065465;
	latitudeRotRad = 0.436333;
	translation = glm::dvec3(-1.72643, 6.78067, -1.10433);

	longitudeRotRad = 1.0472;
	latitudeRotRad = 1.570796;
	translation = glm::dvec3(1.4, 0.0, -0.8);

	// 200 skewed
	longitudeRotRad = 0.820304;
	latitudeRotRad = 0.654498;
	translation = glm::dvec3(2.40908, 6.76859, -2.30227);

	longitudeRotRad = 0.322886;
	latitudeRotRad = 0.558505;
	translation = glm::dvec3(0.578206, 6.45964, -3.31619);

	longitudeRotRad = 0.218166;
	latitudeRotRad = 0.383972;
	translation = glm::dvec3(0.216706, 5.39621, -3.07148);

	// 20 skewed
	longitudeRotRad = 0.733038;
	latitudeRotRad = 0.427605;
	translation = glm::dvec3(1.71190, 5.94776, -2.26126);

	// 200 side 
	longitudeRotRad = -1.0472;
	latitudeRotRad = -0.253072;
	translation = glm::dvec3(0.0, 6.17, 0.0);

	// 20 side 
	longitudeRotRad = -1.0472;
	latitudeRotRad = -0.261799;
	translation = glm::dvec3(0.0, 6.46, 0.0);
}