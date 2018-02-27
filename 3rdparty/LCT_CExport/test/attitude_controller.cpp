#include "attitude_controller.h"

Ballbot::AttitudeController::AttitudeController() :
	_Ts(1.000000e-02),
	_algorithm(-1)
{
	//do nothing
}

void Ballbot::AttitudeController::load(int algorithm) {
	switch(algorithm) {
		case CONTROLLER1:{
			load_controller1();
			break;}

		case CONTROLLER2:{
			load_controller2();
			break;}

		default:{
			return;}
	}

	_algorithm = algorithm;
	reset();
}

float* Ballbot::AttitudeController::update(float* measurements) {
	switch(_algorithm) {
		case CONTROLLER1:{
			update_controller1(measurements);
			break;}

		case CONTROLLER2:{
			update_controller2(measurements);
			break;}

	}

	return _y;
}

float* Ballbot::AttitudeController::actuation() {
	return _y;
}

float Ballbot::AttitudeController::Ts() {
	return _Ts;
}

void Ballbot::AttitudeController::reset() {
	for(int k = 0; k<3; k++)
		_x[k] = 0.0f;

	for(int k = 0; k<1; k++)
		_y[k] = 0.0f;
}

// private inline functions
// controller1 functions
inline void Ballbot::AttitudeController::load_controller1() {
	//do nothing for LTI here
}

inline void Ballbot::AttitudeController::update_controller1(float *measurements) {
	float x[2];
	for(int k = 0; k<2; k++)
		x[k] = _x[k];

	_y[0] = -1.37900e-01f*x[0] +6.19900e-01f*x[1] ;
	_x[0] = -1.91600e-01f*x[0] -1.81400e-01f*x[1] ;
	_x[1] = +1.81400e-01f*x[0] -1.91600e-01f*x[1] +5.44700e-01f*measurements[0] ;
}

// controller2 functions
inline void Ballbot::AttitudeController::load_controller2() {
	//do nothing for LTI here
}

inline void Ballbot::AttitudeController::update_controller2(float *measurements) {
	float x[3];
	for(int k = 0; k<3; k++)
		x[k] = _x[k];

	_y[0] = -1.66600e+00f*x[0] -4.15900e-01f*x[1] -8.42200e-02f*x[2] ;
	_x[0] = -3.57300e-01f*x[0] +2.65900e-01f*x[1] -3.17200e-01f*x[2] ;
	_x[1] = +2.52300e-01f*x[0] -9.95200e-02f*x[1] -5.24300e-01f*x[2] +6.92200e-02f*measurements[0] ;
	_x[2] = -3.28100e-01f*x[0] -5.17500e-01f*x[1] +1.31000e-01f*x[2] +2.48700e+00f*measurements[0] ;
}

