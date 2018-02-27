#include "attitude_controller.h"
#include <iostream>

int main() {
	float input[1];

	Ballbot::AttitudeController controller;
	std::cout << "controller1" << std::endl << "=======" << std::endl;
	controller.load(Ballbot::AttitudeController::CONTROLLER1);
	input[0] = 1.0f/controller.Ts();
	controller.update(input);

	input[0] = 0.0;
	for(int k=0;k<50;k++){
		std::cout << controller.actuation()[0] << std::endl;
		controller.update(input);
	}

	std::cout << "controller2" << std::endl << "=======" << std::endl;
	controller.load(Ballbot::AttitudeController::CONTROLLER2);
	input[0] = 1.0f/controller.Ts();
	controller.update(input);

	input[0] = 0.0;
	for(int k=0;k<50;k++){
		std::cout << controller.actuation()[0] << std::endl;
		controller.update(input);
	}

}