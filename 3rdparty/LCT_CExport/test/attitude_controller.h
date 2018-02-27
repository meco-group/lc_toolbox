#ifndef ATTITUDECONTROLLER_H
#define ATTITUDECONTROLLER_H

namespace Ballbot {

class AttitudeController
{
public:
	typedef enum algorithm_t {
		CONTROLLER1,
		CONTROLLER2,
	} algorithm_t;

private:
	const float _Ts;
	int _algorithm;

	float _x[3];
	float _y[1];

	//controller1 functions
	inline void load_controller1();
	inline void update_controller1(float* measurements);

	//controller2 functions
	inline void load_controller2();
	inline void update_controller2(float* measurements);

public:
	AttitudeController();

	void load(int algorithm);
	float* update(float* measurements);

	float* actuation();
	float Ts();
	void reset();
}; //class

}; //namespace

#endif //ATTITUDECONTROLLER_H
