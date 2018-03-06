#include "Simulator.h"

#pragma once
class MCSim_application
{
private:
	static double lambda; 
	static double field; //xi module. vector is in simulator obj
	static double targetVolumeDensity; //phi
	static double particleMagneticMoment;

	static int particleAmount;
	static double particleDiameter;
	static double aspect;
	static double targetTubeR, targetTubeL;
public:
	MCSim_application();
	//static double GetParticleDiameter() { return particleDiameter; }
	static double GetLambda() { return lambda; }
	static double GetField() { return field; }
	static double GetTargetVolumeDensity() { return targetVolumeDensity; }
	static double GetTargetRadius() { return targetTubeR; }
	static double GetTubeAspect() { return aspect; }
	static double GetParticleMagneticMoment() { return particleMagneticMoment; }

	static void setField(double module) { field = module; }

	static void showSystem(Simulator simulator);
private:
	void GenerateTube();
};
