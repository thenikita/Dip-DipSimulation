#include "Simulator.h"

#pragma once
class MCSim_application
{
private:
	static double lambda; 
	static double field; //xi module. vector is in simulator obj
	static double targetVolumeDensity; //phi

	static int particleAmount;
	static double particleDiameter;
	static double aspect;
	static double targetTubeR, targetTubeL;
public:
	MCSim_application();
	//static double GetParticleDiameter() { return particleDiameter; }
	static double getLambda() { return lambda; }
	static double getField() { return field; }
	static double getTargetVolumeDensity() { return targetVolumeDensity; }
	static double getTargetRadius() { return targetTubeR; }
	static double getTubeAspect() { return aspect; }

	static void setField(double module) { field = module; }

	static void showSystem(Simulator simulator);
private:
	void GenerateTube();
};
