#include "stdafx.h"
#include "Theory.h"
#include "MCSim_application.h"
#include <math.h>

double Theory::CalculateMagnetizationForSingleParticle(double field)
{
	return 1 / tanh(field) - 1 / field;
}

double Theory::CalculateFieldMultiplier()
{
	double l = MCSim_application::GetLambda();
	double p = MCSim_application::GetTargetVolumeDensity();
	double m = MCSim_application::GetParticleMagneticMoment();

	return 4 * l * p / m;
}

double Theory::CalculateMagnetizationForSystem(double field)
{
	double m = MCSim_application::GetParticleMagneticMoment();
	double n = MCSim_application::GetTargetVolumeDensity();
	double M =
		CalculateMagnetizationForSingleParticle(field)
		* m
		* n;

	double H = field / CalculateFieldMultiplier();
	double effectiveM = H + M / 3;

	return m * n * CalculateMagnetizationForSingleParticle(effectiveM);
}
