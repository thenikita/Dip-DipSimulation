//#include "MCSim_application.h"
#pragma once
#include <vector>
#include <string>
class Particle
{
public:
	Particle(
		double x,
		double y,
		double z,
		double mx,
		double my,
		double mz);
	static double CalculateDipoleEnergy(double lambda, Particle first, Particle second);
	static double CalculateInFieldEnergy(double field);
	static double ScalarProduction(std::vector<double> first, std::vector<double> second);
public:
	double x, y, z;
	double mx, my, mz;
	//static const double d = MCSim_application::GetParticleDiameter();

public:
	std::string toString();
	static double calculateCosinus(std::vector<double> first, std::vector<double> second);
	static double calculateVectorModule(std::vector<double> first);
	bool CheckForErrors(double R, double L);
};