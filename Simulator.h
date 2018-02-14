#pragma once
#include <vector>
#include <random>
#include <ctime>
#include <chrono>

#include "Particle.h"

class Simulator
{
public:
private:
	std::vector<Particle> particles;
	std::vector<double> field;
	double tubeR;
	double tubeL;
	double muliplyer = 100;
	std::mt19937_64 generator;
	double deltaCoordinate = 0.1;
	int stepsAmount = 100000;
	double particleDiameter = 1;
	int particleAmount;
public:
	Simulator(
		double tubeR, 
		double tubeL, 
		int particleAmount);

	std::vector<Particle> getParticles() { return particles; }
	void MakeIterations(int particleAmount);

	double GetTubeR() { return tubeR; }
	double GetTubeL() { return tubeL; }
private:
	void GenerateParticles(int particleAmount);
	Particle GenerateDeltaState(const Particle particle);
	double CalculateParticleEnergy(const Particle particle);
	bool CheckParticleForCollisions(const Particle particle,int ignored);
	static double GenerateRandom(double min, double max, std::mt19937_64 &generator);
	void SetGeneratotRandomSeed();
	std::vector<double> GetVectorField(std::vector<double> point, bool mode);
	bool ChechSystemForErrors();
};

