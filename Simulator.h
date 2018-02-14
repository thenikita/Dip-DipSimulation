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
	double tubeR, tempR;
	double tubeL, tempL;
	double muliplyer = 100;
	std::mt19937_64 generator;
	double deltaCoordinate = 0.1;
	int stepsAmount = 10;
	double particleDiameter = 1;
	int particleAmount;
public:
	Simulator(
		double tubeR, 
		double tubeL, 
		int particleAmount);

	std::vector<Particle> getParticles() { return particles; }
	void MakeIterations(int particleAmount, bool ifNeedResize);

	double GetTubeR() { return tubeR; }
	double GetTubeL() { return tubeL; }

	double SetTubeR(double R) { this->tubeR = R; }
	double SetTubeL(double L) { this->tubeL = L; }
private:
	void GenerateParticles(int particleAmount);
	Particle GenerateDeltaState(const Particle particle);
	double CalculateParticleEnergy(const Particle particle);
	bool CheckParticleForCollisions(const Particle particle,int ignored);
	static double GenerateRandom(double min, double max, std::mt19937_64 &generator);
	void SetGeneratotRandomSeed();
	std::vector<double> GetVectorField(std::vector<double> point, bool mode);
	bool ChechSystemForErrors();
	void CollectTubeSizes(int num);
	void ResizeTube();
	bool CheckIfNeedResize();
};

