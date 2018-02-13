#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <thread>

#include "Simulator.h"
#include "Particle.h"
#include "MCSim_application.h"
#include "Theory.h"

using std::cin;
using std::cout;
using std::endl;

using std::mt19937_64;
using std::uniform_real_distribution;


//TODO: chech if tube size less than D of particle

Simulator::Simulator(
	double tubeR,
	double tubeL,
	int particleAmount)
{
	cout << "Simulator is launched!" << endl;

	//this->particleAmount = particleAmount;

	this->tubeR = tubeR * muliplyer;	// Here we assigned target values of tube to current values.
	this->tubeL = tubeL * muliplyer;	// Below we will increase this ones because of method of getting
										// target volume density.
	
	cout << "Enter the field (xi) vector (x, y, z): " << endl;
	
	double tempField;
	for (int i = 0; i < 3; i++)
	{
		cin >> tempField;
		field.push_back(tempField);
	}

	cout << "Field is set to: (";
	cout << field.at(0) << ", " << field.at(1) << ", " << field.at(2) << ")" << endl;

	MCSim_application::setField(Particle::calculateVectorModule(field));
	cout << "Field module: " << MCSim_application::getField() << endl;

	SetGeneratotRandomSeed();

	GenerateParticles(particleAmount);

	ChechSystemForErrors();
}


//TODO: found a bug here, it generates bad
//TODO: make a thing whick checks if tubeR > particle's R
void Simulator::GenerateParticles(int particleAmount)
{
	cout << "Generating " << particleAmount << " particles to " << std::flush;
	cout << "TubeR: " << tubeR << ". TubeL: " << tubeL << "." << endl;

	//this->particleAmount = particleAmount;

	//Particle tempPart(0, 0, 0, 1, 0, 0);
	//particles.push_back(tempPart);

	for (int i = 0; i < particleAmount; i++)
	{
		double x, y, z;
		double mx, my, mz;
		double R = tubeR;

		//TODO: generate more than one particle
		x = GenerateRandom(-R + 1, R - 1, generator);
		R = sqrt(R * R - x * x);
		y = GenerateRandom(-R + 1, R - 1, generator);
		z = GenerateRandom(-tubeL / 2 + 1, tubeL / 2 - 1, generator);

		Particle temp(x, y, z, 1, 0, 0);
		if (!CheckParticleForCollisions(temp, 100000))
		{
			this->particleAmount++;
			particles.push_back(temp);
			cout << "Particles ready: " << i + 1 << "\r" << std::flush;
		}
		else
		{
			i--;
		}
	}

	cout << "\nParticles've been generated!" << endl;
}


//
// IMPORTANT MEMO
// Here we calculate energy predivided by kT so later we can
// use it as is in transition calculations without more dividings.
// 

Particle Simulator::GenerateDeltaState(const Particle particle)
{
	double newX = particle.x + GenerateRandom(-1, 1, generator) * deltaCoordinate;
	double newY = particle.y + GenerateRandom(-1, 1, generator) * deltaCoordinate;
	double newZ = particle.z + GenerateRandom(-1, 1, generator) * deltaCoordinate;
	
	double left = 1;

	double newMZ = GenerateRandom(-1, 1, generator);
	left -= newMZ * newMZ;
	double newMY = GenerateRandom(-sqrt(left), sqrt(left), generator);
	left -= newMY * newMY;
	//TODO: check whether we need +- x generation
	double newMX = sqrt(left);

	Particle newState(  newX,
						newY,
						newZ,
						newMX,
						newMY,
						newMZ);
	return newState;
}


//TODO: make iterations
//TODO: add deltaCoordinate decreasing in future iterations to immprove speed
void Simulator::MakeIterations(int particleAmount)
{
	//this->particleAmount = particleAmount;

	cout << "Iteration process is starting..." << endl;

	double averageProjection = 0;

	for (int i = 0; i < stepsAmount; i++)
	{
		cout << "STEP #" << i << " Tube: " << tubeR << " " << tubeL << "              \r" << std::flush;

		for (int j = 0; j < particleAmount; j++)
		{
			Particle temp = GenerateDeltaState(particles.at(j));

			if (!CheckParticleForCollisions(temp, j))
			{
				double energy = CalculateParticleEnergy(temp) - CalculateParticleEnergy(particles.at(j));

				if (exp(energy) > GenerateRandom(0, 1, generator))
				{
					particles.at(j) = temp;
				}

				averageProjection += particles.at(j).mz;
			}
		}
	}
	cout << "\nERRORS: " << particles.at(0).CheckForErrors(tubeR, tubeL) << endl;

	cout << "\nIterations've finished!" << endl;
	cout << "Average is: " << averageProjection / stepsAmount / particleAmount << endl;
	cout << "Theory is:  " << Theory::CalculateForSingleParticle(MCSim_application::getField()) << endl;

	ChechSystemForErrors();
}


double Simulator::CalculateParticleEnergy(const Particle particle)
{
	//here's the field part of energy only!
	double outerField = MCSim_application::getField();

	std::vector<double> moment;
	moment.push_back(particle.mx);
	moment.push_back(particle.my);
	moment.push_back(particle.mz);

	double cosTheta = Particle::calculateCosinus(moment, field);
	//cout << "cos: " << cosTheta << endl;

	double outerFieldEnergy = outerField * cosTheta;

	double dipoleInteractionEnergy = 0;
	//TODO: calculate dipole interaction energy
	//Already have dip-dip function in particle obj

	double totalEnergy = dipoleInteractionEnergy + outerFieldEnergy;

	return totalEnergy;
}


//TODO: IMPORTANT checking for collisions: /walls /particles
bool Simulator::CheckParticleForCollisions(const Particle particle, int ignored)
{
	//w walls
	double R =
		particle.x * particle.x +
		particle.y * particle.y;

	if ((sqrt(R) + particleDiameter / 2) > tubeR) return true;
	if (abs(particle.z) + particleDiameter / 2 > tubeL / 2) return true;

	//w parts
	for (int i = 0; i < this->particleAmount; i++)
	{
		if ((i != ignored))
		{
			double distance =
				(particles.at(i).x - particle.x) * (particles.at(i).x - particle.x) +
				(particles.at(i).y - particle.y) * (particles.at(i).y - particle.y) +
				(particles.at(i).z - particle.z) * (particles.at(i).z - particle.z);

			if (distance < particleDiameter * particleDiameter) return true;
		}
	}

	return false;
}


double Simulator::GenerateRandom(double min, double max, mt19937_64 &generator)
{
	uniform_real_distribution<double> distr(min, max);
	return distr(generator);
}


void Simulator::SetGeneratotRandomSeed()
{
	generator.seed(1);
	std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
	
	int timePast = int(GenerateRandom(0, 1000000000, generator));
	std::this_thread::sleep_for(std::chrono::nanoseconds(timePast));
	std::chrono::steady_clock::time_point then = std::chrono::steady_clock::now();
	unsigned seed = std::chrono::duration_cast<std::chrono::nanoseconds>(then - now).count();
	cout << "RANDOM SEED IS: " << seed << " In time: " << timePast << "nsecs." << endl;
	generator.seed(seed);
}


std::vector<double> Simulator::GetVectorField(std::vector<double> point, bool mode)
{
	// true = normal to z oriented field
	// false = center oriented field
	if (mode == true)
	{
		return field;
	}
	else
	{
		std::vector<double> temp;
		temp.push_back(-point.at(0));
		temp.push_back(-point.at(1)); 
		temp.push_back(-point.at(2));

		return temp;
	}
}


bool Simulator::ChechSystemForErrors()
{
	cout << "Checking system for errors in " << this->particleAmount << endl;
	for (int i = 0; i < this->particleAmount; i++)
	{
		//cout << i << std::flush;
		for (int j = i + 1; j < this->particleAmount; j++)
		{
			//cout << j << std::flush;

			double distance =
				pow(particles.at(i).x - particles.at(j).x, 2) +
				pow(particles.at(i).y - particles.at(j).y, 2) +
				pow(particles.at(i).z - particles.at(j).z, 2);

			if (distance < particleDiameter * particleDiameter)
			{
				cout << "\nERROR on " << i << ", " << j << endl;
				return true;
			}
		}
		//cout << endl;
	}
	cout << "NO ERRORS" << endl;
	return false;
}
