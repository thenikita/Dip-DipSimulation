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

const double PI = 3.14;

Simulator::Simulator(
	double tubeR,
	double tubeL,
	int particleAmount)
{
	cout << "Simulator is launched!" << endl;

	this->particleAmount = 0;

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
	cout << "Field module: " << MCSim_application::GetField() << endl;
	cout << "\n********************************************\n";

	SetGeneratorRandomSeed();

	GenerateParticles(particleAmount);

	CheckSystemForErrors();
}


void Simulator::ShowSystem()
{
	cout << "\n################# SYSTEM ###################\n";
	cout << "Lambda:           " << MCSim_application::GetLambda() << endl;
	cout << "Field Module:     " << MCSim_application::GetField() << endl;
	cout << "Target Density:   " << MCSim_application::GetTargetVolumeDensity() << endl;
	cout << "Current Density:  " << CalculateCurrentVolumeDensity() << endl;

	cout << "Particles:        " << particleAmount << endl;
	cout << "  N  |            COORDINATES             |            MAGNETIC MOMENT       " << endl;

	for (int i = 0; i < particleAmount; i++)
	{
		cout << "  " << i << "    ";
		cout << particles.at(i).toString() << endl;
	}
	cout << "\n################# SYSTEM ###################\n";
}

void Simulator::GenerateParticles(int particleAmount)
{
	cout << "Generating " << particleAmount << " particles to " << std::flush;
	cout << "TubeR: " << tubeR << ". TubeL: " << tubeL << "." << endl;

	//Particle tempPart(0, 0, 0, 1, 0, 0);
	//particles.push_back(tempPart);

	for (int i = 0; i < particleAmount; i++)
	{
		double x, y, z;
		double mx, my, mz;
		double R = tubeR;

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


/* Generate Delta State function generates new state of i-th particle
   which is based on current position in the tube and using the
   random number generator*/
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


/* Make Iterations function is making iterations moving particles
   inside the tube. Use after resizing the tube to relevant size
   if wanna get results accounting demagnetisation and other shit*/

//TODO: add deltaCoordinate decreasing in future iterations to immprove speed
void Simulator::MakeIterations(int particleAmount, bool ifNeedResize)
{
	//this->particleAmount = particleAmount;

	cout << "\n****************ITERATIONS******************\n";
	

	double averageProjection = 0;
	
	//check it
	for (int i = 0; i < stepsAmount; i++)
	{
		cout << 100 * i / stepsAmount << " %  Tube: " << tubeR << " " << tubeL << "   \r" << std::flush;

		for (int j = 0; j < particleAmount; j++)
		{
			Particle temp = GenerateDeltaState(particles.at(j));

			if (!CheckParticleForCollisions(temp, j))
			{
				double energy =
					CalculateParticleEnergy(temp, true) -
					CalculateParticleEnergy(particles.at(j), true);

				if (exp(energy) > GenerateRandom(0, 1, generator))
				{
					particles.at(j) = temp;
				}

				averageProjection += particles.at(j).mz;
			}
		}
	}

	cout << "100 %  Tube: " << tubeR << " " << tubeL << "   \r" << std::flush;
	cout << "\nIterations've finished!" << endl;
	cout << "Average is: " << averageProjection / stepsAmount / particleAmount << endl;
	
	cout 
		<< "Theory for isolated:  " 
		<< Theory::CalculateMagnetizationForSingleParticle(MCSim_application::GetField()) 
		<< endl;

	cout 
		<< "Theory for system:    " 
		<< Theory::CalculateMagnetizationForSystem(MCSim_application::GetField())
		<< endl;

	CheckSystemForErrors();

	cout << "****************ITERATIONS******************\n";
}


/* Make Resuzing function performs resizing of the tube to get particles'
   volume density be equal to target*/
void Simulator::MakeResizing()
{
	cout << "\n*****************RESIZING*******************\n";

	bool goOn = true;

	//TODO: finish the resizing procedure
	int k = 0;
	while (k < 2)
	{
		k++;

		tempR = tubeR;
		tempL = tubeL;

		cout
			<< " Tube: "
			<< tubeR << " "
			<< tubeL << " "
			<< " Max: "
			<< tempR << " "
			<< tempL << "          \n" << std::flush;

		for (int i = 0; i < particleAmount; i++)
		{
			Particle temp = GenerateDeltaState(particles.at(i));

			if (!CheckParticleForCollisions(temp, i))
			{
				double energy =
					CalculateParticleEnergy(temp, false) -
					CalculateParticleEnergy(particles.at(i), false);

				if (exp(energy) > GenerateRandom(0, 1, generator))
				{
					particles.at(i) = temp;
					CollectTubeSizes(i);
				}
			}

			ResizeTubeIfPossible();
			goOn = CheckIfTubeNeedResize();
		}

		cout 
			<< " Tube: " 
			<< tubeR << " " 
			<< tubeL << " " 
			<< " Max: " 
			<< tempR << " "
			<< tempL << "          \n\n" << std::flush;
	}

	cout << "Resize finished! Now will check if there's any errors..." << endl;
	CheckSystemForErrors();

	cout << "\n*****************RESIZING*******************\n";
}


/* Calculate Particle Energy function calculates energy of particles 
   inside the outer magnetic field and in dipole-dipole interaction */
double Simulator::CalculateParticleEnergy(const Particle particle, bool mode)
{
	//here's the field part of energy
	std::vector<double> point;
	point.push_back(particle.x);
	point.push_back(particle.y);
	point.push_back(particle.z);

	std::vector<double> outerField = GetVectorField(point, mode);

	std::vector<double> moment;
	moment.push_back(particle.mx);
	moment.push_back(particle.my);
	moment.push_back(particle.mz);

	double outerFieldEnergy = Particle::CalculateProjection(moment, outerField);
	//cout << "cos: " << cosTheta << endl;

	double dipoleInteractionEnergy = 0;


	//TODO: calculate dipole interaction energy
	//Already have dip-dip function in particle obj

	double totalEnergy = dipoleInteractionEnergy + outerFieldEnergy;

	return totalEnergy;
}

double Simulator::CalculateCurrentVolumeDensity()
{
	double particleVolume = particleAmount * 4 / 3 * PI * pow((particleDiameter / 2), 3);
	double tubeVolume = PI * tubeR * tubeR * tubeL;

	return particleVolume / tubeVolume;
}


/* Check For Collisions function checks if particle is collided with
   one of unignored particles or walls */
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


/* Generate Random function just generates random double value */
double Simulator::GenerateRandom(double min, double max, mt19937_64 &generator)
{
	uniform_real_distribution<double> distr(min, max);
	return distr(generator);
}


/* Set Random Seed function sets seed for random generator */
void Simulator::SetGeneratorRandomSeed()
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


/* Get Vector Field function returns outer field vector in two 
   modes: normal and radial */
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


/* Check System For Errors is backup function which checks system
   just one more time for collisions and maybe other errors */
bool Simulator::CheckSystemForErrors()
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


/* Collect Tube Sizes function checks if it's possible to remember
   new smallest tube's parameters */
void Simulator::CollectTubeSizes(int num)
{
	double R = sqrt(
		particles.at(num).x * particles.at(num).x +
		particles.at(num).y * particles.at(num).y);

	double L = abs(particles.at(num).z) * 2;

	if (R < this->tempR) this->tempR = R;
	if (L < this->tempL) this->tempL = L;
}


void Simulator::ResizeTubeIfPossible()
{
	if (tempL < tubeL) tubeL = tempL;
	if (tempR < tubeR) tubeR = tempR;
}


bool Simulator::CheckIfTubeNeedResize()
{
	double R = MCSim_application::GetTargetRadius();
	double A = MCSim_application::GetTubeAspect();

	if ((R > tubeR) && (R * A > tubeL))
	{
		return false;
	}
	else
	{
		return true;
	}
}
