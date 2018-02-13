#include "stdafx.h"
#include <iostream>
#include <math.h>

#include "MCSim_application.h"
#include "Simulator.h"

using std::cin;
using std::cout;
using std::endl;

double MCSim_application::lambda = 1;
double MCSim_application::field = 0;
double MCSim_application::targetVolumeDensity = 0.4;
int MCSim_application::particleAmount = 1;
double MCSim_application::targetTubeL = 0;
double MCSim_application::targetTubeR = 0;
double MCSim_application::aspect = 20;
double MCSim_application::particleDiameter = 1;

//
// IMPORTANT MEMO
// Do not forget that everything is normed to particles' 
// diameter which is 1, so all skalar operations
// should be calculated wuth this remark!
//

//
// IMPORTANT MEMO
// Tube's length is going to axis Z
//

MCSim_application::MCSim_application()
{
	//TODO:set default values for lambda xi etc
	cout << "You are wellcome to Monte-Carlo Simulation...\n";
	cout << "Initialization of program..." << endl;

	cout << "Enter the amount of particles. It should be larger that default amount: " << particleAmount << endl;
	int tempParticleAmount;
	cin >> tempParticleAmount;
	if (tempParticleAmount > particleAmount) particleAmount = tempParticleAmount;
	cout << "Particle Amount is set to: " << particleAmount << endl;

	cout << "Enter volume density of particles in the tube, should be less than " << targetVolumeDensity << endl;
	double tempVolDen;
	cin >> tempVolDen;
	if (tempVolDen < targetVolumeDensity) targetVolumeDensity = tempVolDen;
	cout << "Volume density is set to " << targetVolumeDensity << endl;

	cout << "Now going to generate tube..." << endl;
	GenerateTube();
	cout << "New generated tube is R = " << targetTubeR << " and L = " << targetTubeL << endl;
	cout << "Now going to generate particles..." << endl;

	Simulator simulator(targetTubeR, targetTubeL, particleAmount);

	showSystem(simulator);

	simulator.MakeIterations(particleAmount);

	showSystem(simulator);

	system("pause");
}

const double PI = 3.14159;

void MCSim_application::showSystem(Simulator simulator)
{
	cout << "Current state of system: " << endl;
	cout << "Lambda:    " << lambda << endl;
	cout << "Field:     " << field << endl;
	cout << "Density:   " << targetVolumeDensity << endl;

	cout << "Particles: " << particleAmount << endl;
	cout << "  N  |          COORDINATES           |          MAGNETIC MOMENT       " << endl;
	for (int i = 0; i < particleAmount; i++)
	{
		cout << "  "  << i << "    ";
		cout << simulator.getParticles().at(i).toString() << endl;
	}
}

void MCSim_application::GenerateTube()
{
	// N * 4/3 pi r^3
	double allParticlesVolume = particleAmount * 4.0 / 3.0 * PI * pow(particleDiameter, 3) / 8.0;
	double tubeVolume = allParticlesVolume / targetVolumeDensity;
	// V = pi r^2 * L = pi r^3 * aspect
	targetTubeR = pow(tubeVolume / PI / aspect, 0.3333333);
	targetTubeL = targetTubeR * aspect;
}