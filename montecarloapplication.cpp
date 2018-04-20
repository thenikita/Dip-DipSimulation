#include <iostream>
#include <cmath>

#include "montecarloapplication.h"
#include "simulator.h"

using std::cin;
using std::cout;
using std::endl;

//
// IMPORTANT MEMO
// Do not forget that everything is normed to particles' 
// diameter which is 1, so all skalar operations
// should be calculated with this remark!
//

// TODO include particle diameter and mass to calculations

//
// IMPORTANT MEMO
// Tube's length is going to axis Z
//

MonteCarloApplication::MonteCarloApplication( ) {
    this->lambda = 1;
    this->fieldModule = 0;
    this->targetVolumeDensity = 0.4;
    this->particleCount = 1;
    this->targetTubeL = 0;
    this->targetTubeR = 0;
    this->aspect = 20;
    this->particleDiameter = 1;
    // TODO add particle mass
    this->particleMass = 1;
    this->particleMagneticMoment = 1;

    cout << "\n*************************************************************\n";
    cout << "You are welcome to Monte-Carlo Simulation...\n";
    cout << "Initialization of program..." << endl;
    cout << "\n*************************************************************\n";

    cout    << "Enter the amount of particles. "
            << "Default amount - "
            << particleCount
            << endl;

    unsigned int tempParticleCount;
    cin >> tempParticleCount;
    if ( tempParticleCount > particleCount )
        particleCount = tempParticleCount;
    cout    << "Particle Amount is set to "
            << particleCount
            << endl;

    cout    << "Enter volume density of particles in the tube, "
            << "should be less than "
            << targetVolumeDensity
            << endl;
    double tempVolDen;
    cin >> tempVolDen;
    if (( tempVolDen < targetVolumeDensity ) && ( tempVolDen != 0 ))
        targetVolumeDensity = tempVolDen;
    cout    << "Volume density is set to "
            << targetVolumeDensity
            << endl;
    }

void MonteCarloApplication::Run( ) {
    auto *simulator = new Simulator( &lambda,
                                     &targetVolumeDensity,
                                     &particleMagneticMoment,
                                     &particleCount,
                                     &particleDiameter,
                                     &aspect );

    // In future here can be more than one thread of simulations
    simulator->Run();

    //system( "pause" );
}

const double PI = 3.14159;