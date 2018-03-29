#include <iostream>
#include <cmath>

#include "montecarloapplication.h"

using std::cin;
using std::cout;
using std::endl;

/*
double MonteCarloApplication::lambda = 1;
double MonteCarloApplication::fieldModule = 0;
double MonteCarloApplication::targetVolumeDensity = 0.4;
int MonteCarloApplication::particleCount = 1;
double MonteCarloApplication::targetTubeL = 0;
double MonteCarloApplication::targetTubeR = 0;
double MonteCarloApplication::aspect = 20;
double MonteCarloApplication::particleDiameter = 1;
double MonteCarloApplication::particleMagneticMoment = 1;
*/

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

MonteCarloApplication::MonteCarloApplication( ) {
    this->lambda = 1;
    this->fieldModule = 0;
    this->targetVolumeDensity = 0.4;
    this->particleCount = 1;
    this->targetTubeL = 0;
    this->targetTubeR = 0;
    this->aspect = 20;
    this->particleDiameter = 1;
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

    GenerateTube(particleCount,
                 particleDiameter,
                 targetVolumeDensity,
                 targetTubeR,
                 targetTubeL,
                 aspect);

    auto *simulator = new Simulator( &lambda,
                                     &targetVolumeDensity,
                                     &particleMagneticMoment,
                                     &particleCount,
                                     &particleDiameter,
                                     &aspect,
                                     targetTubeR,
                                     targetTubeL );

    simulator->Run();

    system( "pause" );
}

const double PI = 3.14159;

void MonteCarloApplication::Show( Simulator simulator ) {
    cout << "\n######################### SYSTEM ###########################\n";
    cout << "Lambda:           " << lambda << endl;
    cout << "Field Module:     " << fieldModule << endl;
    cout << "Density:          " << targetVolumeDensity << endl;

    cout << "Particles:        " << particleCount << endl;
    cout << "  N  |"
            "            COORDINATES             |"
            "            MAGNETIC MOMENT       "
         << endl;

    for ( int i = 0; i < particleCount; i++ ) {
        cout << "  " << i << "    ";
        cout << simulator.getParticles( ).at( i ).Show( ) << endl;
    }
    cout << "\n######################### SYSTEM ###########################\n";
}

void MonteCarloApplication::GenerateTube(
        unsigned int particleCount,
        double particleDiameter,
        double targetVolumeDensity,
        double &targetTubeR,
        double &targetTubeL,
        double aspect ) {

    cout << "\n*************************************************************\n";
    cout << "Now going to generate tube..." << endl;

    // N * 4/3 pi r^3
    double allParticlesVolume =
            particleCount * 4.0 / 3.0 * PI * pow( particleDiameter, 3 ) / 8.0;
    double tubeVolume = allParticlesVolume / targetVolumeDensity;

    // V = pi r^2 * L = pi r^3 * aspect
    targetTubeR = pow( tubeVolume / PI / aspect, 0.3333333 );

    if ( targetTubeR < particleDiameter ) {
        targetTubeR = 0.55 * particleDiameter;
    }

    targetTubeL = targetTubeR * aspect;

    cout    << "New generated tube: "
            << "R = " << targetTubeR
            << " and L = " << targetTubeL
            << endl;
}