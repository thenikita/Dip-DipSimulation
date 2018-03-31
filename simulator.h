#pragma once

#include <vector>
#include <random>
#include <ctime>
#include <chrono>

#include "particle.h"

class Simulator {

public:

private:
    std::vector<Particle> particles;
    std::vector<double> field;

    double multiplier = 100;
    std::mt19937_64 generator;
    double deltaCoordinate = 0.1;

    unsigned int stepsAmount = 10;

    double particleDiameter;
    unsigned int particleCount;
    unsigned int currentParticles;
    double lambda;
    double fieldModule;
    double targetVolumeDensity;
    double particleMagneticMoment;
    double aspect;

    double tubeRadius, tubeLength;
    double tempRadius, tempLength;

public:
    Simulator(
            double *lambda,
            double *targetVolumeDensity,
            double *particleMagneticMoment,
            unsigned int *particleCount,
            double *particleDiameter,
            double *aspect,
            double &targetTubeR,
            double &targetTubeL );

    void Run();

    std::vector<Particle> GetParticles( ) { return particles; }

    void MakeIterations( int particleAmount);

    void MakeResizing( );

    double GetTubeR( ) { return tubeRadius; }

    double GetTubeL( ) { return tubeLength; }

    double SetTubeR( double R ) { this->tubeRadius = R; }

    double SetTubeL( double L ) { this->tubeLength = L; }

    double getLambda( ) const { return lambda; }

    double getFieldModule( ) const { return fieldModule; }

    void ShowSystem( );

private:
    static void GenerateTube(
            unsigned int particleCount,
            double particleDiameter,
            double targetVolumeDensity,
            double &targetTubeR,
            double &targetTubeL,
            double aspect );

    void GenerateParticles( int particleCount );

    Particle GenerateDeltaState( Particle particle, bool ifNeedResize );

    double CalculateParticleEnergy( Particle particle, bool mode );

    double CalculateCurrentVolumeDensity( );

    bool CheckParticleForCollisions( Particle particle, int ignored );

    static double GenerateRandom(
            double min,
            double max,
            std::mt19937_64 &generator );

    void SetGeneratorRandomSeed( );

    std::vector<double> GetVectorField(
            std::vector<double> point,
            bool mode );

    bool CheckSystemForErrors( );

    void CollectTubeSizes( unsigned num );

    void ResizeTubeIfPossible( );

    bool CheckIfTubeNeedResize( );

};
