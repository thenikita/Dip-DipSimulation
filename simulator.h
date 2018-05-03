#pragma once

#include <vector>
#include <random>
#include <ctime>
#include <chrono>
#include <string>

#include "particle.h"

class Simulator {

public:
    const std::string CONSTANT = "const";
    const std::string LINEAR = "linear";
    const std::string CENTRAL = "central";

private:
    std::vector<Particle> particles;
    std::vector<double> field;

    double multiplier = 2;
    std::mt19937_64 generator;
    double deltaCoordinate = 0.1;

    unsigned int stepsAmount = 10000;

    double particleDiameter;
    unsigned int particleCount;
    unsigned int currentParticles;
    double lambda;
    double fieldModule;
    double targetVolumeDensity;
    double particleMagneticMoment = 1;
    double aspect;

    double tubeRadius, tubeLength;
    double tempRadius, tempLength;

    time_t timeAtStart;

    double singleTheoryMagnetization;
    double systemTheoryMagnetization;
    double currentMagnetization;

public:
    Simulator(
            double *lambda,
            double *targetVolumeDensity,
            double *particleMagneticMoment,
            unsigned int *particleCount,
            double *particleDiameter,
            double *aspect );

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
            double &tubeRadius,
            double &tubeLength,
            double aspect );

    void GenerateParticles( int particleCount );

    Particle GenerateDeltaState( Particle particle,
                                 double minX = -1.0,
                                 double maxX = 1.0,
                                 double minY = -1.0,
                                 double maxY = 1.0,
                                 double minZ = -1.0,
                                 double maxZ = 1.0 );

    double CalculateParticleEnergy( Particle particle,
                                    std::string mode = "const" );

    double CalculateCurrentVolumeDensity( );

    bool CheckParticleForCollisions( Particle particle, int ignored );

    static double GenerateRandom(
            double min,
            double max,
            std::mt19937_64 &generator );

    void SetGeneratorRandomSeed( );

    std::vector<double> GetVectorField(
            std::vector<double> point,
            std::string mode = "const" );

    bool CheckSystemForErrors( );

    void CollectTubeSizes( unsigned num );

    bool ResizeTubeIfPossible( bool resizeR,
                               bool resizeL);

    bool CheckIfTubeNeedResize( );

    void SetStartingTubeSize( );

    void ShowTimeSpent( time_t currentTime );

    void SetStartingTime( );

    void HelpMakeResizing( );

};
