#include "simulator.h"

#pragma once

class MonteCarloApplication {

private:
    double lambda;
    double fieldModule;                //xi module. vector is in simulator obj
    double targetVolumeDensity;  //phi
    double particleMagneticMoment;

    unsigned int particleCount;
    double particleDiameter;
    double aspect;
    double targetTubeR, targetTubeL;

public:
    MonteCarloApplication( );
    void Run( );

private:
    static void GenerateTube( unsigned int particleCount,
                              double particleDiameter,
                              double targetVolumeDensity,
                              double &targetTubeR,
                              double &targetTubeL,
                              double aspect );

    void Show( Simulator simulator );
};
