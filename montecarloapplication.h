#pragma once

class MonteCarloApplication {
private:
    double lambda;
    double fieldModule;                //xi module. vector is in simulator obj
    double targetVolumeDensity;  //phi
    double particleMagneticMoment;

    unsigned int particleCount;
    double particleDiameter;
    double particleMass;
    double aspect;
    double targetTubeR, targetTubeL;

public:
    MonteCarloApplication( );
    void Run( );
};
