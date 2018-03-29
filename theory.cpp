#include "theory.h"
#include "montecarloapplication.h"
#include <math.h>

double Theory::CalculateMagnetizationForSingleParticle( double field ) {
    return 1 / tanh( field ) - 1 / field;
}

double Theory::CalculateFieldMultiplier( ) {
    double l = MonteCarloApplication::GetLambda( );
    double p = MonteCarloApplication::GetTargetVolumeDensity( );
    double m = MonteCarloApplication::GetParticleMagneticMoment( );

    return 4 * l * p / m;
}

double Theory::CalculateMagnetizationForSystem( double field ) {
    double m = MonteCarloApplication::GetParticleMagneticMoment( );
    double n = MonteCarloApplication::GetTargetVolumeDensity( );
    double M =
            CalculateMagnetizationForSingleParticle( field )
            * m
            * n;

    double H = field / CalculateFieldMultiplier( );
    double effectiveM = H + M / 3;

    return m * n * CalculateMagnetizationForSingleParticle( effectiveM );
}
