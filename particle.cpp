#include "particle.h"
#include "montecarloapplication.h"

#include <iostream>
#include <cmath>

using std::vector;
using std::cout;
using std::endl;

double Particle::d = 0.0;

Particle::Particle(
        double x,
        double y,
        double z,
        double mx,
        double my,
        double mz ) {
    this->x = x;
    this->y = y;
    this->z = z;

    this->mx = mx;
    this->my = my;
    this->mz = mz;
}

double Particle::CalculateDipoleInteractionEnergy(
        double lambda,
        Particle first,
        Particle second ) {
    double d = 1;

    double moduleR = sqrt( pow( first.x - second.x, 2 )
                           + pow( first.y - second.y, 2 )
                           + pow( first.z = second.z, 2 ));

    vector<double> r, e1, e2;
    r.push_back( first.x - second.x );
    r.push_back( first.y - second.y );
    r.push_back( first.z - second.z );

    e1.push_back( first.mx );
    e1.push_back( first.my );
    e1.push_back( first.mz );

    e2.push_back( second.mx );
    e2.push_back( second.mx );
    e2.push_back( second.mx );

    double e1r = ProductScalars( e1, r );
    double e2r = ProductScalars( e2, r );
    double e1e2 = ProductScalars( e1, e2 );

    return -lambda / ( pow( moduleR / d, 3 )) *
           ( 3 * e1r * e2r / pow( moduleR / d, 2 ) - e1e2 );
}

// TODO implement calculating energy in the outer field
double Particle::CalculateInFieldEnergy( double field ) {
    double cos = 0;
    //double lambda = MonteCarloApplication::GetLambda( );
    //double xi = MonteCarloApplication::GetField( );

    return -field * cos;
}

/*Calculate Projection function returns projection of the second
  vector on the first vector's direction*/
double Particle::CalculateProjection(
        std::vector<double> first,
        std::vector<double> second ) {

    double scalar = ProductScalars( first, second );
    double firstModule = CalculateVectorModule( first );
    // double secondmodule = CalculateVectorModule( second );

    //cout << scalar << " " << firstModule << " " << secondmodule << endl;
    return scalar / firstModule;
}

double Particle::CalculateVectorModule( std::vector<double> vector ) {
    double sum = 0;

    for ( unsigned i = 0; i < 3; i++ ) {
        sum += vector.at( i ) * vector.at( i );
    }
    return sqrt( sum );
}

bool Particle::CheckForErrors( double R, double L ) {
    double x = this->x;
    double y = this->y;
    double z = this->z;

    double r =
            x * x +
            y * y;
    if ( r > R ) return true;
    if ( z > L / 2 ) return true;

    double mx = this->mx;
    double my = this->my;
    double mz = this->mz;

    double moment =
            mx * mx +
            my * my +
            mz * mz;

    return moment > 1;
}

double Particle::ProductScalars( vector<double> first, vector<double> second ) {

    //cout << "scalar: " << first.at(0) << first.at(1) << first.at(2) << endl;
    return first.at( 0 ) * second.at( 0 )
           + first.at( 1 ) * second.at( 1 )
           + first.at( 2 ) * second.at( 2 );
}

std::string Particle::ToString( ) {
    return "[" + std::to_string( this->x )
           + ", " + std::to_string( this->y )
           + ", " + std::to_string( this->z ) + "]"
           + "   ["
           + std::to_string( this->mx )
           + ", " + std::to_string( this->my )
           + ", " + std::to_string( this->mz ) + "]";
}