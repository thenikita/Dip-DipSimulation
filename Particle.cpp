#include "Particle.h"
#include "MCSim_application.h"

#include <cmath>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

Particle::Particle(
	double x,
	double y,
	double z,
	double mx,
	double my,
	double mz)
{
	this->x = x;
	this->y = y;
	this->z = z;

	this->mx = mx;
	this->my = my;
	this->mz = mz;

}

double Particle::CalculateDipoleEnergy(
	double lambda, 
	Particle first, 
	Particle second)
{
	double d = 1;

	double moduleR = sqrt(pow(first.x - second.x, 2)
					+ pow(first.y - second.y, 2)
					+ pow(first.z = second.z, 2));

	vector<double> r, e1, e2;
	r.push_back(first.x - second.x);
	r.push_back(first.y - second.y);
	r.push_back(first.z - second.z);

	e1.push_back(first.mx);
	e1.push_back(first.my);
	e1.push_back(first.mz);

	e2.push_back(second.mx);
	e2.push_back(second.mx);
	e2.push_back(second.mx);

	double e1r = ScalarProduction(e1, r);
	double e2r = ScalarProduction(e2, r);
	double e1e2 = ScalarProduction(e1, e2);

	return -lambda / (pow(moduleR / d, 3)) *
			(3 * e1r * e2r / pow(moduleR / d, 2) - e1e2);
}

double Particle::CalculateInFieldEnergy(double field)
{
	double cos = 0;
	double lambda = MCSim_application::GetLambda();
	double xi = MCSim_application::GetField();


	return - field * cos;
}

/*Calculate Projection function returns projection of the second
  vector on the first vector's direction*/
double Particle::CalculateProjection(std::vector<double> first, std::vector<double> second)
{
	double scalar = ScalarProduction(first, second);
	double firstModule = calculateVectorModule(first);
	double secondmodule = calculateVectorModule(second);

	//cout << scalar << " " << firstModule << " " << secondmodule << endl;
	return scalar / firstModule;
}

double Particle::calculateVectorModule(std::vector<double> first)
{
	double sum = 0;

	for (int i = 0; i < 3; i++)
	{
		sum += pow(first.at(i), 2);
		//TODO: check if simpli multiplying will be faster than pow func
	}
	return sqrt(sum);
}

bool Particle::CheckForErrors(double R, double L)
{
	double x = this->x;
	double y = this->y;
	double z = this->z;

	double r =
		x * x +
		y * y;
	if (r > R) return true;
	if (z > L / 2) return true;

	double mx = this->mx;
	double my = this->my;
	double mz = this->mz;

	double moment = 
			mx * mx +
			my * my +
			mz * mz;
	if (moment > 1) return true;

	return false;
}

double Particle::ScalarProduction(vector<double> first, vector<double> second) {

	//cout << "scalar: " << first.at(0) << first.at(1) << first.at(2) << endl;
	return first.at(0) * second.at(0)
			+ first.at(1) * second.at(1)
			+ first.at(2) * second.at(2);
}

std::string Particle::toString()
{
	return "[" + std::to_string(this->x)
		+ ", " + std::to_string(this->y)
		+ ", " + std::to_string(this->z) + "]"
		+ "   ["
			   + std::to_string(this->mx)
		+ ", " + std::to_string(this->my)
		+ ", " + std::to_string(this->mz) + "]";
}