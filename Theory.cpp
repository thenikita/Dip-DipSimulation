#include "stdafx.h"
#include "Theory.h"
#include <math.h>

double Theory::CalculateForSingleParticle(double field)
{
	return 1 / tanh(field) - 1 / field;
}
