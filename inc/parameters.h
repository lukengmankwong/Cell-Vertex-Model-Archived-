#ifndef PARAMETERS_H
#define	PARAMETERS_H

#include <cmath>

namespace param
{
	extern const double dt;
	extern const double a;
	extern const double A_0;
	extern const double K_a;
	
	extern const double l_min;
	extern const double l_new;
	extern const double A_min;
	extern const double A_max;
	
	extern double LAMBDA;
	extern double GAMMA;
	
	void set_LAMBDA(double LAMBDA_);
	void set_GAMMA(double GAMMA_);
}

#endif // PARAMETERS_H
