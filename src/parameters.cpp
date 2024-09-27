#include "parameters.h"

namespace param
{
	const double dt = 1e-6;
	const double a = 0.2;
	const double A_0 = 1.0;
	const double K_a = 10.0;
	
	const double l_min = 0.005*std::sqrt(A_0);
	const double l_new = 0.01*std::sqrt(A_0);
	const double A_min = 0.1*A_0;
	const double A_max = 2.0*A_0;
	
	double LAMBDA = 0;
	double GAMMA = 0.4*K_a*A_0;
	
	void set_LAMBDA(double LAMBDA_) { LAMBDA = LAMBDA_*K_a*std::powf(A_0, 1.5); }
	void set_GAMMA(double GAMMA_) { GAMMA = GAMMA_*K_a*A_0; }
}
