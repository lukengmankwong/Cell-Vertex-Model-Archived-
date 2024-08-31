#ifndef PARAMETERS_H
#define	PARAMETERS_H

namespace param
{

	constexpr double dt = 1e-6;
	
	constexpr double A_0 = 1.0;
	constexpr double K_a = 110.0;
	
	constexpr double LAMBDA_ = 0.12;
	constexpr double GAMMA_  = 0.04;
	constexpr double LAMBDA = LAMBDA_*K_a*std::powf(A_0, 1.5);
	constexpr double GAMMA 	= GAMMA_*K_a*A_0;

	constexpr double a = 0.2;
	
	constexpr double l_min = 0.05*std::sqrt(A_0);
	constexpr double A_min = 0.2*A_0;
	constexpr double A_max = 2.0*A_0;
	
}

#endif // PARAMETERS_H
