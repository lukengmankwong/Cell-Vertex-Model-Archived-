#ifndef PARAMETERS_H
#define	PARAMETERS_H

namespace param
{
	constexpr double HALF_INV_PI = 0.15915494309189535;
	
	constexpr int CELL_COUNT = 10000;
	constexpr double dt = 1e-6;
	
	constexpr double A_0 = 1.0/CELL_COUNT;
	constexpr double K_a = 100.0;
	
	constexpr double LAMBDA_ = -0.85;
	constexpr double GAMMA_  = 0.1;
	constexpr double LAMBDA = LAMBDA_*K_a*std::powf(A_0, 1.5);
	constexpr double GAMMA 	= GAMMA_*K_a*A_0;
	
	
	constexpr double a = 25.0;
	
}

#endif // PARAMETERS_H
