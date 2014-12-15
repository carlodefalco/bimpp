#include <nonlinear_operators.h>

void
backtracking_inexact_newton_option::theta_choice(double a, double b, double c)
{
	double f_theta,f_min,f_max;
	theta=-b/(2*a);
	f_theta=a*theta*theta+b*theta+c;
	f_min=a*theta_min*theta_min+b*theta_min+c;
	f_max=a*theta_min*theta_min+b*theta_min+c;
	if(theta > theta_max || theta < theta_min)
		if(f_min < f_max)
			theta=theta_min;
		else
			theta=theta_max;		
	else
		if(f_theta <= f_min )
			if(f_theta >= f_max)
				theta=theta_max;
		else
			if(f_min <= f_max)
				theta=theta_min;
			else
				theta=theta_max;
}
