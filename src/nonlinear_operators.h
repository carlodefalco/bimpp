#ifndef HAVE_NL_OPERATORS_H
#define HAVE_NL_OPERATORS_H 1
#include <bim_sparse.h>
#include <operators.h>

struct Inexact_Newton_Option
{
	int maxIter; //Max number of Iteration

	double tol; //Iteration stops if ||x_new - x_old||<tol
	double minRes; //Iteration stops if ||F(x)||<minRes
	double forcing; //Linear Iteration stops if linear residual<forcing
	
	normType type; //Type of norm
	
};

struct Backtracking_Inexact_Newton_Option
{
	
	int maxIter; //Max number of Iteration

	double tol; //Iteration stops if ||x_new - x_old||<tol
	double minRes; //Iteration stops if ||F(x)||<minRes
	double forcing; //Linear Iteration stops if linear residual<forcing
	
	double t; 
	double theta_min; 
	double theta_max;
	double theta;

	normType type; //Type of norm
	
	///compute theta that minimizing over [theta_min, theta_max] the quadratic function f(theta)
	///for which f(0)=g(0), f'(0)=g'(0) and f(1)=g(1) where g(theta)=||F(x_k+theta*du_k||^2
	void
	theta_choice(double a, double b, double c)
	{
		double f_theta,f_min,f_max;
		std::cout<<"a "<<a<<std::endl;
		std::cout<<"b "<<b<<std::endl;
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
};


struct Inexact_Newton_Status
{
	int iteration; //Number of Iterations when Newton Method finished
	double residual; //Residual of function when Newton Method finished
	bool converged; //Flag is true when residual is less then minRes or the norm of difference between two iterations is less then tol
};

struct Stream_Option
{
	int verbosity; 
	char* filename;
};

///compute forcing term ||F(x_k)-F(x_k-1)-DF(x_k-1)du_k-1||/||F(x_k-1)||
struct ForcingType1
{
	double eta_max;

	template<typename P>
	double operator()(P& problem, std::vector<double>& uold, std::vector<double>& du, double eta_old, normType norm_type)
	{
		sparse_matrix df;
		std::vector<double> fnew,fold,unew,temp;
		double eta_new,eta_temp;
		unew.resize(du.size());
		for(int i=0;i<du.size();++i)
			unew[i]=uold[i]+du[i];

		problem(fnew,unew);
		problem(df,fold,uold);
		
		bim3a_matrix_vector_product(df,du,temp);
		
		for(int i=0;i<temp.size();++i)
			{
				temp[i]*=1;
				temp[i]+=fold[i]-fnew[i];
			}
		
		bim3a_norm(problem.msh,fnew,eta_new,norm_type);
		bim3a_norm(problem.msh,fold,eta_temp,norm_type);
		
		eta_new/=eta_temp;

		eta_temp=pow(eta_old,(1+sqrt(5))/2);
		if(eta_temp > 0.1)
			return std::min(std::max(eta_new, eta_temp),eta_max);
		else return std::min(eta_new, eta_max);
	}
};

///compute forcing term ||F(x_k)||-||F(x_k-1)+DF(x_k-1)du_k-1||/||F(x_k-1)||
struct ForcingType2
{
	double eta_max;

	template<typename P>
	double operator()(P& problem, std::vector<double>& uold, std::vector<double>& du, double eta_old, normType norm_type)
	{
		sparse_matrix df;
		std::vector<double> fnew,fold,unew,temp;
		double eta_new, eta_temp;
		unew.resize(du.size());
		for(int i=0;i<du.size();++i)
			unew[i]=uold[i]+du[i];

		problem(fnew,unew);
		problem(df,fold,uold);
		
		bim3a_matrix_vector_product(df,du,temp);
		
		for(int i=0;i<temp.size();++i)
			{
				temp[i]-=fold[i];
			}
		
		bim3a_norm(problem.msh,fnew,eta_new,norm_type);
		bim3a_norm(problem.msh,temp,eta_temp,norm_type);
		
		eta_new-=eta_temp;
	
		bim3a_norm(problem.msh,fold,eta_temp,norm_type);
	
		eta_new/=eta_temp;

		eta_temp=pow(eta_old,(1+sqrt(5))/2);
		if(eta_temp > 0.1)
			return std::min(std::max(eta_new, eta_temp),eta_max);
		else return std::min(eta_new, eta_max);
	}
};

///compute forcing term gamma*(||F(x_k)||/||F(x_k-1)||)^alpha
struct ForcingType3
{
	double gamma;
	double alpha;
	double eta_max;
	
	template<typename P>
	double operator()(P& problem, std::vector<double>& uold, std::vector<double>& du, double eta_old, normType norm_type)
	{
		std::vector<double> fnew,fold,unew,temp;
		double eta_new, eta_temp;
		unew.resize(du.size());
		for(int i=0;i<du.size();++i)
			unew[i]=uold[i]+du[i];

		problem(fnew,unew);
		problem(fold,uold);
		bim3a_norm(problem.msh,fnew,eta_new,norm_type);
		bim3a_norm(problem.msh,fold,eta_temp,norm_type);
		
		eta_new/=eta_temp;
		eta_new=gamma*pow(eta_new,alpha);
		
		eta_temp=gamma*pow(eta_old,alpha);

		if(eta_temp > 0.1)
			return std::min(std::max(eta_new, eta_temp),eta_max);
		else return std::min(eta_new, eta_max);
	}
};

///compute costant forcing temp
struct ForcingCostant
{
	template<typename P>
	double operator()(P &problem, std::vector<double>& uold, std::vector<double>& du, double eta_old, normType norm_type)
	{
		return eta_old;
	}
};
#endif
