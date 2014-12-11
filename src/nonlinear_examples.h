#ifndef HAVE_NL_EXAMPLES_H
#define HAVE_NL_EXAMPLES_H 1
#include <bim_sparse.h>
#include <operators.h>

class plaplacian{
	private:
		double p;
		std::vector<double> exactsol; //exact solution
		std::vector<double> f; //right hand side of p-laplacian
		std::vector<double> bsol; //Diriclhet boundary condition
		std::vector<int> bnodes; //Diriclhet boundary nodes
	public:
		mesh msh;

		plaplacian (){};
		~plaplacian(){};

		void read_mesh(const std::string mesh_name)
			{
				msh.read(mesh_name);
				msh.precompute_properties ();
			}

		void import(const double pp,
							 	 const std::vector<double>& esol, 
							 	 const std::vector<double>& ff, 
							 	 const std::vector<double>& bc,
							 	 const std::vector<int>& bn)
			{
				p=pp;
				exactsol=esol;
				f=ff;
				bsol=bc;
				bnodes=bn;
			}
		//compute left hand side and right and side of p-laplacian valued in u
		void operator()(sparse_matrix& lhs, std::vector<double>& rhs,const std::vector<double>& u);
		
		//compute functional in u
		void operator()(std::vector<double>& f,const std::vector<double>& u);

		//get the exact solution of the problem
		void get_solution(std::vector<double>& sol);
};

class equation{
	private:
		std::vector<double> exactsol;
	public:
		mesh msh;		
	
		equation(){};
		~equation(){};

		void import(const std::vector<double> esol)
			{
				exactsol=esol;
			}
		
		void operator()(sparse_matrix& lhs, std::vector<double>& rhs, std::vector<double>& u);
		void operator()(std::vector<double>& f, std::vector<double>& u);
		void get_solution(std::vector<double>& sol);
};
#endif
