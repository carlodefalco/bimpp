

#include <stdio.h>
#include "lis_config.h"
#include "lis.h"

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <fstream>
#include <stdlib.h>
//#include <mpi.h> 
LIS_INT main(LIS_INT argc, char* argv[])
{
	int nnz,gn,is,ie,is_r,ie_r,n;
	int int_rank,rank,size;
	sparse_matrix       sp;
	std::vector<int> jc,ir;
	std::vector<double> xa;
		
	LIS_MATRIX A;
	LIS_INT *row,*col;
	LIS_SCALAR * value;
	lis_initialize(&argc, &argv);
	#ifdef USE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD,&int_rank);
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		rank=int_rank;
	#else
		rank=0;
		size=1;
	#endif
	
	if(rank==0)
		{
			std::cout << "\n\n*****\nLis Test 1\n*****\n";			
			sp.resize(6);
			sp[0][0]=10;
			sp[0][4]=-2;
			sp[1][0]=3;
			sp[1][1]=9;
			sp[1][5]=3;
			sp[2][1]=7;
			sp[2][2]=8;
			sp[2][3]=7;
			sp[3][0]=3;
			sp[3][2]=8;
			sp[3][3]=7;
			sp[3][4]=5;
			sp[4][1]=8;
			sp[4][3]=9;
			sp[4][4]=9;
			sp[4][5]=13;
			sp[5][1]=4;
			sp[5][4]=2;
			sp[5][5]=-1;
	
			sp.csr(xa,jc,ir,0);	
			
			std::cout<<"row: ";
			for(int i=0;i<ir.size();i++)
				std::cout<<ir[i]<<" ";
			std::cout<<std::endl;
			std::cout<<"col: ";
			for(int i=0;i<jc.size();i++)
				std::cout<<jc[i]<<" ";
			std::cout<<std::endl;
			std::cout<<"val: ";
			for(int i=0;i<xa.size();i++)
				std::cout<<xa[i]<<" ";
			std::cout<<std::endl;
			
			sp[0][0]=5;
			sp[4][5]=2;
			std::cout<<std::endl<<"Matrix Modified"<<std::endl;
			sp.csr_update(xa,jc,ir,0);	
	
			nnz=xa.size();
			gn=sp.size();
			
		}
	
  #ifdef USE_MPI
		MPI_Bcast(&nnz,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&gn,1,MPI_INT,0,MPI_COMM_WORLD);
	
		
		if(rank==0)
			{
				for(int k=1;k<size;++k)
					{
						for(int i=0;i<nnz;++i)
							{
								MPI_Send(&xa[i],1,MPI_DOUBLE,k,0,MPI_COMM_WORLD);
								MPI_Send(&jc[i],1,MPI_INT,k,0,MPI_COMM_WORLD);
							}
						for(int i=0;i<gn;++i)
							{						
								MPI_Send(&ir[i],1,MPI_INT,k,0,MPI_COMM_WORLD);							
							}
						MPI_Send(&ir[gn],1,MPI_INT,k,0,MPI_COMM_WORLD);
					}
			}
		else
			{
				xa.resize(nnz);
				jc.resize(nnz);
				ir.resize(gn+1);

				MPI_Status *status;
				
						for(int i=0;i<nnz;++i)
							{
								MPI_Recv(&xa[i],1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,status);
								MPI_Recv(&jc[i],1,MPI_INT,0,0,MPI_COMM_WORLD,status);
							}
			
						for(int i=0;i<gn;++i)
							{						
								MPI_Recv(&ir[i],1,MPI_INT,0,0,MPI_COMM_WORLD,status);							
							}
						
						MPI_Recv(&ir[gn],1,MPI_INT,0,0,MPI_COMM_WORLD,status);
					
			}
	#endif

	int counter=0;
	int counter_R=0;
	
	if(rank==0)
		{
			n=gn/size+gn%size;
			nnz=ir[n];
			is=0;
			ie=is+nnz;
			is_r=0;
			ie_r=n+1;
		}
	else
		{
			n=gn/size;
			nnz=ir[n*(rank+1)+gn%size]-ir[n*rank+gn%size];
			is=ir[n*rank+gn%size];
			ie=is+nnz;
			is_r=gn%size+n*rank;
			ie_r=is_r+n+1;
		}
	row=(LIS_INT *)malloc((n+1)*sizeof(LIS_INT));
	col=(LIS_INT *)malloc(nnz*sizeof(LIS_INT));
	value=(LIS_SCALAR *)malloc(nnz*sizeof(LIS_SCALAR));
	
	for(int i=is;i<ie;i++)
		{
			col[counter]=jc[i];
			value[counter]=xa[i];
			counter++;
		}
	for(int i=is_r;i<ie_r;++i)
		{
			if(rank!=0)
				row[counter_R]=ir[i]-ir[is_r];
			else
				row[counter_R]=ir[i];
			counter_R++;
		}
	
	lis_matrix_create(LIS_COMM_WORLD, &A);
	lis_matrix_set_size(A,n,0);
	lis_matrix_set_csr(nnz,row,col,value,A);
	
	LIS_INT error=lis_matrix_assemble(A);
	if(rank==0) std::cout<<"Matrix Assembled"<<std::endl;
	if(rank==0)
		{
			std::cout<<"row: ";
			for(int i=0;i<n+1;i++)
				std::cout<<row[i]<<" ";
			std::cout<<std::endl;
			std::cout<<"col: ";
			for(int i=0;i<nnz;i++)
				std::cout<<col[i]<<" ";
			std::cout<<std::endl;
			std::cout<<"val: ";
			for(int i=0;i<nnz;i++)
				std::cout<<value[i]<<" ";
			std::cout<<std::endl;
			
		}
	
	lis_matrix_destroy(A);
	lis_finalize();
  return 0;
}

