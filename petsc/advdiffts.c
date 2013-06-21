/*
 *  advdiffts.c
 *  
 *
 *  Created by Carlo de Falco on 23/10/05.
 *
 *   Copyright (C) 2005  Carlo de Falco ( carlo@mate.polimi.it )
 * this program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * this program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MNME2D; if not, write to the Free Software
 * Foundation, Inc.,\n 59 Temple Place, Suite 330, Boston, MA  02111-1307
 *         USA
 *
 */

/* 
   Program usage:  mpirun -np <procs>  advdiffts -prefix "prefix"\
   -deltat dt -maxts mts -epsilon e -ts_type TS_TYPE [-help] \
   [all PETSc options] 
*/

static char help[] = " \
experiment with parallel unstructured meshes:\n \
Solve a linear Advection-Diffusion problem with \
exponentially fitted finite elements\n";

#include "petscao.h" 
#include "petscts.h" 
#include "petscmath.h"

#define __pre__

#undef __FUNCT__
#define __FUNCT__ "main"

#define MAX_IF_NODES 1000

/******************************************************
Application Data Structures
*******************************************************/
typedef struct {		
  PetscInt  vertex[3],elnum;		
  PetscReal area,shgx[3],shgy[3];
  PetscInt  bccode[3];
} Ttriangle;

typedef struct{
  Ttriangle               *triangles;
  AO                      reordering;
  ISLocalToGlobalMapping  mapping;
  VecScatter              scatter;
  PetscReal               *locx,*locy;
  Mat                     P;
  int                     *egen,toten,*egnn;
  int                     *startnn,totnn,totlocnn;
}Tmesh;

typedef struct{
  PetscReal               *source;
  PetscReal               *potential;
  PetscReal               *density;
}Tdata;

typedef struct{
  Tmesh                   *mesh;
  char                    pref[50];
  KSP                     ksp;
}Tmonitor;
/*******************************************************/


/******************************************************
Application Routines
*******************************************************/
PetscErrorCode permutationmatrix(MPI_Comm com,int egnn, int totnn, int* petscordering, 
				 int* appordering, Mat* M );
PetscErrorCode readmesh         (Tmesh *mesh,Tdata *data, Vec *D,const char prefix[] );
PetscErrorCode destroymesh      (Tmesh *mesh);
PetscErrorCode destroydata      (Tdata *data);
PetscErrorCode compshg          (PetscReal x0,PetscReal y0,PetscReal x1,PetscReal y1,
				 PetscReal x2,PetscReal y2, PetscReal *shgx, PetscReal *shgy);
PetscErrorCode compsg           (Tmesh *mesh, Tdata *data, Mat *SG, Mat *M, PetscReal epsilon);
PetscReal      bern             (PetscReal x);
PetscErrorCode monitor          (TS ts,PetscInt steps,PetscReal time,Vec x,void *mctx);
/*******************************************************/

/******************************************************
Global vriables for profiling Info
******************************************************/
PetscEvent             USER_EVENT;
/*******************************************************/

int main(int argc,char **argv)
{
  PetscErrorCode         ierr;
  Tmesh                  mesh;
  Tdata                  data;
  Mat                    SG,M;
  Vec                    sol,Msol;
  PetscMPIInt            rank,size;
  char                   prefix[50],filename[50];
  PetscTruth             flg;
  TS                     ts;
  PetscReal              zero=0.0,mone=-1.0;
  Tmonitor               ctxmonitor;
  PetscViewer            viewer;
   
  PetscInitialize(&argc,&argv,PETSC_NULL,help);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);CHKERRQ(ierr);
	
  ierr = PetscLogEventRegister(&USER_EVENT,"* Bernoulli *",PETSC_VIEWER_COOKIE);CHKERRQ(ierr);
  
  /*
    Read Command Line Options
  */

  PetscOptionsGetString(PETSC_NULL,"-prefix",prefix,49,&flg);
  if(!flg)
    sprintf(prefix,"nproc2/");

  
  PetscReal deltat;
  PetscOptionsGetReal(PETSC_NULL,"-deltat",&deltat,&flg);
  if(!flg)
    deltat = (PetscReal) 0.1; 
  
  PetscInt maxts;
  PetscOptionsGetInt(PETSC_NULL,"-maxts",&maxts,&flg);
  if(!flg)
    maxts = (PetscInt) 10;

  PetscReal epsilon;
  PetscOptionsGetReal(PETSC_NULL,"-epsilon",&epsilon,&flg);
  if(!flg)
    epsilon = (PetscReal) 0.000001;

  /*
    Read mesh
  */

  PreLoadBegin(PETSC_TRUE,"Read Data");
  readmesh(&mesh,&data,&sol,prefix);

  
  /*
    Build linear system matrices
  */

  PreLoadStage("Build System");
  ierr = PetscPrintf(MPI_COMM_WORLD,
		     "\n\n****\nBuild linear system matrices\n****\n");
  CHKERRQ(ierr);
  
  compsg(&mesh,&data,&SG,&M,epsilon);
  
  VecDuplicate(sol,&Msol);
  PetscObjectSetName((PetscObject) sol,"sol");
  PetscObjectSetName((PetscObject) Msol,"Msol");
  
  MatMult(M,sol,Msol);

  ierr = MatScale(SG,mone*epsilon);

  PetscObjectSetName((PetscObject) SG,"SG");
  PetscObjectSetName((PetscObject) mesh.P,"P");
  PetscObjectSetName((PetscObject) M,"M");

#undef __printmat__
#ifdef __printmat__
  PetscViewerASCIIOpen(MPI_COMM_WORLD,"matrici.m",&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  
  MatView(SG,viewer);
  MatView(mesh.P,viewer);
  MatView(M,viewer);
  VecView(sol,viewer);
  VecView(Msol,viewer);

  PetscViewerDestroy(viewer);
#endif

  /*
    Step time
  */	

  PreLoadStage("Step Time");
  
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nStep time \n****\n");
  CHKERRQ(ierr);
  
  ierr = TSCreate(MPI_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_LINEAR); CHKERRQ(ierr); 
  ierr = TSSetInitialTimeStep(ts,zero,deltat); CHKERRQ(ierr);
  ierr = TSSetDuration(ts,maxts,(PetscReal)maxts * deltat ); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,Msol); CHKERRQ(ierr);
  ierr = TSSetRHSMatrix(ts,SG,SG,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
  ierr = TSSetType(ts,TS_BEULER); CHKERRQ(ierr);


  ctxmonitor.mesh = &mesh;
  sprintf(ctxmonitor.pref,"%s",prefix);
  KSPCreate(MPI_COMM_WORLD,&(ctxmonitor.ksp));
  KSPSetOperators(ctxmonitor.ksp,M,M,SAME_NONZERO_PATTERN);
  KSPSetFromOptions(ctxmonitor.ksp);

  ierr = TSSetMonitor(ts,monitor,(void *)(&ctxmonitor),PETSC_NULL); 
  CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

  PetscInt steps; PetscReal ftime;
  ierr = TSStep(ts,&steps,&ftime); CHKERRQ(ierr);
  
  
  /*
    Clean up
  */
  
    PreLoadStage("Clean up");
    
    ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nClean up\n****\n");CHKERRQ(ierr);
    
    ierr=destroymesh(&mesh);CHKERRQ(ierr);
    ierr=destroydata(&data);CHKERRQ(ierr);
    
    ierr = 	MatDestroy(SG);CHKERRQ(ierr);
    ierr = 	MatDestroy(M);CHKERRQ(ierr);
    ierr = 	VecDestroy(sol);CHKERRQ(ierr);
    ierr = 	VecDestroy(Msol);CHKERRQ(ierr);
	ierr =  KSPDestroy(ksp);CHKERRQ(ierr);
    PreLoadEnd();
    
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
    
}

/************************************************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "readmesn"
/************************************************************************************************************************/
PetscErrorCode readmesh   (Tmesh *mesh,Tdata *data,Vec *D,const char prefix[] ){
	
	
  Vec                    x,y,V,S;
  Vec                    locvx,locvy,locvv,locvs,locvd;
  PetscMPIInt            rank,size;
  PetscErrorCode         ierr;
  PetscInt               locind[3],*bcnodes,nbc;
	
  char                   filename[50];
  double                 rx,ry,rv,rs,rd;
  PetscReal              *locx,*locy,*locv,*locs,*locd;
	
  int                    ii,jj,kk,ll[3],mm,*appordering,*petscordering;
	
	
  FILE                   *fptr1,*fptr2,*fptr3;
  IS                     from, to; 
	
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);CHKERRQ(ierr);
	
	
  ierr = PetscMalloc(size * sizeof(PetscInt),&mesh->egen); CHKERRQ(ierr);
  ierr = PetscMalloc(size * sizeof(PetscInt),&mesh->egnn); CHKERRQ(ierr);
  ierr = PetscMalloc(size * sizeof(PetscInt),&mesh->startnn); CHKERRQ(ierr);
	
  /***********************************/
  /*
    Read node partitioning info
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nRead node partitioning info\n****\n");CHKERRQ(ierr);
	
  sprintf(filename,"%st.npart.%d",prefix,size);
  fptr1 = fopen(filename,"r"); 
  if (!fptr1) 
    SETERRQ(0,"Could not open file");
	
  mm=0;
  mesh->totnn=0;
	
  /*count noodes*/
  PetscMemzero(mesh->egnn,size * sizeof(PetscInt));
  while(mm!=EOF) {
    mm=fscanf(fptr1,"%d",&kk);
    if(mm!=EOF){
      mesh->egnn[kk]++;
      mesh->totnn++;
    }
  }
	
  PetscMemzero(mesh->startnn,size * sizeof(PetscInt));
  for (ii=0;ii<size;ii++){
    if(ii)
      for(jj=0;jj<ii;jj++) 
	mesh->startnn[ii]+=mesh->egnn[jj];
  }
	
  PetscPrintf(MPI_COMM_WORLD,"The mesh has a total of %d nodes\n",mesh->totnn);
  for(ii=0;ii<size;ii++)
    PetscPrintf(MPI_COMM_WORLD,"Process %d owns %d nodes\n",ii,mesh->egnn[ii]);
	
  fclose(fptr1);
	
  fptr1 = fopen(filename,"r"); 
  if (!fptr1) 
    SETERRQ(0,"Could not open file");
	
  sprintf(filename,"%sp",prefix);
  fptr2 = fopen(filename,"r"); 
  if (!fptr2) 
    SETERRQ(0,"Could not open file");
	
  sprintf(filename,"%sdata",prefix);
  fptr3 = fopen(filename,"r"); 
  if (!fptr3) 
    SETERRQ(0,"Could not open file");
	
  ierr = PetscMalloc((mesh->egnn[rank] + MAX_IF_NODES) * sizeof(PetscInt),&appordering); CHKERRQ(ierr);
  ierr = PetscMalloc((mesh->egnn[rank] + MAX_IF_NODES) * sizeof(PetscInt),&petscordering); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->egnn[rank] * sizeof(PetscReal),&mesh->locx); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->egnn[rank] * sizeof(PetscReal),&mesh->locy); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->egnn[rank] * sizeof(PetscReal),&data->potential); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->egnn[rank] * sizeof(PetscReal),&data->density); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->egnn[rank] * sizeof(PetscReal),&data->source); CHKERRQ(ierr);
	
	
  VecCreateMPI(MPI_COMM_WORLD,mesh->egnn[rank],mesh->totnn,&x);
  VecDuplicate(x,&y);VecDuplicate(x,&V);VecDuplicate(x,&S);VecDuplicate(x,D);
	
  jj=0;
  for(ii=0;ii<mesh->totnn;ii++)
    {
		
      fscanf(fptr1,"%d",&kk);
      fscanf(fptr2,"%lg %lg",&rx,&ry);
      fscanf(fptr3,"%lg %lg %lg",&rv,&rs,&rd);
		
      if(kk==(int)rank){
			
	petscordering[jj]=jj+mesh->startnn[rank];
	appordering[jj]=ii;
			
	mesh->locx[jj]=(PetscReal)rx;
	mesh->locy[jj]=(PetscReal)ry;				
			
	data->potential[jj]=(PetscReal)rv;
	data->density[jj]=(PetscReal)rd;
	data->source[jj]=(PetscReal)rs;				
			
	jj++;
      }
    }
	
  fclose(fptr1);
  fclose(fptr2);		
  fclose(fptr3);		
	
  VecSetValues(x,mesh->egnn[rank],petscordering,mesh->locx,INSERT_VALUES);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
	
  VecSetValues(y,mesh->egnn[rank],petscordering,mesh->locy,INSERT_VALUES);		
  VecAssemblyBegin(y);
  VecAssemblyEnd(y);
	
  VecSetValues(V,mesh->egnn[rank],petscordering,data->potential,INSERT_VALUES);		
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
	
  VecSetValues(*D,mesh->egnn[rank],petscordering,data->density,INSERT_VALUES);		
  VecAssemblyBegin(*D);
  VecAssemblyEnd(*D);	
	
  VecSetValues(S,mesh->egnn[rank],petscordering,data->source,INSERT_VALUES);		
  VecAssemblyBegin(S);
  VecAssemblyEnd(S);
	
  ierr = PetscFree(mesh->locx); CHKERRQ(ierr);
  ierr = PetscFree(mesh->locy); CHKERRQ(ierr);
  ierr = PetscFree(data->potential); CHKERRQ(ierr);
  ierr = PetscFree(data->density); CHKERRQ(ierr);
  ierr = PetscFree(data->source); CHKERRQ(ierr);
	
  /***********************************/
  /*
    Build AO object
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nBuild AO object\n****\n");CHKERRQ(ierr);
  ierr = AOCreateBasic(MPI_COMM_WORLD,mesh->egnn[rank],appordering,petscordering,&mesh->reordering);  CHKERRQ(ierr);
	
  /***********************************/
  /*
    Build permutation matrix
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nBuild permutation matrix\n****\n");CHKERRQ(ierr);
  ierr = permutationmatrix(MPI_COMM_WORLD,mesh->egnn[rank],mesh->totnn,petscordering,appordering,&mesh->P);CHKERRQ(ierr);
	
  /***********************************/
  /*
    Read element partitioning info
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nRead element partitioning info\n****\n");CHKERRQ(ierr);
	
  sprintf(filename,"%st.epart.%d",prefix,size);
  fptr1 = fopen(filename,"r"); 
  if (!fptr1) 
    SETERRQ(0,"Could not open file");
	
	
  mm=0;
  mesh->toten=0;
	
  PetscMemzero(mesh->egen,size * sizeof(PetscInt));
	
  /*count elements*/
  while(mm!=EOF) {
    mm=fscanf(fptr1,"%d",&kk);
    if(mm!=EOF){
      mesh->egen[kk]++;
      mesh->toten++;
    }
  }
	
  fclose(fptr1);
	
  PetscPrintf(MPI_COMM_WORLD,"The mesh has a total of %d elements\n",mesh->toten);
  for(ii=0;ii<size;ii++)
    PetscPrintf(MPI_COMM_WORLD,"Process %d owns %d elements\n",ii,mesh->egen[ii]);
	
  fptr1 = fopen(filename,"r"); 
  if (!fptr1) 
    SETERRQ(0,"Could not open file");
	
  sprintf(filename,"%st",prefix);
  fptr2 = fopen(filename,"r"); 
  if (!fptr2) 
    SETERRQ(0,"Could not open file");
	
  fscanf(fptr2,"%*d %*d\n");
	
  ierr = PetscMalloc(mesh->egen[rank] * sizeof(Ttriangle),&mesh->triangles); CHKERRQ(ierr);
	
  jj=0;
  for(ii=0;ii<mesh->toten;ii++){
		
    fscanf(fptr1,"%d",&kk);
    fscanf(fptr2,"%d %d %d",&ll[0],&ll[1],&ll[2]);			
		
    if(kk==(int)rank){
      mesh->triangles[jj].elnum=ii;
      mesh->triangles[jj].vertex[0]=ll[0]-1;
      mesh->triangles[jj].vertex[1]=ll[1]-1;
      mesh->triangles[jj].vertex[2]=ll[2]-1;
      jj++;		
			
    }
  }
	
  fclose(fptr1);
  fclose(fptr2);
	
  /***********************************/
  /*
    Cycle through elements to build local2global mapping
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nCycle through elements to build local2global mapping\n****\n");CHKERRQ(ierr);
	
  mesh->totlocnn = mesh->egnn[rank];
	
  for(ii=0;ii<mesh->egen[rank];ii++){	
    for(jj=0;jj<3;jj++){
      for(kk=0;kk<mesh->totlocnn;kk++){
	if (mesh->triangles[ii].vertex[jj]==appordering[kk])
	  goto found;
      }
    found:
      if(kk==mesh->totlocnn){
	mesh->totlocnn++;
	appordering[kk]=mesh->triangles[ii].vertex[jj];
      }
			
    }
  }
	
	
  PetscSynchronizedFlush(MPI_COMM_WORLD);
	
  ierr = 	VecCreateSeq(PETSC_COMM_SELF,mesh->totlocnn,&locvx);  CHKERRQ(ierr);
  ierr = 	VecDuplicate(locvx,&locvy);  CHKERRQ(ierr);
  ierr = 	VecDuplicate(locvx,&locvv);  CHKERRQ(ierr);
  ierr = 	VecDuplicate(locvx,&locvs);  CHKERRQ(ierr);
  ierr = 	VecDuplicate(locvx,&locvd);  CHKERRQ(ierr);	
	
  ISLocalToGlobalMappingCreate(MPI_COMM_WORLD,mesh->totlocnn,appordering,&mesh->mapping);
	
  AOApplicationToPetsc(mesh->reordering,mesh->totlocnn,appordering);
	
  ierr = 	ISCreateGeneral(PETSC_COMM_SELF,mesh->totlocnn,appordering,&from);  CHKERRQ(ierr);
  ierr = 	ISCreateStride(PETSC_COMM_SELF,mesh->totlocnn,(PetscInt)0,(PetscInt)1,&to);  CHKERRQ(ierr);
	
  ierr = 	VecScatterCreate(x,from,locvx,to,&mesh->scatter);  CHKERRQ(ierr);
	
	
  ierr = 	VecScatterBegin(x,locvx,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
  ierr = 	VecScatterEnd(x,locvx,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
	
  ierr = 	VecScatterBegin(y,locvy,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
  ierr = 	VecScatterEnd(y,locvy,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
	
  ierr = 	VecScatterBegin(V,locvv,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
  ierr = 	VecScatterEnd(V,locvv,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
	
  ierr = 	VecScatterBegin(S,locvs,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
  ierr = 	VecScatterEnd(S,locvs,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
	
  ierr = 	VecScatterBegin(*D,locvd,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
  ierr = 	VecScatterEnd(*D,locvd,INSERT_VALUES,SCATTER_FORWARD,mesh->scatter);  CHKERRQ(ierr);
	
	
  ierr = PetscMalloc(mesh->totlocnn * sizeof(PetscReal),&mesh->locx); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->totlocnn * sizeof(PetscReal),&mesh->locy); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->totlocnn * sizeof(PetscReal),&data->potential); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->totlocnn * sizeof(PetscReal),&data->density); CHKERRQ(ierr);
  ierr = PetscMalloc(mesh->totlocnn * sizeof(PetscReal),&data->source); CHKERRQ(ierr);
	
  ierr = 	VecGetArray(locvx,&locx); CHKERRQ(ierr);
  ierr = 	VecGetArray(locvy,&locy); CHKERRQ(ierr);
  ierr = 	VecGetArray(locvv,&locv); CHKERRQ(ierr);
  ierr = 	VecGetArray(locvs,&locs); CHKERRQ(ierr);
  ierr = 	VecGetArray(locvd,&locd); CHKERRQ(ierr);
	
  PetscMemcpy(mesh->locx,locx,mesh->totlocnn * sizeof(PetscReal));
  PetscMemcpy(mesh->locy,locy,mesh->totlocnn * sizeof(PetscReal));
  PetscMemcpy(data->potential,locv,mesh->totlocnn * sizeof(PetscReal));
  PetscMemcpy(data->density,locd,mesh->totlocnn * sizeof(PetscReal));
  PetscMemcpy(data->source,locs,mesh->totlocnn * sizeof(PetscReal));
	
  ierr = 	VecRestoreArray(locvx,&locx); CHKERRQ(ierr);
  ierr = 	VecRestoreArray(locvy,&locy); CHKERRQ(ierr);
  ierr = 	VecRestoreArray(locvv,&locv); CHKERRQ(ierr);
  ierr = 	VecRestoreArray(locvs,&locs); CHKERRQ(ierr);
  ierr = 	VecRestoreArray(locvd,&locd); CHKERRQ(ierr);
	
  AOPetscToApplication(mesh->reordering,mesh->totlocnn,appordering);
	
  /***********************************/
  /*
    Read BC information
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nRead BC information\n****\n");CHKERRQ(ierr);
	
  sprintf(filename,"%sbc",prefix);
  fptr1 = fopen(filename,"r"); 
  if (!fptr1) 
    SETERRQ(0,"Could not open file");
	
  mm=0;
  nbc=0;
	
  /*count BC noodes*/
  while(mm!=EOF) {
    mm=fscanf(fptr1,"%d",&kk);
    if(mm!=EOF){
      nbc++;
    }
  }
	
	
  fclose(fptr1);
	
  ierr =	PetscPrintf(MPI_COMM_WORLD,"%d Boundary nodes\n",nbc);CHKERRQ(ierr);
	
  fptr1 = fopen(filename,"r"); 
  if (!fptr1) 
    SETERRQ(0,"Could not open file");
	
  ierr = PetscMalloc(nbc * sizeof(PetscInt),&bcnodes); CHKERRQ(ierr);
	
  for(ii=0;ii<nbc;ii++){
    fscanf(fptr1,"%d",&bcnodes[ii]);
    bcnodes[ii]=bcnodes[ii]-1;
  }
	
  fclose(fptr1);
	
  /***********************************/
  /*
    Cycle through elements and compute areas and shape function gradients
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,
		     "\n\n****\nCycle through elements and compute areas and shape function gradients\n****\n");
  CHKERRQ(ierr);
	
  
  for(ii=0;ii<mesh->egen[rank];ii++)
    {
      PetscMemzero(mesh->triangles[ii].bccode,3 * sizeof(PetscInt));
      for(mm=0;mm<3;mm++){
	jj=0;
	while((jj<nbc)&&(bcnodes[jj]<=mesh->triangles[ii].vertex[mm])){
	  if(bcnodes[jj]==mesh->triangles[ii].vertex[mm])
	    mesh->triangles[ii].bccode[mm]=1;
	  jj++;
	}
      }
      
      ISGlobalToLocalMappingApply(mesh->mapping,IS_GTOLM_DROP,
				  3,mesh->triangles[ii].vertex,&kk,locind);
		
      if (kk!=3){
	sprintf(filename,"i nodi %d %d %d hanno numerazione locale %d %d %d",
		mesh->triangles[ii].vertex[0],mesh->triangles[ii].vertex[1],mesh->triangles[ii].vertex[2],
		locind[0],locind[1],locind[2]);	
	SETERRQ(1,filename)
	  }
		
      PetscReal x0,y0,x1,y1,x2,y2;
		
      x0 = mesh->locx[locind[0]];
      y0 = mesh->locy[locind[0]];
      x1 = mesh->locx[locind[1]];
      y1 = mesh->locy[locind[1]];
      x2 = mesh->locx[locind[2]];
      y2 = mesh->locy[locind[2]];
		
      mesh->triangles[ii].area = PetscAbs((x1-x0)*(y2-y0)-(y1-y0)*(x2-x0))/2.0;
		
      compshg( x0,y0,x1,y1,x2,y2,mesh->triangles[ii].shgx,mesh->triangles[ii].shgy);
		
    }
	
  /***********************************/
  /*
    Write output
  */
#undef __DOIT__
#ifdef __DOIT__
	
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nWrite output\n****\n");CHKERRQ(ierr);
	
  sprintf(filename,"%sbcs",prefix);
  if(rank==0)
    fptr1 = fopen(filename,"w"); 
  else
    fptr1 = fopen(filename,"a"); 
	
	
  if (!fptr1) 
    SETERRQ(0,"Could not open file");
	
  for(ii=0;ii<mesh->egen[rank];ii++)
    {
      PetscSynchronizedFPrintf(MPI_COMM_WORLD,fptr1,"%d: (%d,%d),(%d,%d),(%d,%d)\n",mesh->triangles[ii].elnum,
			       mesh->triangles[ii].vertex[0],mesh->triangles[ii].bccode[0],
			       mesh->triangles[ii].vertex[1],mesh->triangles[ii].bccode[1],
			       mesh->triangles[ii].vertex[2],mesh->triangles[ii].bccode[2]);
    }
	
  PetscSynchronizedFlush(MPI_COMM_WORLD);
  fclose(fptr1);
#endif  //__DOIT__
	
	/***********************************/
	/*
	  Free Memory
	*/
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nFree memory\n****\n");CHKERRQ(ierr);
	
  ierr = PetscFree(appordering); CHKERRQ(ierr);
  ierr = PetscFree(petscordering); CHKERRQ(ierr);
	
  ierr = VecDestroy(locvx); CHKERRQ(ierr);
  ierr = VecDestroy(locvy); CHKERRQ(ierr);
	
  ierr = VecDestroy(x); CHKERRQ(ierr);
  ierr = VecDestroy(y); CHKERRQ(ierr);	
	
  ierr = VecDestroy(locvv); CHKERRQ(ierr);
  ierr = VecDestroy(locvs); CHKERRQ(ierr);
  ierr = VecDestroy(locvd); CHKERRQ(ierr);
	
  ierr = VecDestroy(V); CHKERRQ(ierr);
  ierr = VecDestroy(S); CHKERRQ(ierr);	
  /***********************************/
  /*
    End
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,
		     "\n\n****************************\nEnd of Data Input!!\n*******************************\n");
  CHKERRQ(ierr);
	
  return ierr;
	
}

/************************************************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "destroymesh"
/************************************************************************************************************************/

PetscErrorCode destroymesh(Tmesh *mesh ){
	
  PetscErrorCode         ierr;
	
  ierr=PetscFree(mesh->triangles);CHKERRQ(ierr);
  ierr=AODestroy(mesh->reordering);CHKERRQ(ierr);
  ierr=ISLocalToGlobalMappingDestroy(mesh->mapping);CHKERRQ(ierr);
  ierr=VecScatterDestroy(mesh->scatter);CHKERRQ(ierr);
  ierr=PetscFree(mesh->locx);CHKERRQ(ierr);
  ierr=PetscFree(mesh->locy);CHKERRQ(ierr);
  ierr=MatDestroy(mesh->P);CHKERRQ(ierr);
  ierr=PetscFree(mesh->egen);CHKERRQ(ierr);
  ierr=PetscFree(mesh->egnn);CHKERRQ(ierr);
  ierr=PetscFree(mesh->startnn);CHKERRQ(ierr);
	
  return ierr;
}

/************************************************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "destroydata"
/************************************************************************************************************************/

PetscErrorCode destroydata(Tdata *data ){
	
  PetscErrorCode         ierr;
	
  ierr=PetscFree(data->potential);CHKERRQ(ierr);
  ierr=PetscFree(data->source);CHKERRQ(ierr);
  ierr=PetscFree(data->density);CHKERRQ(ierr);
	
  return ierr;
}

/************************************************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "permutationmatrix"
/************************************************************************************************************************/

PetscErrorCode permutationmatrix(MPI_Comm com,int egnn, int totnn, int* petscordering, int* appordering, Mat* M )
{
	
  int                    ii,jj,kk,ll;
  PetscReal              one=1.0;
  PetscErrorCode         ierr;
	
  ii=petscordering[0];
  jj=petscordering[egnn-1];
	
  MatCreateMPIAIJ(com,egnn,egnn,
		  totnn,totnn,
		  1,PETSC_NULL,
		  1,PETSC_NULL,M);
	
  ll=0;
  for (kk=ii;kk<=jj;kk++){
    ierr = 	MatSetValue(*M,kk,appordering[ll++],one,INSERT_VALUES); CHKERRQ(ierr);
  }
	
  ierr = 		MatAssemblyBegin(*M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = 		MatAssemblyEnd(*M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	
  return ierr;
	
}

/************************************************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "compshg"
/************************************************************************************************************************/

PetscErrorCode compshg(PetscReal x0,PetscReal y0,PetscReal x1,PetscReal y1,
		       PetscReal x2,PetscReal y2, PetscReal *shgx, PetscReal *shgy){
	
  PetscReal denom;
	
  denom = (-(x1*y0) 
	   + (x2*y0) 
	   + (x0*y1) 
	   - (x2*y1) 
	   - (x0*y2) 
	   + (x1*y2));
	
  shgx[0]  =  (y1-y2)/denom;
  shgy[0]  = -(x1-x2)/denom;
  shgx[1]  = -(y0-y2)/denom;
  shgy[1]  =  (x0-x2)/denom;
  shgx[2]  =  (y0-y1)/denom;
  shgy[2]  = -(x0-x1)/denom;
	
  return 0;
}



/************************************************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "compsg"
/************************************************************************************************************************/
PetscErrorCode compsg(Tmesh *mesh, Tdata *data, Mat *SG,Mat *M, PetscReal epsilon)
{
  int                    iel,inode,jnode,kk;
  PetscInt               lv[3],LV[3];
  PetscReal              locmat[9],locdiffmat[9],locvec[3];
#define locmat(inode,jnode)  locmat[jnode + 3 * inode] 
#define locdiffmat(inode,jnode)  locdiffmat[jnode + 3 * inode] 
  PetscErrorCode         ierr;
  PetscMPIInt            rank,size;
	
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);CHKERRQ(ierr);
	
  /*
    Create matrix
  */	
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nCreate matrix\n****\n");CHKERRQ(ierr);
	
  ierr = MatCreateMPIAIJ(MPI_COMM_WORLD,mesh->egnn[rank],mesh->egnn[rank],
			 mesh->totnn,mesh->totnn,
			 1,PETSC_NULL,
			 1,PETSC_NULL,SG);CHKERRQ(ierr);
  
  ierr = MatCreateMPIAIJ(MPI_COMM_WORLD,mesh->egnn[rank],mesh->egnn[rank],
			 mesh->totnn,mesh->totnn,
			 1,PETSC_NULL,
			 1,PETSC_NULL,M);CHKERRQ(ierr); CHKERRQ(ierr);

  /*
    Add matrix and vector entries
  */
	
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nAdd matrix entries\n****\n");CHKERRQ(ierr);
	
  for (iel=0;iel<mesh->egen[rank];iel++){
		
    /*Local diffusion matrix*/
    for (inode=0;inode<3;inode++){
      for (jnode=0;jnode<=inode;jnode++){
	locdiffmat(inode,jnode) = mesh->triangles[iel].area * ( mesh->triangles[iel].shgx[inode] *
								mesh->triangles[iel].shgx[jnode] +
								mesh->triangles[iel].shgy[inode] *
								mesh->triangles[iel].shgy[jnode] );
      }
    }
		
    PetscLogFlops(3*6);
		
    PetscMemcpy (lv,mesh->triangles[iel].vertex,3*sizeof(PetscInt));
    AOApplicationToPetsc(mesh->reordering,3,lv);
		
    ISGlobalToLocalMappingApply(mesh->mapping,IS_GTOLM_DROP,
				3,mesh->triangles[iel].vertex,&kk,LV);
		
    /* Add transport term to local matrix*/
    locmat(0,1) = locdiffmat(1,0) * bern((data->potential[LV[1]]-data->potential[LV[0]])/epsilon);
    locmat(1,0) = locdiffmat(1,0) * bern((data->potential[LV[0]]-data->potential[LV[1]])/epsilon);
    locmat(0,2) = locdiffmat(2,0) * bern((data->potential[LV[2]]-data->potential[LV[0]])/epsilon);
    locmat(2,0) = locdiffmat(2,0) * bern((data->potential[LV[0]]-data->potential[LV[2]])/epsilon);
    locmat(1,2) = locdiffmat(2,1) * bern((data->potential[LV[2]]-data->potential[LV[1]])/epsilon);
    locmat(2,1) = locdiffmat(2,1) * bern((data->potential[LV[1]]-data->potential[LV[2]])/epsilon);
    locmat(0,0) = -(locmat(1,0) + locmat(2,0));
    locmat(1,1) = -(locmat(0,1) + locmat(2,1));
    locmat(2,2) = -(locmat(0,2) + locmat(1,2));
		
    PetscLogFlops(9);
		
    /*Local RHS*/		
    for (inode=0;inode<3;inode++){
      locvec[inode] = mesh->triangles[iel].area /3.0;//* data->source[LV[inode]] / 3.0 ;
    }
		
    //Impose BCs without symmetrization 
    for (inode=0;inode<3;inode++){
      if(mesh->triangles[iel].bccode[inode]==1){
	for(jnode=0;jnode<3;jnode++)
	  if(inode!=jnode){
	    locmat(inode,jnode)=0.0;
	  }
						
	//locvec[inode]=locmat(inode,inode)*data->density[LV[inode]];
	PetscLogFlops(1);
      }
    }
		
		
    for (inode=0;inode<3;inode++){      
      ierr = 	MatSetValue(*M,lv[inode],lv[inode],locvec[inode],ADD_VALUES); CHKERRQ(ierr);
      for (jnode=0;jnode<3;jnode++){
	ierr = 	MatSetValue(*SG,lv[inode],lv[jnode],locmat(inode,jnode),ADD_VALUES); CHKERRQ(ierr);
      }
    }
		
  }
  /*
    Assemble matrix
  */
  ierr = PetscPrintf(MPI_COMM_WORLD,"\n\n****\nAssemble matrices\n****\n");CHKERRQ(ierr);
	
  ierr = 		MatAssemblyBegin(*SG,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = 		MatAssemblyEnd(*SG,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	
	
  ierr = 		MatAssemblyBegin(*M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = 		MatAssemblyEnd(*M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);	
	
  return ierr;
}

/************************************************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "bern"
/************************************************************************************************************************/

PetscReal      bern             (PetscReal x){
	
  PetscReal              xlim=0.01,ax=PetscAbs(x),
    zero=0.0,bp,bignumber=80.0,one=1.0,fp,
    df,segno,ii;
	
  PetscErrorCode         ierr;
	
	
  ierr = PetscLogEventBegin(USER_EVENT,0,0,0,0);CHKERRQ(ierr);
	
  /*
    Compute Bernoulli function
    for 0 argument
  */	
	
  if (ax == zero){
    bp=one;
  }
	
  /*
    Compute Bernoulli function
    for big arguments using 
    asymptotic values
  */		
	
  else if (ax > bignumber){
    if (x > zero){
      bp=zero;
    } else {
      bp=-x;
    }
  } 
	
  /*
    Compute Bernoulli function
    for intermediate arguments using 
    the definition
  */
	
  else if (ax > xlim){
    bp=x/(exp(x)-one);
    PetscLogFlops(10);
  }
	
  /*
    Compute Bernoulli function
    for small arguments using Taylor
    expansion
  */
	
  else {
    ii=one;fp=one;df=one;segno=one;
		
    while ((PetscAbs(df) + one) > one)
      {
	ii=ii+one;
	segno=-segno;
	df=df*ax/ii;
	fp=fp+df;
	bp=one/fp;
	PetscLogFlops(1);
      }
    if (x < 0)
      bp = x + bp;
  }
	
	
  ierr = PetscLogEventEnd(USER_EVENT,0,0,0,0);CHKERRQ(ierr);
	
  return bp;
}	

/************************************************************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "monitor"
/************************************************************************************************************************/

PetscErrorCode monitor          (TS ts,PetscInt steps,PetscReal time,Vec x,void *ctx)
{
  char        filestep[50];
  PetscViewer viewer;
  Vec         out,out2;
  Tmesh       *mesh;

  PetscPrintf(MPI_COMM_WORLD,
	      "timestep %d\n",steps);

  mesh = ((Tmonitor*)ctx)->mesh;

  sprintf (filestep,"%sstep_%5.5d.m",((Tmonitor*)ctx)->pref,steps);
  PetscViewerASCIIOpen(MPI_COMM_WORLD,filestep,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);

  VecDuplicate(x,&out);
  VecDuplicate(x,&out2);

  sprintf (filestep,"step%5.5d",steps);
  PetscObjectSetName((PetscObject) out2,filestep);
  KSPSolve(((Tmonitor*)ctx)->ksp,x,out);
  MatMultTranspose( ((Tmesh*)mesh)->P,out,out2);

  VecView(out2,viewer);
  VecDestroy(out);
  VecDestroy(out2);
  PetscViewerDestroy(viewer);

  return 0;
}
