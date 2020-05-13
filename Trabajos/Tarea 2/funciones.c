#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funciones.h"

//Se usa para notacion matricial en 4 dimensiones
int tensor(int i, int j, int k, int l){
	return l +  USIZE * i + USIZE * POINTS * j + USIZE * POINTS2 * k;
}

//Notación para matrices en 3 dimensiones
int matrix(int i, int j, int k){
	return  i + POINTS * j + POINTS2 * k;
	
}

//calcula el termino e de la energia
double calc_e(double *grid, int i, int j, int k){
	return (1.0/pow(grid[tensor(i,j,k,0)],2.0))*(grid[tensor(i,j,k,0)]*grid[tensor(i,j,k,4)]-0.5*(pow(grid[tensor(i,j,k,1)],2.0)+pow(grid[tensor(i,j,k,2)],2.0)+pow(grid[tensor(i,j,k,3)],2.0)));
}

void init(double *U, double *UB){
	
	int i,j,k,l;
	
	for	(i=0; i<POINTS;i++){
	
		for (j=0; j<POINTS; j++){
		
			for (k=0; k<POINTS; k++){
				
				U[tensor(i,j,k,0)]=1.0;
				U[tensor(i,j,k,1)]=0.0;
				U[tensor(i,j,k,2)]=0.0;
				U[tensor(i,j,k,3)]=0.0;
				U[tensor(i,j,k,4)]=25E-6;
			
			}
		}
	}
	
	U[tensor(64,64,64,4)]=10;
	
	for	(j=0; j<USIZE;j++){
	
		for (k=0; k<USIZE; k++){
		
			for (l=0; l<USIZE; l++){
			
				UB[tensor(0,j,k,l)]=0.0;
				UB[tensor(127,j,k,l)]=0.0;
							
			}
		}
	}
}

void calc_Fx(double *U, double *F){
	
	int i,j,k,l;
	
	for(i=0;i<POINTS;i++){
		
		for(j=0;j<POINTS;j++){
			
			for(k=0;k<POINTS;k++){
			
				F[tensor(i,j,k,0)] = U[tensor(i,j,k,1)];
				F[tensor(i,j,k,1)] = (pow(U[tensor(i,j,k,1)],2.0)/U[tensor(i,j,k,0)])+U[tensor(i,j,k,0)]*calc_e(U,i,j,k)*(GAMMA+1.0);
				F[tensor(i,j,k,2)] = U[tensor(i,j,k,1)]*U[tensor(i,j,k,2)]/U[tensor(i,j,k,0)];
				F[tensor(i,j,k,3)] = U[tensor(i,j,k,1)]*U[tensor(i,j,k,3)]/U[tensor(i,j,k,0)];
				F[tensor(i,j,k,4)] = U[tensor(i,j,k,1)]*((U[tensor(i,j,k,4)]/U[tensor(i,j,k,0)])+calc_e(U,i,j,k)*(GAMMA+1.0));
			
			}
		}
	}
}

void calc_Fy(double *U, double *F){
	int i,j,k,l;
	
	for (k = 0; k < POINTS; k++)
	{
		for (j = 0; j <POINTS ; j++)
		{
			for (i = 0; i < POINTS; i++)
			{
				F[tensor(i,j,k,0)] = U[tensor(i,j,k,2)];
				F[tensor(i,j,k,1)] = U[tensor(i,j,k,1)]*U[tensor(i,j,k,2)]/U[tensor(i,j,k,0)];
				F[tensor(i,j,k,2)] = (pow(U[tensor(i,j,k,2)],2.0)/U[tensor(i,j,k,0)])+U[tensor(i,j,k,0)]*calc_e(U,i,j,k)*(GAMMA+1.0);
				F[tensor(i,j,k,3)] = U[tensor(i,j,k,2)]*U[tensor(i,j,k,3)]/U[tensor(i,j,k,0)];
				F[tensor(i,j,k,3)] = U[tensor(i,j,k,2)]*((U[tensor(i,j,k,4)]/U[tensor(i,j,k,0)])+calc_e(U,i,j,k)*(GAMMA+1.0));
			}
			
		}
		
	}
	
	
	
}

//Calcula funcion F en coordenada Z
void calc_Fz(double *U, double *F){
	
	
	double e;
	int i;
	int j;
	int k;
	int l;
	//U[tensor(i,j,k,0)]
	
	for (k = 0; k < POINTS ; k++)
	{
		for (j = 0; j < POINTS; j++)
		{
			for (i = 0; i < POINTS	; i++)
			{
				F[tensor(i,j,k,0)] = U[tensor(i,j,k,3)]  ;
				F[tensor(i,j,k,1)] = U[tensor(i,j,k,3)] * U[tensor(i,j,k,1)] / U[tensor(i,j,k,0)];
				F[tensor(i,j,k,2)] = U[tensor(i,j,k,3)] * U[tensor(i,j,k,2)] / U[tensor(i,j,k,0)];
				F[tensor(i,j,k,3)] = pow(U[tensor(i,j,k,3)],2) / U[tensor(i,j,k,0)] + U[tensor(i,j,k,0)] * calc_e(U,i,j,k) * U[tensor(i,j,k,0)] * ( GAMMA + 1) ;
				F[tensor(i,j,k,4)] = U[tensor(i,j,k,3)] * (U[tensor(i,j,k,4)]/U[tensor(i,j,k,0)] + calc_e(U,i,j,k)*(GAMMA -1)) ;
				
			}
			
		}
		
	}
	
}

//calcula Fi+1/2
void calc_FA(double *U, double *UB, double *FAx, double  *FAy, double  *FAz ){
	
	int i;
	int j;
	int k;	
	int l;

	for (k = 0; k < POINTS; k++)
	{
		for (j = 0; j < POINTS; j++)
		{
			for (i = 0; i < POINTS - 1 ; i++)
			{
				for (l = 0; l < USIZE ; l++)
				{
					UB[tensor(i,j,k,l)] = 0.5  *( U[tensor(i+1,j,k,l)] + U[tensor(i,j,k,l)]);
				}
				
			}
			
		}
		
	}
	calc_Fx(UB,FAx);
	
	for (k = 0; k < POINTS; k++)
	{
		for (j = 0; j < POINTS-1; j++)
		{
			for (i = 0; i < POINTS ; i++)
			{
				for (l = 0; l < USIZE; l++)
				{
					UB[tensor(i,j,k,l)] = 0.5  *( U[tensor(i,j+1,k,l)] + U[tensor(i,j,k,l)]);
				}
				
			}
			
		}
		
	}
	
	calc_Fy(UB,FAy);
	
	for (k = 0; k < POINTS-1; k++)
	{
		for (j = 0; j < POINTS; j++)
		{
			for (i = 0; i < POINTS ; i++)
			{
				for (l = 0; l < USIZE; l++)
				{
					UB[tensor(i,j,k,l)] = 0.5  *( U[tensor(i,j,k+1,l)] + U[tensor(i,j,k,l)]);
				}
				
			}
			
		}
		
	}
	calc_Fz(UB,FAz);	
}

//calcula Fi-1/2
void calc_FB(double *U, double *UB, double *FBx, double  *FBy, double  *FBz ){
	
	
	int i;
	int j;
	int k;	
	int l;

	for (k = 0; k < POINTS; k++)
	{
		for (j = 0; j < POINTS; j++)
		{
			for (i = 1; i < POINTS; i++)
			{
				for (l = 0; l < USIZE; l++)
				{
					UB[tensor(i,j,k,l)] = 0.5  *( U[tensor(i,j,k,l)] + U[tensor(i-1,j,k,l)]);
				}
				
			}
			
		}
		
	}
	calc_Fx(UB,FBx);
	for (k = 0; k < POINTS; k++)
	{
		for (j = 1; j < POINTS; j++)
		{
			for (i = 0; i < POINTS; i++)
			{
				for (l = 0; l < USIZE; l++)
				{
					UB[tensor(i,j,k,l)] = 0.5  *( U[tensor(i,j,k,l)] + U[tensor(i,j-1,k,l)]);
				}
				
			}
			
		}
		
	}
	calc_Fy(UB,FBy);
	for (k = 1; k < POINTS; k++)
	{
		for (j = 0; j < POINTS; j++)
		{
			for (i = 0; i < POINTS; i++)
			{
				for (l = 0; l < USIZE; l++)
				{
					UB[tensor(i,j,k,l)] = 0.5  *( U[tensor(i,j,k,l)] + U[tensor(i,j,k-1,l)]);
				}
				
			}
			
		}
		
	}
	calc_Fz(UB,FBz);	
}

double calc_vmax(double *U,double *c){
  int i;
  int j;
  int k;
  double max1;
  double max2;
  double tempx=0.0;
  double tempy=0.0;
  double tempz=0.0; 
  double temp2x=0.0; 
  double umax=0.0;

  for ( k = 0; k < POINTS; k++)
    {
      for (j = 0; j <  POINTS ; j++)
	  {
		for (i = 0; i < POINTS ; i++)
		{
			tempx = U[tensor(i,j,k,1)]/U[tensor(i,j,k,0)];
			temp2x = c[matrix(i,j,k)];
			tempy = U[tensor(i,j,k,2)]/U[tensor(i,j,k,0)];
			tempz = U[tensor(i,j,k,3)]/U[tensor(i,j,k,0)];
			if (tempx>tempy)
			{
				max1=tempx;
			}
			else
			{
				max1 = tempy;
			}

			if (max1>tempz)
			{
				max2=tempz;
			}
			else
			{
				max2 = tempz;
			}			
	
			
			if (max2+temp2x > umax)
			{
				umax = max2+temp2x;	
			}
		}		  
	  }
    }
  return umax;
}

void cspeed(double *U, double *c){
  
  double rh;
  double pres;
  int i;
  int j;
  int k;
  
  for (k = 0; k < POINTS; k++)
  {
	  for (j = 0; j < POINTS; j++)
	  {
		  for (i = 0; i < POINTS; i++)
		  {
			  rh = U[tensor(i,j,k,0)];
			  pres = (GAMMA - 1.0) * U[tensor(i,j,k,0)] * calc_e(U,i,j,k) ;
			  c[matrix(i,j,k)] = pow(GAMMA * pres/rh ,0.5);	
		  }  
	  }  
  }
}

void calc_UA(double *U, double *UB, double *Fbx, double *Fby, double *Fbz,double *Fax, double *Fay, double *Faz, double dt){
		 
	int i;
	int j;
	int k; 
	int l;
	for (k = 0; k < POINTS ; k++)
	{
		for (j = 0; j < POINTS ; j++)
		{
			for (i = 0; i < POINTS	; i++)
			{
				for (l = 0; l < USIZE; l++)
				{
					UB[tensor(i,j,k,l)] = U[tensor(i,j,k,l)] + (dt/DELTA)*( (Fbx[tensor(i,j,k,l)]-Fax[tensor(i,j,k,l)])+(Fby[tensor(i,j,k,l)]-Fay[tensor(i,j,k,l)])+(Fbz[tensor(i,j,k,l)]-Faz[tensor(i,j,k,l)]) )  ;
				}
				
			}
			
		}
		
	}
	
	for (k = 0; k < POINTS ; k++)
	{
		for (j = 0; j < POINTS ; j++)
		{
			for (i = 0; i < POINTS	; i++)
			{
				for (l = 0; l < USIZE; l++)
				{
					U[tensor(i,j,k,l)] = UB[tensor(i,j,k,l)];
				}
				
			}
			
		}
		
	}
	
	
	
}

//funcion usada para parar lax
int take_data(double *U, int a, double h){
	int i,j,k,l;
	
	if (U[tensor(64+a,64,64,1)] > 0.0 || U[tensor(64-a,64,64,1)] > 0.0)
	{	
		if ( a<120 )
		{
			h = 9.0;
		}
		else
		{
			h=12.0;
		}
				
	}
	else if (U[tensor(64,64+a,64,2)] > 0.0 || U[tensor(64,64-a,64,2)] > 0.0)
	{
		if ( a<60 )
		{
			h = 9.0;
		}
		else
		{
			h=12.0;
		}			
	}
	else if(U[tensor(64,64,64-a,2)] > 0.0 || U[tensor(64,64+a,64,2)] > 0.0)
	{
		if ( a<60 )
		{
			h = 9.0;
		}
		else
		{
			h=12.0;
		}
	}
	return h;
}

double lax(double *U, double *UB, double *Fbx, double *Fby, double *Fbz,double *Fax, double *Fay, double *Faz, double *C){
	
	double velmax;
	double dt = 1E-4;
	
	calc_FB(U,UB,Fbx,Fby,Fbz);
	calc_FA(U,UB,Fax,Fay,Faz);
	calc_UA(U,UB,Fbx,Fbx,Fbz,Fax,Fay,Faz,dt);
	cspeed(U,C);
	velmax = calc_vmax(U,C);
	dt = 0.5 * DELTA/velmax;
	return dt;
}

void calc_var(double *U, double *rho){
	
	int i;
	int j;
	int k;
	
	for (k = 0; k < POINTS; k++)
	{
		for (j = 0; j < POINTS; j++)
		{
			for (i = 0; i < POINTS; i++)
			{
				rho[matrix(i,j,k)] = U[tensor(i,j,k,0)];
			}	
		}
	}
}
