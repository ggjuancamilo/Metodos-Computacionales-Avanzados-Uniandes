#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funciones.h"




int main(int argc, char **argv)
{
	
	double *U = malloc(POINTS3*USIZE*sizeof(double));
	double *UB = malloc(POINTS3*USIZE*sizeof(double));
	
	double *rho = malloc(POINTS3*sizeof(double));
	
	double *Fbx = malloc(POINTS3*USIZE*sizeof(double));
	double *Fby = malloc(POINTS3*USIZE*sizeof(double));
	double *Fbz = malloc(POINTS3*USIZE*sizeof(double));
	double *Fax = malloc(POINTS3*USIZE*sizeof(double));
	double *Fay = malloc(POINTS3*USIZE*sizeof(double));
	double *Faz = malloc(POINTS3*USIZE*sizeof(double));
	double *C = malloc(POINTS3*USIZE*sizeof(double));
	
	
	double h=0.0;
	int a =5;
	int p=0;
	init(U, UB);
	double dt;
	while (h<11.0)
	{	
		dt = lax(U, UB, Fbx, Fby, Fbz, Fax, Fay, Faz, C);
		//printf("%d\n", 5);
		//de aqui en adelante imprime para los 3 casos
		h = take_data(U,a,h);
		if (h!=0.0)
		{	
			int i,j,k;
			double rhoff;
			calc_var(U, rho);
			
			// p se usa para saber en quee radio voy
			p++;
			
			
			FILE *salida = fopen("data1", "w+");
					
			
			if (p==2){
				FILE *salida = fopen("data2", "w+");
			}
			if (p==3){
				FILE *salida = fopen("data3", "w+");
			}
			for(k=0; k < POINTS; k++)
			{
				for (j = 0; j < POINTS; j++)
				{
					for (i = 0; i < POINTS; i++)
					{
						rhoff = rho[matrix(i,j,k)];
						
						fprintf(salida, "%d %d %d %f \n", i, j, k, rhoff);
					}
					
				}
				
			}
			fclose(salida);
			if (p==1)
			{
				a = 30;
			}
			else if (p==2)
			{
				a = 60;
			}	
		}
		
		//printf("%f \n", dt);
	}
	
	
	
	
	return 0;
	
	
		
}

