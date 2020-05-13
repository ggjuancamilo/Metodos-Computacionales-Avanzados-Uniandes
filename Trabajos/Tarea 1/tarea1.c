#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"


/*Constantes*/
#define N 64
#define T 5.0 * pow(N,2.2) 
#define BETA 1.0
#define PI 3.14159265359
#define DT 0.005
#define DX 1.0

double Ak(int k, double *x);
double Ek(int k, double Ak, double dAk);
double F(double u, double uf, double ub);


int main(int argc, char **argv){
  
  int procesadores = 8;
  if(argc == 2){
    procesadores = atoi(argv[1]);
  }

  double *u = malloc(N*sizeof(double));/*posiciones y temporal de posciciones para hacer calculos*/
  double *utemp = malloc(N*sizeof(double));
  double *v = malloc(N*sizeof(double));/*derivada de las posiciones, utilizadas para el metodo leap frog*/
  double *vmed = malloc(N*sizeof(double));
  double *A1 = malloc(1000*sizeof(double));
  double *A2 = malloc(1000*sizeof(double));
  double *A3 = malloc(1000*sizeof(double));
  double *E1 = malloc(1000*sizeof(double));/*energia para modos de oscilacion k=01,2,3*/
  double *E2 = malloc(1000*sizeof(double));
  double *E3 = malloc(1000*sizeof(double));
  double *tiempo = malloc(1000*sizeof(double));
  

  /*Condiciones iniciales*/
  int i;
  for(i=0; i<N; i++){
    u[i] = sin(PI*i/(N-1));
    utemp[i] = sin(PI*i/(N-1));
    v[i]=0.0;
 
  }

  /*realiza las iteraciones, Falta guardar cada xxxx iteraciones y revisar las V inicial y final*/
  double t = 0.0;
  int k=(int)floor((T/DT)/999.0);
  int h=0;
  

  for(t=0.0;t < T; t = t+DT){
    if((int)(t/DT) == k*h && h<1000){
  
      A1[h]=Ak(1,u);
      A2[h]=Ak(2,u);
      A3[h]=Ak(3,u);      
      E1[h]=Ek(1,Ak(1,u),Ak(1,v));
      E2[h]=Ek(2,Ak(2,u),Ak(2,v));
      E3[h]=Ek(3,Ak(3,u),Ak(3,v));
      tiempo[h]=t;
      h+=1;
      }    
 
  omp_set_num_threads(procesadores);
#pragma omp parallel for private(i), shared(vmed,u,utemp)
    for(i=1; i<(N-2); i++){
      vmed[i] = v[i]+(F(u[i], u[i+1], u[i-1])*DT/2.0);
      utemp[i]=u[i]+vmed[i]*DT;
     
    }

  omp_set_num_threads(procesadores);
#pragma omp parallel for private(i), shared(v,utemp)
    for(i=1; i<(N-2); i++){
      v[i]=vmed[i]+(F(utemp[i], utemp[i+1], utemp[i-1])*DT/2.0);

    }
    /*Actualizo de los temp a los finales*/

    for(i=1; i<(N-2); i++){
      u[i]=utemp[i];
    }
    /*    printf("%f\n", t);*/

  }
 

  for(i=0; i < 1000;i++){
    printf("%f %f %f %f\n", tiempo[i], E1[i], E2[i], E3[i]);

  }
  return 0;

}
/*ver wiki*/
double Ak(int k, double *x){
  int i;
  double sum = 0.0;
  for(i=0;i<N;i++ ){
    sum += x[i]*sin(k*i*PI/N);
  }
  return pow(2.0/N,0.5)*sum;
}

/*Energias del modo k*/
double Ek(int k, double Ak, double dAk){
  double omega2 = 4.0*pow(sin(k*PI/(2*N)), 2.0);
  return (pow(dAk, 2.0) + (omega2*pow(Ak, 2.0)))/2;
}

/*segunda derivada de la posicion respecto al tiempo bajo la def del problema, uf=uforward, ub=ubackward*/
double F(double u, double uf, double ub){
  return (uf-(2.0*u)+ub)+(BETA*(pow(uf-u, 2.0)-pow(u-ub, 2.0)));
}



/*https://en.wikipedia.org/wiki/Leapfrog_integration,http://www.scholarpedia.org/article/Fermi-Pasta-Ulam_nonlinear_lattice_oscillations*/
