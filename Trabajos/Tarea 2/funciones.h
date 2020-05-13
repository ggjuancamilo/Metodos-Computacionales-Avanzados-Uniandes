#define USIZE 5
#define GAMMA  1.4
#define POINTS 128
#define POINTS2 128*128
#define POINTS3 128*128*128
#define DELTA 2.0


int tensor(int i, int j, int k, int l);
int matrix(int i, int j, int k);
double calc_e(double *grid, int i, int j, int k);
void init(double *U, double *UB);
void calc_Fx(double *U, double *F);
void calc_Fy(double *U, double *F);
void calc_Fz(double *U, double *F);
void calc_FA(double *U, double *UB, double *FAx, double  *FAy, double  *FAz );
void calc_FB(double *U, double *UB, double *FBx, double  *FBy, double  *FBz );
double calc_vmax(double *U,double *c);
void cspeed(double *U, double *c);
void calc_UA(double *U, double *UB, double *Fbx, double *Fby, double *Fbz,double *Fax, double *Fay, double *Faz, double dt);
int take_data(double *U, int a, double h);
double lax(double *U, double *UB, double *Fbx, double *Fby, double *Fbz,double *Fax, double *Fay, double *Faz, double *C);
void calc_var(double *U, double *rho);
