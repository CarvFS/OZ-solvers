/*
Program which solves Orstein-Zernike equation using the Picard iteration method
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

double* Epot(double x[], double t);
double* FT(double x[], double y[], double fun[]);
double* FTinv(double x[], double y[], double fun[]);
double* Bridge(double x[]);
double* ozs(double r[], double q[], double rho, double gammai[], double V[]);

// defining global constants and parameters
double pi = 4.0*atan(1.0);
#define N 800 // Generates N points
#define rmax 15.0 // Defining maximum coordinate values
#define dr (rmax/N) // r intervals
#define dk (pi/rmax) // q intervals

int main(int argc, char* argv[]){

    int i, j;
    // generating coordinates r and q
    double r[N], q[N];

    for(i = 0; i < N; i++){
        r[i] = (i+1)*dr;
        q[i] = (i+1)*dk;
    }
    
    // defining liquid parameters (temperature and density)
    double rho, T;

    rho = 0.9;
    T = 2.0;

    // getting the potential energy function at r values
    double *V;
    V = Epot(r, T);

    // initial guess for gamma function
    double gammai[N];
    memset(gammai, 0, sizeof gammai);

    // iterating
    double tol = 2, *g_new, alfa, sum;
    while(tol > 1e-4){
        sum = 0.0;
        g_new = ozs(r, q, rho, gammai, V); // acquiring the new value for gamma(r)

        for(i = 0; i < N; i++){
            alfa = 1.0; // mixing parameter
            g_new[i] = alfa*g_new[i] + (1.0-alfa)*gammai[i];
            sum += (g_new[i] - gammai[i])*(g_new[i] - gammai[i]);
            gammai[i] = g_new[i];
        }

        tol = sqrt(sum);
        printf("tol = %f\n", tol);
    }

    return 0;
}

//////////////////////////////////// Functions to be used

double* Epot(double x[], double t){
    static double E[N];
    int i;
    
    for(i = 0; i < N; i++){
        E[i] = (4.0/t)*((1.0/pow(x[i],12.0)) - (1.0/pow(x[i], 6.0)));
    }
    
    return E;
}

double* FT(double x[], double y[], double fun[]){
    static double fq[N];
    int i, j;

    for(j = 0; j < N; j++){
        fq[j] = 0;
        for(i = 0; i < N; i++){
            fq[j] += (4.0*pi/y[j])*sin(y[j]*x[i])*fun[i]*x[i]*dr;
        }
    }
    
    return fq;
}

double* FTinv(double x[], double y[], double fun[]){
    static double f[N];
    int i, j;

    for(j = 0; j < N; j++){
        f[j] = 0;
        for(i = 0; i < N; i++){
            f[j] += 1.0/(x[j]*2.0*pi*pi)*sin(x[j]*y[i])*fun[i]*y[i]*dk;
        }
    }
    
    return f;
}

double* Bridge(double x[]){
    static double Br[N];
    int i;

    // HNC closure relation
    for(i = 0; i < N; i++){
        Br[i] = 0.0*x[i];
    }

    return Br;
}

double* ozs(double r[], double q[], double rho, double gammai[], double V[]){
    double C[N], *C_hat, g_hat[N], *g_novo, *Br;
    int i;
    // Calculating C
    Br = Bridge(r);
    for(i = 0; i < N; i++){
        C[i]=exp(gammai[i] - V[i] + Br[i]) - gammai[i] - 1.0;
    }
    
    // Calculating Fourier transform of C, C_hat
    C_hat = FT(r,q,C);
    
    // Acquring Fourier transform of gamma function, g_hat, from OZ equation
    for(i = 0; i < N; i++){
        g_hat[i]=(rho*pow(C_hat[i],2))/(1.0-rho*C_hat[i]);
    }
    
    // Calculating inverse Fourier transform of g_hat to retrieve gamma(r)
    g_novo = FTinv(r,q,g_hat);

    return g_novo;
}
