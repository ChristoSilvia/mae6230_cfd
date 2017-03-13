#include <iostream>
#include <math.h>


using namespace std;

/* Global variables */
int nx, ny, N;// Dimensions of the mesh
double rho0;// Density scale
double U;// Velocity scale
double H;//Length scale


/* FUNCTION: initial conditions.
    Defines the initial conditions for density, velocity and momentum.
*/
void setInitialConditions(double *rho, double *u, double *v, double *rhou, double *rhov)
{
    for(int k = 0; k < N; k++)
    {
        // Uniform density and zero initial momentum.
        rho[k] = rho0;
        u[k] = 0.;
        v[k] = 0.;
        rhou[k] = 0.;
        rhov[k] = 0.;
    }
}


/* FUNCTION: mesh generator.
*/
void createMesh(double *x, double dx, double *y, double dy)
{
    for(int i = 0; i<nx; i++)
        x[i] = i*dx;

    for(int j = 0; j<ny; j++)
        y[j] = j*dy;
}

/* FUNCTION: index translator.
    Auxiliar function that maps indices from a 2D matrix to the
    index of a 1D array formed by appending the rows of the matrix.

    i = matrix row index
    j = matrix column index
*/
int ind(int i, int j)
{
    return i + nx*j;
}

int main()
{
    // Define parameters for the mesh and create it
    nx = 3, ny = 3, N = nx*ny;
    rho0 = 1.;
    H = 1.0;

    double dx = H/(nx-1);
    double dy = H/(ny-1);
    double *x = new double[nx];
    double *y = new double[ny];

    createMesh(x, dx, y, dy);

    // Physical variables at time t_n, stored as 1D arrays
    double *rho = new double[N];
    double *u = new double[N];
    double *v = new double[N];
    double *rhou = new double[N];
    double *rhov = new double[N];

    // Give initial values to start the simulation
    setInitialConditions(rho, u, v, rhou, rhov);

    // Start simulation
    int nt = 100; /* This number should be chosen so that the flow does not change anymore */

    // Memory to save the predictor step (star variables) (p for predictor)
    double *rhop = new double[N];
    double *up = new double[N];
    double *vp = new double[N];
    double *rhoup = new double[N];
    double *rhovp = new double[N];

    for(int n = 0; n < nt; n++)
    {
        // From the initial conditions, compute the predictor for the interior points
        for(int i = 1; i < nx-1; i++)
        {
            for(int j = 1; j < ny-1; j++)
            {
                rhop[ind(i, j)] = rho[ind(i, j)];// formulas in construction
                rhoup[ind(i, j)] = rhou[ind(i, j)];
                rhovp[ind(i, j)] = rhov[ind(i, j)];
            }
        }

        // and for the boundary points

        // get the predicted velocity from the predicted momentum

        // With the predictor calculated, compute the corrector for the interior points

        // and for the boundary points

        // get the corrected velocity from the corrected momentum
    }

    // Release heap memory...
    delete[] x;
    delete[] y;
    delete[] rho;
    delete[] u;
    delete[] v;
    delete[] rhou;
    delete[] rhov;
    delete[] rhop;
    delete[] up;
    delete[] vp;
    delete[] rhoup;
    delete[] rhovp;

    // ... and exit program.
    return 0;
}
