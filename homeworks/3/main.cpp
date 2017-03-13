#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>// NEEDED FOR PADDING
#include <sstream>// NEEDED FOR PADDING

using namespace std;

/* Global variables */
int nx, ny, N;// Dimensions of the mesh
double rho0;// Density scale
double U;// Velocity scale
double H;//Length scale


std::ofstream data;// Data for each time step will be saved here.
std::string filename;


/* FUNCTION: fill filenames with zeroes
*/
std::string ZeroPadNumber(int num)
{
    std::ostringstream ss;
    ss << std::setw(5) << std::setfill('0') << num;
    return ss.str();
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

/* FUNCTION: initial conditions.
    Defines the initial conditions for density and momentum.
    This includes the boundaries.
*/
void setInitialConditions(double *rho, double *rhou, double *rhov, double *u, double *v)
{
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            rho[ind(i, j)] = rho0;

            if(j == ny-1) rhou[ind(i, j)] = rho0*U;
            else rhou[ind(i, j)] = 0.;

            rhov[ind(i, j)] = 0;

            u[ind(i, j)] = rhou[ind(i, j)]/rho[ind(i, j)];
            v[ind(i, j)] = rhov[ind(i, j)]/rho[ind(i, j)];
        }
    }
}


/* FUNCTION: mesh generator.
*/
void createMesh(double *x, double dx, double *y, double dy)
{
    for(int i = 0; i < nx; i++)
        x[i] = i*dx;

    for(int j = 0; j < ny; j++)
        y[j] = j*dy;
}

/* MAIN FUNCTION
*/
int main()
{
    // Define parameters for the mesh and create it
    nx = 8, ny = 8, N = nx*ny;
    rho0 = 1.0;
    H = 1.0;
    U = 1.0;
    double Re = 100*rho0*H*U;
    double Ma = 0.1*U;

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
    filename = "data/data_" + ZeroPadNumber(0) + ".csv";
    data.open(filename.c_str());
    for(int i = 0; i < nx; i++)
        for(int j = 0; j < ny; j++)
            data << x[i] << "\t" << y[j] << "\t" << rho[ind(i, j)] << "\t" << u[ind(i, j)] << "\t" << v[ind(i, j)] << "\n";
    data.close();

    // Start simulation
    int nt = 3; /** This number should be chosen so that the flow does not change anymore **/
    double sigma = 0.5;
    double dt = sigma/( (U/dx) + (U/dy) + (1/Ma)*sqrt( pow(dx,-2) + pow(dy,-2) ) )/( 1. + 2.*Re*min(dx,dy) );

    double a1 = dt/dx;
    double a2 = dt/dy;
    double a3 = dt/(Re*dx*dx);
    double a4 = dt/(Re*dy*dy);
    double a5 = dt/(12.*Re*dx*dy);

    // Memory to store the predictor step (star variables)
    double *rhop = new double[N];
    double *up = new double[N];
    double *vp = new double[N];
    double *rhoup = new double[N];
    double *rhovp = new double[N];

    for(int n = 0; n < nt; n++)
    {
        // From the initial conditions, compute the predictor for the interior points,...
        for(int i = 1; i < nx-1; i++)
        {
            for(int j = 1; j < ny-1; j++)
            {
                rhop[ind(i, j)] = rho[ind(i, j)] - a1*( rhou[ind(i+1, j)] - rhou[ind(i, j)] ) - a2*( rhov[ind(i,j+1)] - rhov[ind(i, j)] );

                rhoup[ind(i, j)] = rhou[ind(i, j)];
                rhoup[ind(i, j)] -= a1*( rhou[ind(i+1, j)]*u[ind(i+1, j)] - rhou[ind(i, j)]*u[ind(i, j)] + (1./pow(Ma,2))*( rho[ind(i+1, j)] - rho[ind(i, j)] ) );
                rhoup[ind(i, j)] -= a2*( rhou[ind(i, j+1)]*v[ind(i, j+1)] - rhou[ind(i, j)]*v[ind(i, j)] );
                rhoup[ind(i, j)] += (4./3.)*a3*( u[ind(i+1, j)] - 2.*u[ind(i, j)] + u[ind(i-1, 1)] ) + a4*( u[ind(i, j+1)] - 2.*u[ind(i, j)] + u[ind(i, j-1)] );
                rhoup[ind(i, j)] += a5*( v[ind(i+1, j+1)] - v[ind(i+1, j-1)] - v[ind(i-1, j+1)] + v[ind(i-1, j-1)] );

                rhovp[ind(i, j)] = rhov[ind(i, j)];
                rhovp[ind(i, j)] -= a1*( rhov[ind(i+1, j)]*u[ind(i+1, j)] - rhov[ind(i, j)]*u[ind(i, j)] );
                rhovp[ind(i, j)] -= a2*( rhov[ind(i, j+1)]*v[ind(i, j+1)] - rhov[ind(i, j)]*v[ind(i, j)] + (1./pow(Ma,2))*( rho[ind(i, j+1)] - rho[ind(i, j)] ) );
                rhovp[ind(i, j)] += a3*( v[ind(i+1, j)] - 2.*v[ind(i, j)] + v[ind(i-1, j)] ) + (4./3.)*a4*( v[ind(i, j+1)] - 2.*v[ind(i, j)] + v[ind(i, j-1)] );
                rhovp[ind(i, j)] += a5*( u[ind(i+1, j+1)] - u[ind(i+1, j-1)] - u[ind(i-1, j+1)] + u[ind(i-1, j-1)] );
            }
        }

        // ... the left and right walls (not including the corners),...
        for(int j = 1; j < ny-1; j++)
        {
            rhop[ind(0, j)] = rho[ind(0, j)] - 0.5*a1*( - 3.*rhou[ind(0, j)] + 4.*rhou[ind(1, j)] - rhou[ind(2, j)] );
            rhoup[ind(0, j)] = 0.;
            rhovp[ind(0, j)] = 0.;

            rhop[ind(nx-1, j)] = rho[ind(nx-1, j)] - 0.5*a1*( rhou[ind(nx-3, j)] - 4.*rhou[ind(nx-2, j)] + 3.*rhou[ind(nx-1, j)] );
            rhoup[ind(nx-1, j)] = 0.;
            rhovp[ind(nx-1, j)] = 0.;
        }

        // ...and the bottom and top boundaries (not including the corners).
        for(int i = 1; i < nx-1; i++)
        {
            rhop[ind(i, 0)] = rho[ind(i, 0)] - 0.5*a2*( - 3.*rhov[ind(i, 0)] + 4.*rhov[ind(i, 1)] - rhov[ind(i, 2)] );
            rhoup[ind(i, 0)] = 0.;
            rhovp[ind(i, 0)] = 0.;

            rhop[ind(i, ny-1)] = rho[ind(i, ny-1)] - 0.5*a2*( rhov[ind(i, ny-3)] - 4.*rhov[ind(i, ny-2)] + 3.*rhov[ind(i, ny-1)] );
            rhop[ind(i, ny-1)] -= 0.5*a1*U*( rho[ind(i+1, ny-1)] - rho[ind(i-1, ny-1)] );
            rhoup[ind(i, ny-1)] = rhop[ind(i, ny-1)]*U;
            rhovp[ind(i, ny-1)] = 0.;
        }

        // For the predictor density at the corners, use weighted averages.
        rhop[ind(0, 0)] = 0.5*( rho[ind(0, 0)] - 0.5*a1*( - 3.*rhou[ind(0, 0)] + 4.*rhou[ind(1, 0)] - rhou[ind(2, 0)] ) );
        rhop[ind(0, 0)] += 0.5*( rho[ind(0, 0)] - 0.5*a2*( - 3.*rhov[ind(0, 0)] + 4.*rhov[ind(0, 1)] - rhov[ind(0, 2)] ) );

        rhop[ind(0, ny-1)] = 0.5*( rho[ind(0, ny-1)] - 0.5*a1*( - 3.*rhou[ind(0, ny-1)] + 4.*rhou[ind(1, ny-1)] - rhou[ind(2, ny-1)] ) );
        rhop[ind(0, ny-1)] += 0.5*( rho[ind(0, ny-1)] - 0.5*a2*( rhov[ind(0, ny-3)] - 4.*rhov[ind(0, ny-2)] + 3.*rhov[ind(0, ny-1)] ) );
        rhop[ind(0, ny-1)] -= 0.25*a1*U*( -3.*rho[ind(0, ny-1)] + 4.*rho[ind(1, ny-1)] - rho[ind(2, ny-1)] );

        rhop[ind(nx-1, 0)] = 0.5*( rho[ind(nx-1, 0)] - 0.5*a1*( rhou[ind(nx-3, 0)] - 4.*rhou[ind(nx-2, 0)] + 3.*rhou[ind(nx-1, 0)] ) );
        rhop[ind(nx-1, 0)] += 0.5*( rho[ind(nx-1, 0)] - 0.5*a2*( - 3.*rhov[ind(nx-1, 0)] + 4.*rhov[ind(nx-1, 1)] - rhov[ind(nx-1, 2)] ) );

        rhop[ind(nx-1, ny-1)] = 0.5*( rho[ind(nx-1, ny-1)] - 0.5*a2*( rhov[ind(nx-1, ny-3)] - 4.*rhov[ind(nx-1, ny-2)] + 3.*rhov[ind(nx-1, ny-1)] ) );
        rhop[ind(nx-1, ny-1)] -= 0.25*a1*U*( rho[ind(nx-3, ny-1)] - 4.*rho[ind(nx-2, ny-1)] + 3.*rho[ind(nx-1, ny-1)] );
        rhop[ind(nx-1, ny-1)] += 0.5*( rho[ind(nx-1, ny-1)] - 0.5*a1*( rhou[ind(nx-3, ny-1)] - 4.*rhou[ind(nx-2, ny-1)] + 3.*rhou[ind(nx-1, ny-1)] ) );

        // For the predictor momentum at the corners, jupt use the velocity
        rhoup[ind(0, 0)] = 0.;
        rhoup[ind(0, ny-1)] = rhop[ind(0, ny-1)]*U;
        rhoup[ind(nx-1, 0)] = 0.;
        rhoup[ind(nx-1, ny-1)] = rhop[ind(nx-1, ny-1)]*U;

        rhovp[ind(0, 0)] = 0.;
        rhovp[ind(0, ny-1)] = 0.;
        rhovp[ind(nx-1, 0)] = 0.;
        rhovp[ind(nx-1, ny-1)] = 0.;

        // Get the predictor velocity
        for(int i = 0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                up[ind(i, j)] = rhoup[ind(i, j)]/rhop[ind(i, j)];
                vp[ind(i, j)] = rhovp[ind(i, j)]/rhop[ind(i, j)];
            }
        }

        // With the predictor calculated, compute the corrector for the interior points...
        for(int i = 1; i < nx-1; i++)
        {
            for(int j = 1; j < ny-1; j++)
            {
                rho[ind(i, j)] = 0.5*( rho[ind(i, j)] + rhop[ind(i, j)] ) - 0.5*a1*(rhoup[ind(i, j)] - rhoup[ind(i, j-1)]) - 0.5*a2*( rhovp[ind(i, j)] - rhovp[ind(i, j-1)] );

                rhou[ind(i, j)] = 0.5*( rhou[ind(i, j)] + rhoup[ind(i, j)] );
                rhou[ind(i, j)] -= 0.5*a1*( rhoup[ind(i, j)]*up[ind(i, j)] - rhoup[ind(i-1, j)]*up[ind(i-1, j)] + (1./pow(Ma,2))*( rhop[ind(i, j)] - rhop[ind(i-1, j)] ) );
                rhou[ind(i, j)] -= 0.5*a2*( rhovp[ind(i, j)]*up[ind(i, j)] - rhovp[ind(i, j-1)]*up[ind(i, j-1)] );
                rhou[ind(i, j)] += (2./3.)*a3*( up[ind(i+1, j)] - 2.*up[ind(i, j)] + up[ind(i-1, j)] ) + 0.5*a4*( up[ind(i, j+1)] - 2.*up[ind(i, j)] + up[ind(i, j-1)] );
                rhou[ind(i, j)] += 0.5*a5*( vp[ind(i+1, j+1)] - vp[i+1, j-1] - vp[ind(i-1, j+1)] + vp[ind(i-1, j-1)] );

                rhov[ind(i, j)] = 0.5*( rhov[ind(i, j)] + rhovp[ind(i, j)] );
                rhov[ind(i, j)] -= 0.5*a1*( rhovp[ind(i, j)]*up[ind(i, j)] - rhovp[ind(i-1, j)]*up[ind(i-1, j)] );
                rhov[ind(i, j)] -= 0.5*a2*( rhovp[ind(i, j)]*vp[ind(i, j)] - rhovp[ind(i, j-1)]*vp[ind(i, j-1)] + (1./pow(Ma,2))*( rhop[ind(i, j)] - rhop[ind(i, j-1)] ) );
                rhov[ind(i, j)] += 0.5*a3*( vp[ind(i+1, j)] - 2.*vp[ind(i, j)] + vp[ind(i-1, j)] ) + (2./3.)*a4*( vp[ind(i, j+1)] - 2.*vp[ind(i, j)] + vp[ind(i, j-1)] );
                rhov[ind(i, j)] += 0.5*a5*( up[ind(i+1, j+1)] - up[ind(i+1, j-1)] - up[ind(i-1, j+1)] + up[ind(i-1, j-1)] );
            }
        }

        // ... the left and right walls (not including the corners),...
        for(int j = 1; j < ny-1; j++)
        {
            rho[ind(0, j)] = 0.5*( rho[ind(0, j)] + rhop[ind(0, j)] ) - 0.25*a1*( - 3.*rhoup[ind(0, j)] + 4.*rhoup[ind(1, j)] - rhoup[ind(2, j)] );
            rhou[ind(0, j)] = 0.;
            rhov[ind(0, j)] = 0.;

            rho[ind(nx-1, j)] = 0.5*( rho[ind(nx-1, j)] + rhop[ind(nx-1, j)] ) - 0.25*a1*( rhoup[ind(nx-3, j)] - 4.*rhoup[ind(nx-2, j)] + 3.*rhoup[ind(nx-1, j)] );
            rhou[ind(nx-1, j)] = 0.;
            rhov[ind(nx-1, j)] = 0.;
        }

        // ...and the bottom and top boundaries (not including the corners).
        for(int i = 1; i < nx-1; i++)
        {
            rho[ind(i, 0)] = 0.5*( rho[ind(i, 0)] + rhop[ind(i, 0)] ) - 0.25*a2*( - 3.*rhovp[ind(i, 0)] + 4.*rhovp[ind(i, 1)] - rhovp[ind(i, 2)] );
            rhou[ind(i, 0)] = 0.;
            rhov[ind(i, 0)] = 0.;

            rho[ind(i, ny-1)] = 0.5*( rho[ind(i, ny-1)] + rhop[ind(i, ny-1)] ) - 0.25*a2*( rhovp[ind(i, ny-3)] - 4.*rhovp[ind(i, ny-2)] + 3.*rhovp[ind(i, ny-1)] );
            rho[ind(i, ny-1)] -= 0.25*a1*U*( rhop[ind(i+1, ny-1)] - rhop[ind(i-1, ny-1)] );
            rhou[ind(i, ny-1)] = rho[ind(i, ny-1)]*U;
            rhov[ind(i, ny-1)] = 0.;
        }

        // For the corrector density at the corners, use weighted averages.
        rho[ind(0, 0)] = 0.5*( 0.5*( rho[ind(0, 0)] + rhop[ind(0, 0)] ) - 0.25*a1*( - 3.*rhoup[ind(0, 0)] + 4.*rhoup[ind(1, 0)] - rhoup[ind(2, 0)] ) );
        rho[ind(0, 0)] += 0.5*( 0.5*( rho[ind(0, 0)] + rhop[ind(0, 0)] ) - 0.25*a2*( - 3.*rhovp[ind(0, 0)] + 4.*rhovp[ind(0, 1)] - rhovp[ind(0, 2)] ) );

        rho[ind(0, ny-1)] = 0.5*( 0.5*( rho[ind(0, ny-1)] + rhop[ind(0, ny-1)] ) - 0.25*a1*( - 3.*rhoup[ind(0, ny-1)] + 4.*rhoup[ind(1, ny-1)] - rhoup[ind(2, ny-1)] ) );
        rho[ind(0, ny-1)] += 0.5*( 0.5*( rho[ind(0, ny-1)] + rhop[ind(0, ny-1)] ) - 0.25*a2*( rhovp[ind(0, ny-3)] - 4.*rhovp[ind(0, ny-2)] + 3.*rhovp[ind(0, ny-1)] ) );
        rho[ind(0, ny-1)] -= 0.125*a1*U*( -3.*rhop[ind(0, ny-1)] + 4.*rhop[ind(1, ny-1)] - rhop[ind(2, ny-1)] );

        rho[ind(nx-1, 0)] = 0.5*( 0.5*( rho[ind(nx-1, 0)] + rhop[ind(nx-1, 0)] ) - 0.25*a1*( rhoup[ind(nx-3, 0)] - 4.*rhoup[ind(nx-2, 0)] + 3.*rhoup[ind(nx-1, 0)] ) );
        rho[ind(nx-1, 0)] += 0.5*( 0.5*( rho[ind(nx-1, 0)] + rhop[ind(nx-1, 0)] ) - 0.25*a2*( - 3.*rhovp[ind(nx-1, 0)] + 4.*rhovp[ind(nx-1, 1)] - rhovp[ind(nx-1, 2)] ) );

        rho[ind(nx-1, ny-1)] = 0.5*( 0.5*( rho[ind(nx-1, ny-1)] + rhop[ind(nx-1, ny-1)] ) - 0.25*a2*( rhovp[ind(nx-1, ny-3)] - 4.*rhovp[ind(nx-1, ny-2)] + 3.*rhovp[ind(nx-1, ny-1)] ) );
        rho[ind(nx-1, ny-1)] -= 0.125*a1*U*( rhop[ind(nx-3, ny-1)] - 4.*rhop[ind(nx-2, ny-1)] + 3.*rhop[ind(nx-1, ny-1)] );
        rho[ind(nx-1, ny-1)] += 0.5*( 0.5*( rho[ind(nx-1, ny-1)] + rhop[ind(nx-1, ny-1)] ) - 0.25*a1*( rhoup[ind(nx-3, ny-1)] - 4.*rhoup[ind(nx-2, ny-1)] + 3.*rhoup[ind(nx-1, ny-1)] ) );

        // For the corrector momentum at the corners, just use the velocity.
        rhou[ind(0, 0)] = 0.;
        rhou[ind(0, ny-1)] = rho[ind(0, ny-1)]*U;
        rhou[ind(nx-1, 0)] = 0.;
        rhou[ind(nx-1, ny-1)] = rho[ind(nx-1, ny-1)]*U;

        rhov[ind(0, 0)] = 0.;
        rhov[ind(0, ny-1)] = 0.;
        rhov[ind(nx-1, 0)] = 0.;
        rhov[ind(nx-1, ny-1)] = 0.;

        // Get the corrector velocity
        for(int i = 0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                u[ind(i, j)] = rhou[ind(i, j)]/rho[ind(i, j)];
                v[ind(i, j)] = rhov[ind(i, j)]/rho[ind(i, j)];
            }
        }

        // Print some variable on console
        for(int j = ny-1; j >= 0; j--)
        {
            for(int i = 0; i < nx-1; i++)
                cout << u[ind(i, j)] << " ";

            cout << endl;
        }
        cout << endl;

        filename = "data/data_" + ZeroPadNumber(n) + ".csv";
        data.open(filename.c_str());
        for(int i = 0; i < nx; i++)
            for(int j = 0; j < ny; j++)
                data << x[i] << "\t" << y[j] << "\t" << rho[ind(i, j)] << "\t" << u[ind(i, j)] << "\t" << v[ind(i, j)] << "\n";
        data.close();
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
