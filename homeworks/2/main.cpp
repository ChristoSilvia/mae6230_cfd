#include <iostream>
#include <math.h>

using namespace std;
const double PI = 4*atan(1.);

/*  FUNCTION: 1st derivative of f(x).
    Given a stencil, it generates the finite diff coefficients

    x = point at which dfdx is evaluated
    ns = number of points on the stencil
    xs = array of locations of stencil points
    D = array to save coefficients
*/
void dfdx(double x, int ns, double *xs, double *D){
    // For the k-th stencil point, calculate its 
    //     coefficient and store it in D[k].
    double aux;
    for(int i=0; i<ns; i++){
        D[i] = 0;
        for(int k=0; k<ns; k++){
            aux = 1/(xs[i]-xs[k]);
            for(int j=0; j<ns; m++){
                if(j!=i && j!=k)
                    aux *= (x - xs[j])/(xs[i] - xs[j]);
            }
            D[k] += aux;
        }
    }
}

/*  FUNCTION: Mesh generator.

    n = number of cells on the mesh
    x = array of mesh locations
*/
void mesh(int n, double *x)
{
    for(int i=0; i<=n; i++)
        x[i] = (double)i/n + sin(2*PI*i/n)/10.;
}

/* FUNCTION: analytical test function.
*/
double y(double x)
{
    return exp(exp(x));
}

/* FUNCTION: First derivative of y(x).
*/
double dydx(double x)
{
    return exp(x)*exp(exp(x));
}

int main()
{
    // Generate a mesh of n+1 locations
    int n = 30;
    double *x = new double[n+1];
    mesh(n, x);
    
    // Generate the coefficients at each point on the mesh using a stencil of ns locations.
    // End points without enough info to calculate the derivative must be excluded.
    int ns = 5;
    for(int i=ns/2; i<=n-ns/2; i++)// Limits valid for a centered (ns odd) stencil.
    {
        // Extract the stencil points and save them in temp array xs.
        double* xs = new double[ns];
        for(int l=0; l<ns; l++)
            xs[l] = x[i-ns/2+l];
        
        // Compute the ns derivative coeffs at x[i] and save them in temp array D.
        double *D = new double[ns];
        dfdx(x[i], ns, xs, D);

        // With the coefficients obtained, compute the approximation for f(x_i) (Eq. 6)
        double dfnumdx = 0;
        for(int k=0; k<ns; k++)
            dfnumdx = dfnumdx + y(xs[k])*D[k];

        // Compare with the analytical expression and print it.
        double error = dydx(x[i]) - dfnumdx;
        cout << error << endl;

        // Once the stencil is no longer needed, release memory
        delete[] xs;
        delete[] D;
    }

    // Release memory storing the mesh...
    delete[] x;

    // ...and exit program.
    return 0;

}
