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
void dfdx(double x, int ns, double *xs, double *D)
{
    /* PRINT xs
    for(int k=0; k<ns; k++)
        cout << xs[k] << " ";
    cout << endl; */

    // Since the finite diff coefficients depend on the differences (xs[l]-xs[m]),
    // we generate an auxiliar matrix "dx" containing these differences.
    double** dxs = new double*[ns];
    for(int l=0; l<ns; l++)
    {
        dxs[l] = new double[ns];
        for(int m=0; m<ns; m++)
            dxs[l][m] = xs[l] - xs[m];
    }
    /* PRINT dxs
    for(int l=0; l<ns; l++)
    {
        for(int m=0; m<ns; m++)
            cout << dxs[l][m] << " ";
        cout << endl;
    }
    cout << endl; */


    // For the k-th stencil point, calculate its coefficient and store it in D[k].

    double aux;// auxiliar variable to store the inner product in dfk/dx (Eq. 4)
    for(int k=0; k<ns; k++)
    {
        D[k] = 0;
        for(int l=0; l<ns; l++)
            if(l!=k)
            {
                aux = 1;
                for(int m=0; m<ns; m++)
                {
                    if(m!=k && m!=l)
                        aux = aux*(x - xs[m])/dxs[k][m];
                }
                aux = aux/dxs[k][l];
                D[k] = D[k] + aux;
            }
    }

    /* PRINT D[k]
    for(int k=0; k<ns; k++)
        cout << D[k] << " ";
    cout << endl; */


    // Release memory allocated for dxs
    for(int l=0; l<ns; l++)
        delete[] dxs[l];
    delete[] dxs;
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
    /* PRINT x
    for(int i=0; i<=n; i++)
        cout << x[i] << endl;
    cout << endl; */


    // Generate the coefficients at each point on the mesh using a stencil of ns locations.
    // End points without enough info to calculate the derivative must be excluded.
    int ns = 5;
    for(int i=ns/2; i<=n-ns/2; i++)// Limits valid for a centered (ns odd) stencil.
    {
        // Extract the stencil points and save them in temp array xs.
        double* xs = new double[ns];
        for(int l=0; l<ns; l++)
            xs[l] = x[i-ns/2+l];
        // Passing a new matrix makes notation in dfdx() coincide with that on the pdf.

        /* PRINT stencil locations (this should coincide with xs in dfdx() )
        for(int l=0; l<ns; l++)
            cout << xs[l] << " ";
        cout << endl; */

        // Compute the ns derivative coeffs at x[i] and save them in temp array D.
        double *D = new double[ns];
        dfdx(x[i], ns, xs, D);

        /* Print finite diff coeffs (this should coincide with D in dfdx() )
        for(int k=0; k<ns; k++)
            cout << D[k] << " ";
        cout << endl; */


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
