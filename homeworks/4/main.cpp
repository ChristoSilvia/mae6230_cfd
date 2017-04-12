#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int nx, ny, N, Nu, Nv;
float *p, *u, *up, *v, *vp;

/* FUNCTION: maps indices of 2D matrices to 1D
    arrays for matrices Aij with i = 0...n.
*/
int ind(int i, int j, int n)
{
    return i + j*(n+1);
}

/* Auxiliary discretization on the mesh */
float u_cen(int i, int j)// u at the center of the cell, u_{i, j}
{
    return 0.5*( u[ind(i+1,j,nx)] + u[ind(i,j,nx)] );
}

float u_lbc(int i, int j)// u at the left bottom corner of the cell, u_{i-1/2, j-1/2}
{
    return 0.5*( u[ind(i,j,nx)] + u[ind(i,j-1,nx)] );
}

float v_cen(int i, int j)// v at the center of the cell, v_{i, j}
{
    return 0.5*( v[ind(i,j+1,nx-1)] + v[ind(i,j,nx-1)] );
}

float v_lbc(int i, int j)// v at the left bottom corner of the cell, v_{i-1/2, j-1/2}
{
    return 0.5*( v[ind(i,j,nx-1)] + v[ind(i-1,j,nx-1)] );
}
/* End of Auxiliary Functions on the mesh */


/* FUNCTION: main */
int main()
{
    nx = 32;
    ny = 32;
    N = nx*ny;// # of pressure variables
    Nu = (nx+1)*ny;// # of u-vel variables
    Nv = nx*(ny+1);// # of v-vel variables

    int gs_tol_interval = 100; // Every 100 steps, check the error
    float l2err; // Float for keeping track of error in pressure poisson solution
    float l2tol = 1e-8; // If error is less than this, stop the iteration

    float Re = 100;
    float H = 1.0;
    float U = 1.0;
    float dx = H/nx;
    float dy = H/ny;
    
    float sigma = 0.1;
    float Re_d = Re*min(dx,dy);
    float dt = sigma/( 1.0 + 2.0/Re_d )/( U/dx + U/dy + sqrt( 1.0/dx/dx + 1.0/dy/dy ) );

    float a1 = dt/dx;
    float a2 = dt/dy;
    float a3 = dt/(Re*dx*dx);
    float a4 = dt/(Re*dy*dy);

    p = new float[N];// p_{i, j}
    u = new float[Nu];// u_{i-1/2, j}
    v = new float[Nv];// v_{i, j-1/2}
    up = new float[Nu];// u^\star_{i-1/2, j}
    vp = new float[Nv];// v^star_{i, j-1/2}


    /* Initial conditions for all the cells */
    for(int j=0; j<ny; j++)
    {// start from the bottom wall, from left to right.
        for(int i=0; i<=nx; i++)
        {
            u[ind(nx,j,nx)] = 0.0;
            up[ind(nx,j,nx)] = 0.0;
            if(i<nx)
            {
                v[ind(i,j,nx-1)] = 0.0;
                vp[ind(i,j,nx-1)] = 0.0;
                p[ind(i,j,nx-1)] = 0.0;
            }
        }
    }
    for(int i=0; i<nx; i++)// top wall
        v[ind(i,ny,nx-1)] = 0.0;
    /* End of Initial Conditions */


    /* Boundary conditions at the walls */
    float u_tw = U;// top
    float u_bw = 0.0;// bottom
    float v_lw = 0.0;// left
    float v_rw = 0.0;// right
    /* End of boundary conditions at the walls */


    int Nt = 10000;// # of iterations in the simulation


    /* Simulation begins */
    for(int n=0; n<=Nt; n++)
    {
        /* interior cells */
        for(int j=1; j<ny-1; j++)
        {// start from the bottom inner row, from left to right
            for(int i=1; i<=nx-1; i++)
            {
                up[ind(i,j,nx)] = u[ind(i,j,nx)]
                                    - a1*( u_cen(i,j)*u_cen(i,j) - u_cen(i-1,j)*u_cen(i-1,j) )
                                    - a2*( u_lbc(i,j+1)*v_lbc(i,j+1) - u_lbc(i,j)*v_lbc(i,j) )
                                    + a3*( u[ind(i+1,j,nx)] - 2.0*u[ind(i,j,nx)] + u[ind(i-1,j,nx)] )
                                    + a4*( u[ind(i,j+1,nx)] - 2.0*u[ind(i,j,nx)] + u[ind(i,j-1,nx)] );
                if(i<nx-1)
                    vp[ind(i,j,nx-1)] = v[ind(i,j,nx-1)]
                                        - a1*( v_lbc(i+1,j)*u_lbc(i+1,j) - v_lbc(i,j)*u_lbc(i,j) )
                                        - a2*( v_cen(i,j)*v_cen(i,j) - v_cen(i,j-1)*v_cen(i,j-1) )
                                        + a3*( v[ind(i+1,j,nx-1)] - 2.0*v[ind(i,j,nx-1)] + v[ind(i-1,j,nx-1)] )
                                        + a4*( v[ind(i,j+1,nx-1)] - 2.0*v[ind(i,j,nx-1)] + v[ind(i,j-1,nx-1)] );
            }
        }
        for(int i=1; i<nx-1; i++)// top inner row
        {
            int j = ny-1;
            vp[ind(i,j,nx-1)] = v[ind(i,j,nx-1)]
                                    - a1*( v_lbc(i+1,j)*u_lbc(i+1,j) - v_lbc(i,j)*u_lbc(i,j) )
                                    - a2*( v_cen(i,j)*v_cen(i,j) - v_cen(i,j-1)*v_cen(i,j-1) )
                                    + a3*( v[ind(i+1,j,nx-1)] - 2.0*v[ind(i,j,nx-1)] + v[ind(i-1,j,nx-1)] )
                                    + a4*( v[ind(i,j+1,nx-1)] - 2.0*v[ind(i,j,nx-1)] + v[ind(i,j-1,nx-1)] );
        }
        /* End of interior cells */

        /* Left and right walls */
        for(int j=1; j<ny; j++)
        {
            int i = 0;// left
            vp[ind(i,j,nx-1)] = v[ind(i,j,nx-1)]
                                    - a1*( v_lbc(i+1,j)*u_lbc(i+1,j) - v_lbc(i,j)*u_lbc(i,j) )
                                    - a2*( v_cen(i,j)*v_cen(i,j) - v_cen(i,j-1)*v_cen(i,j-1) )
                                    + a3*( v[ind(1,j,nx-1)] - 3.0*v[ind(0,j,nx-1)] + 2.0*v_lw )// Term that is different from inner cells
                                    + a4*( v[ind(i,j+1,nx-1)] - 2.0*v[ind(i,j,nx-1)] + v[ind(i,j-1,nx-1)] );
            i = nx-1;// right
            vp[ind(i,j,nx-1)] = v[ind(i,j,nx-1)]
                                    - a1*( v_lbc(i+1,j)*u_lbc(i+1,j) - v_lbc(i,j)*u_lbc(i,j) )
                                    - a2*( v_cen(i,j)*v_cen(i,j) - v_cen(i,j-1)*v_cen(i,j-1) )
                                    + a3*( 2.0*v_rw - 3.0*v[ind(nx-1,j,nx-1)] + v[ind(nx-2,j,nx-1)] )// Term that is different from inner cells
                                    + a4*( v[ind(i,j+1,nx-1)] - 2.0*v[ind(i,j,nx-1)] + v[ind(i,j-1,nx-1)] );
        }

        /* Top and bottom walls */
        for(int i = 1; i<nx; i++)
        {
            int j = 0;//bottom
            up[ind(i,j,nx)] = u[ind(i,j,nx)]
                                - a1*( u_cen(i,j)*u_cen(i,j) - u_cen(i-1,j)*u_cen(i-1,j) )
                                - a2*( u_lbc(i,j+1)*v_lbc(i,j+1) - u_lbc(i,j)*v_lbc(i,j) )
                                + a3*( u[ind(i+1,j,nx)] - 2.0*u[ind(i,j,nx)] + u[ind(i-1,j,nx)] )
                                + a4*( u[ind(i,1,nx)] - 3.0*u[ind(i,0,nx)] + 2.0*u_bw );// Term that is different from inner cells
            j = ny-1;// top
            up[ind(i,j,nx)] = u[ind(i,j,nx)]
                                - a1*( u_cen(i,j)*u_cen(i,j) - u_cen(i-1,j)*u_cen(i-1,j) )
                                - a2*( u_lbc(i,j+1)*v_lbc(i,j+1) - u_lbc(i,j)*v_lbc(i,j) )
                                + a3*( u[ind(i+1,j,nx)] - 2.0*u[ind(i,j,nx)] + u[ind(i-1,j,nx)] )
                                + a4*( 2.0*u_tw - 3.0*u[ind(i,ny-1,nx)] + u[ind(i,ny-2,nx)] );// Term that is different from inner cells
        }
        // Page 7 on assignment: the RHSs of Eq. (16) should be zero.
        float volume_integral = 0.0;// sum of all RHSs of the Poisson Equation.
        for(int j=0; j<ny; j++)
            for(int i=0; i<nx; i++)
                volume_integral += (1/dx)*( up[ind(i+1,j,nx)] - up[ind(i,j,nx)] ) + (1/dy)*( vp[ind(i,j+1,nx-1)] - vp[ind(i,j,nx-1)] );

        /* Here we make the source term of p */


        /* Here we do our Gauss-Seidel Iterations */

        // loop control variables        
        bool not_over = true;
        int nits = 0;
        float error = 0.0;
        float last_error = 100.0;

        // common fractions
        float half = 0.5;
        float third = 1.0/3.0;
        float fourth = 0.25;

        // ratios
        float xfactor = dx/dt;
        float yfactor = dy/dx;
        float invdeltasq = 1.0/dx/dx;
        float xerrfactor = 1.0/dx/dt;
        float yerrfactor = 1.0/dy/dt;

        // indices
        int i,j;
        
        cout << "writing source term\n";
        ofstream sourcef;
        char source_buffer_fname[23];
        sprintf(source_buffer_fname, "source/source_%05d.csv", n);
        cout << source_buffer_fname << "\n";
        sourcef.open(source_buffer_fname);
        // write source term
        for(j=0; j<ny; j++){
            for(i=0; i<nx; i++){
                sourcef << xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                         + yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]) << "\n";
            }
        }
        sourcef.close();
        

        while(not_over) {
            // lower left corner
            i=0; j=0;
            p[ind(i,j,nx)] = half*(p[ind(i+1,j,nx-1)] 
                                 + p[ind(i,j+1,nx-1)] 
                                 - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                 - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
            // lower boundary
            j=0;
            for(i=1; i<nx-1; i++){
                p[ind(i,j,nx)] = third*(p[ind(i-1,j,nx-1)] 
				      + p[ind(i+1,j,nx-1)] 
                                      + p[ind(i,j+1,nx-1)] 
                                      - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                      - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
            }
            // lower right corner
            i=nx-1; j=0;
            p[ind(i,j,nx)] = half*(p[ind(i-1,j,nx-1)] 
                                 + p[ind(i,j+1,nx-1)] 
                                 - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                 - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
            // main loop
            for(j=1; j<ny-1; j++){
                // left wall
                i=0;
                p[ind(i,j,nx)] = third*(p[ind(i,j-1,nx-1)] 
				      + p[ind(i+1,j,nx-1)] 
                                      + p[ind(i,j+1,nx-1)] 
                                      - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                      - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                // interior
                for(i=1; i<nx-1; i++){
                    p[ind(i,j,nx)] = fourth*(p[ind(i,j-1,nx-1)] 
	            		          + p[ind(i-1,j,nx-1)] 
	            		          + p[ind(i+1,j,nx-1)] 
                                          + p[ind(i,j+1,nx-1)] 
                                          - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                          - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                }
                // right wall
                i=nx-1;
                p[ind(i,j,nx)] = third*(p[ind(i,j-1,nx-1)] 
				      + p[ind(i-1,j,nx-1)] 
                                      + p[ind(i,j+1,nx-1)] 
                                      - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                      - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
            }
            // upper left corner
            i=0; j=ny-1;
            p[ind(i,j,nx)] = half*(p[ind(i+1,j,nx-1)] 
                                 + p[ind(i,j-1,nx-1)] 
                                 - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                 - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
            // upper boundary
            j=ny-1;
            for(i=1; i<nx-1; i++){
                p[ind(i,j,nx)] = third*(p[ind(i,j-1,nx-1)] 
				      + p[ind(i-1,j,nx-1)] 
                                      + p[ind(i+1,j,nx-1)] 
                                      - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                      - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
            }
            // upper right corner
            i=nx-1; j=ny-1;
            p[ind(i,j,nx)] = half*(p[ind(i-1,j,nx-1)] 
                                 + p[ind(i,j-1,nx-1)] 
                                 - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                 - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
            /*
            error = 0.0;
            if(nits % 100 == 0){
                // calculate l-infinity error
                // lower left corner
                i=0; j=0;
                error = max(abs( p[ind(i+1,j,nx)] + p[ind(i,j+1,nx)] - 2*p[ind(i,j,nx)]
half*(p[ind(i+1,j,nx)] 
                                     + p[ind(i,j+1,nx)] 
                                     - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                     - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                // lower boundary
                j=0;
                for(i=1; i<nx-1; i++){
                    p[ind(i,j,nx)] = third*(p[ind(i-1,j,nx)] 
            			      + p[ind(i+1,j,nx)] 
                                          + p[ind(i,j+1,nx)] 
                                          - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                          - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                }
                // lower right corner
                i=nx-1; j=0;
                p[ind(i,j,nx)] = half*(p[ind(i-1,j,nx)] 
                                     + p[ind(i,j+1,nx)] 
                                     - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                     - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                // main loop
                for(j=1; j<ny-1; j++){
                    // left wall
                    i=0;
                    p[ind(i,j,nx)] = third*(p[ind(i,j-1,nx)] 
            			      + p[ind(i+1,j,nx)] 
                                          + p[ind(i,j+1,nx)] 
                                          - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                          - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                    // interior
                    for(i=1; i<nx-1; i++){
                        p[ind(i,j,nx)] = fourth*(p[ind(i,j-1,nx)] 
                        		          + p[ind(i-1,j,nx)] 
                        		          + p[ind(i+1,j,nx)] 
                                              + p[ind(i,j+1,nx)] 
                                              - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                              - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                    }
                    // right wall
                    i=nx-1;
                    p[ind(i,j,nx)] = third*(p[ind(i,j-1,nx)] 
            			      + p[ind(i-1,j,nx)] 
                                          + p[ind(i,j+1,nx)] 
                                          - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                          - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                }
                // upper left corner
                i=0; j=ny-1;
                p[ind(i,j,nx)] = half*(p[ind(i+1,j,nx)] 
                                     + p[ind(i,j-1,nx)] 
                                     - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                     - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                // upper boundary
                j=ny-1;
                for(i=1; i<nx-1; i++){
                    p[ind(i,j,nx)] = third*(p[ind(i,j-1,nx)] 
            			      + p[ind(i-1,j,nx)] 
                                          + p[ind(i+1,j,nx)] 
                                          - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                          - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                }
                // upper right corner
                i=nx-1; j=ny-1;
                p[ind(i,j,nx)] = half*(p[ind(i-1,j,nx)] 
                                     + p[ind(i,j-1,nx)] 
                                     - xfactor*(u[ind(i+1,j,nx)]-u[ind(i,j,nx)])
                                     - yfactor*(v[ind(i,j+1,nx-1)]-v[ind(i,j,nx-1)]));
                */
            if(true){
                cout << "writing file\n";
                ofstream pf;
                char p_buffer_fname[22];
                sprintf(p_buffer_fname, "press%05d/p_%05d.csv", n, nits);
                cout << p_buffer_fname << "\n";
                pf.open(p_buffer_fname);
                // write p
                for(int k=0; k<N; k++)
                    pf << p[k] << "\n";
                pf.close();
            }
            nits++;
            if(nits > 1000)
                not_over=false;
        }


                


        /* OK, we have the new p so we can compute the new velocities. */

        for(int j=0; j<ny-1; j++)
        {
            for(int i=0; i<=nx-1; i++)
            {
                v[ind(i,j+1,nx-1)] = vp[ind(i,j+1,nx-1)] - a2*( p[ind(i,j+1,nx-1)] - p[ind(i,j,nx-1)] );
                if(i<nx-1)
                    u[ind(i+1,j,nx)] = up[ind(i+1,j,nx)] - a1*( p[ind(i+1,j,nx-1)] - p[ind(i,j,nx-1)] );
            }
        }
        for(int i=1; i<nx; i++)
            u[ind(i,ny-1,nx)] = up[ind(i,ny-1,nx)] - a1*( p[ind(i,ny-1,nx-1)] - p[ind(i-1,ny-1,nx-1)] );


        if(true){
            ofstream uf;
            ofstream vf;
            ofstream pf;
            char u_buffer_fname[15];
            char v_buffer_fname[15];
            char p_buffer_fname[15];
            sprintf(u_buffer_fname, "run/u_%05d.csv", n);
            sprintf(v_buffer_fname, "run/v_%05d.csv", n);
            sprintf(p_buffer_fname, "run/p_%05d.csv", n);
            uf.open(u_buffer_fname);
            vf.open(v_buffer_fname);
            pf.open(p_buffer_fname);
            // write u
            for(int k=0; k<Nu; k++)
                uf << u[k] << "\n";
            // write v
            for(int k=0; k<Nv; k++)
                vf << v[k] << "\n";
            // write p
            for(int k=0; k<N; k++)
                pf << p[k] << "\n";
            uf.close();
            vf.close();
            pf.close();
        }

    cout << "\t Iteration " << n << " done. Volume integral :" << volume_integral << "\n\n";
    }/* Simulation ends */


    /* Release allocated memory */
    delete [] p;
    delete [] u;
    delete [] v;
    delete [] up;
    delete [] vp;

    /* exit program */
    return 0;

}
