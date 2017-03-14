#include <iostream>
#include <math.h>
#include <ctime>
#include <stdio.h>
#include <fstream>

using namespace std;

// Global Mesh Variables
int nx;
int ny;

int ind(int i, int j) {
	return nx*j + i;
}

int main(int argc, char* argv[]) {
	// Number of Timesteps
	int nsteps = 400000;

	// Mesh Parameters
	nx = 256;
	ny = 256;
	int N = nx*ny;
	float H = 1.0;
	float dx = H/(nx-1);
	float dy = H/(ny-1);

	// Specified Parameters
	float Ma = 0.1;
	float Re = 100.0;

	// Dimensional Parameters
	float rho_0 = 1.0;
	float c = 1.0;
	float c_squared = c*c;
	float U = Ma*c;
	float mu = rho_0*U*H/Re;

	// Time Parameters
	float sigma = 0.3;
	float Re_delta = rho_0*U/mu*min(dx, dy);
	float dt = sigma/((1.0 + 2.0/Re_delta)*(U/dx + U/dy + c*sqrt(1/dx/dx + 1/dy/dy)));

	// Product Coefficients
	float a1 = dt/dx;
	float a2 = dt/dy;
	float a3 = mu*dt/dx/dx;
	float a343 = (4.0/3.0)*a3;
	float a4 = mu*dt/dy/dy;
	float a443 = (4.0/3.0)*a4;
	float a5 = mu*dt/(12.0*dx*dy);
	float b1 = 0.5*dt/dx;
	float b2 = 0.5*dt/dy;

	// Initializing Arrays
	float rho[N];
	float rho_u[N];
	float rho_v[N];
	float rho_p[N];
	float rho_u_p[N];
	float rho_v_p[N];
	float u[N];
	float v[N];

	// Giving Arrays Initial Values
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny-1; j++){
			rho[ind(i,j)] = rho_0;
			rho_u[ind(i,j)] = 0.0;
			rho_v[ind(i,j)] = 0.0;
			rho_p[ind(i,j)] = rho_0;
			rho_u_p[ind(i,j)] = 0.0;
			rho_v_p[ind(i,j)] = 0.0;
		}
		rho[ind(i,ny-1)] = rho_0;
		rho_u[ind(i,ny-1)] = U*rho_0;
		rho_v[ind(i,ny-1)] = 0.0;
		rho_p[ind(i,ny-1)] = rho_0;
		rho_u_p[ind(i,ny-1)] = U*rho_0;
		rho_v_p[ind(i,ny-1)] = 0.0;
	}

	for(int step_index=0; step_index<nsteps; step_index++){
		
		// Calculate Velocities
		for(int j=0; j<ny; j++){
			for(int i=0; i<nx; i++){
				u[ind(i,j)] = rho_u[ind(i,j)]/rho[ind(i,j)];	
				v[ind(i,j)] = rho_v[ind(i,j)]/rho[ind(i,j)];	
			}
		}


		// Bottom Boundary
		for(int i=1; i<nx-1; i++){
			// Bottom Boundary
			rho_p[ind(i,0)] = rho[ind(i,0)]
					- b2*(	-3*rho_v[ind(i,0)]
						+4*rho_v[ind(i,1)]
						-1*rho_v[ind(i,2)]);
		}	
	
		// Interior First
		for(int j=1; j<ny-1; j++){
			// Left Boundary
			rho_p[ind(0,j)] = rho[ind(0,j)] 
					- b1*(	-3*rho_u[ind(0,j)]
						+4*rho_u[ind(1,j)]
						-1*rho_u[ind(2,j)]);

			for(int i=1; i<nx-1; i++){

				// Mass equation
				rho_p[ind(i,j)] = rho[ind(i,j)]

				- a1*(rho_u[ind(i+1,j)] - rho_u[ind(i,j)])

				- a2*(rho_v[ind(i,j+1)] - rho_v[ind(i,j)]);

				// x momentum equation
				rho_u_p[ind(i,j)] = rho_u[ind(i,j)]

				- a1*(	rho_u[ind(i+1,j)]*u[ind(i+1,j)]
					+ c_squared*rho[ind(i+1,j)]
					- rho_u[ind(i,j)]*u[ind(i,j)]
					- c_squared*rho[ind(i,j)])

				- a2*(	rho_u[ind(i,j+1)]*v[ind(i,j+1)]
					- rho_u[ind(i,j)]*v[ind(i,j)])

				+ a343*(u[ind(i+1,j)]
					- 2*u[ind(i,j)]
					+ u[ind(i-1,j)])

				+ a4*(	u[ind(i,j+1)]
					- 2*u[ind(i,j)]
					+ u[ind(i,j-1)] )

				+ a5*(	v[ind(i+1,j+1)]
					+ v[ind(i-1,j-1)]
					- v[ind(i+1,j-1)]
					- v[ind(i-1,j+1)] );

				// y momentum equation
				rho_v_p[ind(i,j)] = rho_v[ind(i,j)]

				- a1*(	rho_v[ind(i+1,j)]*u[ind(i+1,j)]
					- rho_v[ind(i,j)]*u[ind(i,j)])

				- a2*(	rho_v[ind(i,j+1)]*v[ind(i,j+1)]
					+ c_squared*rho[ind(i,j+1)]
					- rho_v[ind(i,j)]*v[ind(i,j)]
					- c_squared*rho[ind(i,j)])

				+ a3*(	v[ind(i+1,j)]
					- 2*v[ind(i,j)]
					+ v[ind(i-1,j)])

				+ a443*(v[ind(i,j+1)]
					- 2*v[ind(i,j)]
					+ v[ind(i,j-1)] )

				+ a5*(	u[ind(i+1,j+1)]
					+ u[ind(i-1,j-1)]
					- u[ind(i+1,j-1)]
					- u[ind(i-1,j+1)] );
			}
			// Right Boundary
			rho_p[ind(nx-1,j)] = rho[ind(nx-1,j)]
					+ b1*(	-3*rho_u[ind(nx-1,j)]
						+4*rho_u[ind(nx-2,j)]
						-1*rho_u[ind(nx-3,j)]);
		}


		// Top Boundary
		for(int i=1; i<nx-1; i++){
			rho_p[ind(i,ny-1)] = rho[ind(i,ny-1)]
					+ b2*(	-3*rho_v[ind(i,ny-1)]
						+4*rho_v[ind(i,ny-2)]
						-1*rho_v[ind(i,ny-3)])
					- b1*U*(rho[ind(i+1,ny-1)]
						-rho[ind(i-1,ny-1)]);
			rho_u_p[ind(i,ny-1)] = U*rho_p[ind(i,ny-1)];
		}
	
		//////////////////////////////////////////////////////
		// CORRECTOR	
		//////////////////////////////////////////////////////
		// Calculate Velocities
		for(int j=0; j<ny; j++){
			for(int i=0; i<nx; i++){
				u[ind(i,j)] = rho_u_p[ind(i,j)]/rho_p[ind(i,j)];	
				v[ind(i,j)] = rho_v_p[ind(i,j)]/rho_p[ind(i,j)];	
			}
		}
		
		for(int i=1; i<nx-1; i++){
			// Bottom Boundary
			rho[ind(i,0)] = 0.5*(rho[ind(i,0)] + rho_p[ind(i,0)]
					- b2*(	-3*rho_v_p[ind(i,0)]
						+4*rho_v_p[ind(i,1)]
						-1*rho_v_p[ind(i,2)]));
		}
	
		// Interior First
		for(int j=1; j<ny-1; j++){
			// Left Boundary
			rho[ind(0,j)] = 0.5*(rho[ind(0,j)] + rho_p[ind(0,j)]
					- b1*(	-3*rho_u_p[ind(0,j)]
						+4*rho_u_p[ind(1,j)]
						-1*rho_u_p[ind(2,j)]));
			for(int i=1; i<nx-1; i++){

				// Mass Equation
				rho[ind(i,j)] = 0.5*(rho[ind(i,j)]
					+ rho_p[ind(i,j)] 

				- a1*(rho_u_p[ind(i,j)] - rho_u_p[ind(i-1,j)])

				- a2*(rho_v_p[ind(i,j)] - rho_v_p[ind(i,j-1)]));

				// x momentum equation
				rho_u[ind(i,j)] = 0.5*(rho_u[ind(i,j)]
					+ rho_u_p[ind(i,j)]

				- a1*(	rho_u_p[ind(i,j)]*u[ind(i,j)]
					+ c_squared*rho_p[ind(i,j)]
					- rho_u_p[ind(i-1,j)]*u[ind(i-1,j)]
					- c_squared*rho_p[ind(i-1,j)])

				- a2*(	rho_u_p[ind(i,j)]*v[ind(i,j)]
					- rho_u_p[ind(i,j-1)]*v[ind(i,j-1)])

				+ a343*(u[ind(i+1,j)]
					- 2*u[ind(i,j)]
					+ u[ind(i-1,j)])

				+ a4*(	u[ind(i,j+1)]
					- 2*u[ind(i,j)]
					+ u[ind(i,j-1)] )

				+ a5*(	v[ind(i+1,j+1)]
					+ v[ind(i-1,j-1)]
					- v[ind(i+1,j-1)]
					- v[ind(i-1,j+1)] ));

				// y momentum equation
				rho_v[ind(i,j)] = 0.5*(rho_v[ind(i,j)]
					+ rho_v_p[ind(i,j)]

				- a1*(	rho_v_p[ind(i,j)]*u[ind(i,j)]
					- rho_v_p[ind(i-1,j)]*u[ind(i-1,j)])

				- a2*(	rho_v_p[ind(i,j)]*v[ind(i,j)]
					+ c_squared*rho_p[ind(i,j)]
					- rho_v_p[ind(i,j-1)]*v[ind(i,j-1)]
					- c_squared*rho_p[ind(i,j-1)])

				+ a3*(	v[ind(i+1,j)]
					- 2*v[ind(i,j)]
					+ v[ind(i-1,j)])

				+ a443*(v[ind(i,j+1)]
					- 2*v[ind(i,j)]
					+ v[ind(i,j-1)] )

				+ a5*(	u[ind(i+1,j+1)]
					+ u[ind(i-1,j-1)]
					- u[ind(i+1,j-1)]
					- u[ind(i-1,j+1)] ));
			}
			// Right Boundary
			rho[ind(nx-1,j)] = 0.5*(rho[ind(nx-1,j)] + rho_p[ind(nx-1,j)]
					+ b1*(	-3*rho_u_p[ind(nx-1,j)]
						+4*rho_u_p[ind(nx-2,j)]
						-1*rho_u_p[ind(nx-3,j)]));
		}

		// Top Boundary
		for(int i=1; i<nx-1; i++){
			rho[ind(i,ny-1)] = 0.5*(rho[ind(i,ny-1)] + rho_p[ind(i,ny-1)]
					+ b2*(	-3*rho_v[ind(i,ny-1)]
						+4*rho_v[ind(i,ny-2)]
						-1*rho_v[ind(i,ny-3)])
					- b1*U*(rho[ind(i+1,ny-1)]
						-rho[ind(i-1,ny-1)]));
			rho_u[ind(i,ny-1)] = U*rho[ind(i,ny-1)];
		}

		/*
		// Printer
		ofstream rho_csv;
		ofstream rho_u_csv;
		ofstream rho_v_csv;
		char rho_buffer[16];
		char rho_u_buffer[20];
		char rho_v_buffer[20];
		sprintf(rho_buffer, "rho/rho_%05d.csv", step_index);
		sprintf(rho_u_buffer, "rho_u/rho_u_%05d.csv", step_index);
		sprintf(rho_v_buffer, "rho_v/rho_v_%05d.csv", step_index);
		rho_csv.open(rho_buffer);
		rho_u_csv.open(rho_u_buffer);
		rho_v_csv.open(rho_u_buffer);
		
		for(int j=0; j<ny; j++){
			for(int i=0; i<nx-1; i++){
				rho_csv << rho[ind(i,j)] << ",";
				rho_u_csv << rho_u[ind(i,j)] << ",";
				rho_v_csv << rho_v[ind(i,j)] << ",";
			}
			rho_csv << rho[ind(nx-1,j)] << "\n";
			rho_u_csv << rho_u[ind(nx-1,j)] << "\n";
			rho_v_csv << rho_v[ind(nx-1,j)] << "\n";
		}
		rho_csv.close();
		rho_u_csv.close();
		rho_v_csv.close();
		cout << "Completed Iteration: " << step_index << "\n";
		*/
	}

	
	cout << "Completed Simulation\n";

	ofstream rho_csv;
	ofstream rho_u_csv;
	ofstream rho_v_csv;
	rho_csv.open("rho.csv");
	rho_u_csv.open("rho_u.csv");
	rho_v_csv.open("rho_v.csv");
	
	for(int j=0; j<ny; j++){
		for(int i=0; i<nx-1; i++){
			rho_csv << rho[ind(i,j)] << ",";
			rho_u_csv << rho_u[ind(i,j)] << ",";
			rho_v_csv << rho_v[ind(i,j)] << ",";
		}
		rho_csv << rho[ind(nx-1,j)] << "\n";
		rho_u_csv << rho_u[ind(nx-1,j)] << "\n";
		rho_v_csv << rho_v[ind(nx-1,j)] << "\n";
	}
	rho_csv.close();
	rho_u_csv.close();
	rho_v_csv.close();
}	
