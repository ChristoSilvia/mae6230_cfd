#include <iostream>
#include <math.h>
#include <ctime>
#include <stdio.h>
#include <fstream>

using namespace std;

int main() {
	// Mesh Parameters
	int nx = 256;
	int ny = 256;
	float H = 1.0;
	float dx = H/(nx - 1);
	float dy = H/(ny - 1);
	int nsteps = 3;
	printf("Steps in x direction: %d\n", nx);
	printf("Steps in y direction: %d\n", ny);
	printf("Space Step dx: %f\n", dx);
	printf("Space Step dy: %f\n", dy);
	
	// Specified Parameters
	float Ma = 0.1;
	float Re = 100.0;
	printf("Mach Number: %f\n", Ma);
	printf("Reynolds Number: %f\n", Re);
	
	// Simulation Parameters
	float rho_0 = 1.0;
	float c = 1.0;
	float c_squared = c*c;
	float U = Ma*c;
	float mu = rho_0*U*H/Re;
	printf("U: %f\n", U);
	printf("mu: %f\n", mu);

	// Timestep Parameters
	// stability is experimentally guarenteed when:
	// dt <= sigma/(1 + 2/Re_delta)*1/(abs(u)/dx + abs(v)/dy +
	//	c sqrt(1/dx^2 + 1/dy^2) )
	//
	// Re_delta = min(rho dx abs(u)/mu, rho dy abs(v)/ mu)
	//
	// We note that, abs(u) and abs(v) should be less than U
	//	in all cases.
	// Thus, we can use Re_delta = Re.
	float sigma = 0.5;
	float dt = sigma/((1.0+2.0/Re)*(U/dx + U/dy + c*sqrt(1/dx/dx+1/dy/dy)));
	printf("Timestep dt: %f\n", dt);
	
	// Product coeffs
	float a1 = dt/dx;
	float a2 = dt/dy;
	float a3 = mu*dt/dx/dx;
	float four_thirds_a3 = 4.0*a3/3.0;
	float a4 = mu*dt/dy/dy;
	float four_thirds_a4 = 4.0*a4/3.0;
	float a5 = mu*dt/(12.0*dx*dy);
	printf("a1: %f\n",a1);
	printf("a2: %f\n",a2);
	printf("a3: %f\n",a3);
	printf("a4: %f\n",a4);
	printf("a5: %f\n",a5);

	float b1 = 0.5*dt/dx;
	float b2 = 0.5*dt/dy;
	printf("b1: %f\n", b1);
	printf("b2: %f\n", b2);

	// time tracking
	clock_t begin_sim;
	clock_t end_sim;
	float sim_elapsed_time;
	clock_t begin_update;
	clock_t end_update;
	float elapsed_time;

	// Initialize arrays
	// a[i,j] = a[i*nx + j], 0 <= i < ny, 0 <= j < nx
	// Array to store predictor, and to store output of prediction
	float rho[nx*ny];
	float rho_u[nx*ny];
	float rho_v[nx*ny];
	float rho_p[nx*ny];
	float rho_u_p[nx*ny];
	float rho_v_p[nx*ny];
	float rho_q[nx*ny];
	float rho_u_q[nx*ny];
	float rho_v_q[nx*ny];
	float u[nx*ny];
	float v[nx*ny];
	
	// relative indices
	int h;
	int n;
	int s;
	int e;
	int w;
	int ne;
	int se;
	int sw;
	int nw;
	int nn;
	int ss;
	int ee;
	int ww;

	printf("Beginning Initialization\n");
	begin_update = clock();
	// Initialize the initial conditions
	for(int j=0; j<nx; j++){
		for(int i=0; i<ny-1; i++){
			h = i*nx + j;	
			rho[h] = rho_0;
			rho_p[h] = rho_0;
			rho_q[h] = rho_0;
			rho_u[h] = 0.0;
			rho_u_p[h] = 0.0;
			rho_u_q[h] = 0.0;
			rho_v[h] = 0.0;
			rho_v_p[h] = 0.0;
			rho_v_q[h] = 0.0;
		}
		h = (ny-1)*nx + j;	
		rho[h] = rho_0;
		rho_p[h] = rho_0;
		rho_q[h] = rho_0;
		rho_u[h] = U*rho_0;
		rho_u_p[h] = U*rho_0;
		rho_u_q[h] = U*rho_0;
		rho_v[h] = 0.0;
		rho_v_p[h] = 0.0;
		rho_v_q[h] = 0.0;
	}

	end_update = clock();
	elapsed_time = float(end_update-begin_update)/CLOCKS_PER_SEC;
	printf("Setup took: %f sec\n", elapsed_time);

	
	begin_sim = clock();
	for(int step=0; step<nsteps; step++){
		// Print step number
		/*
		if(step%100 == 0){
			begin_update = clock();
		}*/
		begin_update = clock();
	
		// PREDICTOR	

		// Calculate Velocities
		for(int i=0; i < ny; i++){
			for(int j=0; j<nx; j++) {
				u[i*nx+j] = rho_u[i*nx+j]/rho[i*nx+j];
				v[i*nx+j] = rho_v[i*nx+j]/rho[i*nx+j];
			}
		}

		// x=0 row
		for(int j=1; j < nx-1; j++) {
			h = j;
			n = h + nx;
			nn = n + nx;
			rho_p[j] = rho[j] - b2*(4*rho_v[n]-rho_v[nn]);
		}

		// interior
		for(int i=1; i < ny-1; i++) {
			// Update left side
			h = i*nx;
			e = h + 1;
			ee = h + 2;
			rho_p[h] = rho[h] - b1*(4*rho_u[e]-rho_u[ee]);
			
			// Loop over interior
			for(int j=1; j < nx-1; j++){
				h = i*nx+j;
				n = h + nx;
				e = h + 1;
				s = h - nx;
				w = h - 1;
				ne = n + 1;
				se = s + 1;
				sw = s - 1;
				nw = n - 1;
				// calculate predictors
				rho_p[h] = rho[h] - a1*(rho_u[e]-rho_u[h])-a2*(rho_v[n]-rho_v[h]);
				rho_u_p[h] = rho_u[h] - a1*(rho_u[e]*u[e] + c_squared*rho[e]-rho_u[h]*u[h]-c_squared*rho[h]) - a2*(rho_u[n]*v[n] - rho_u[h]*v[h]) + four_thirds_a3*(u[w]-2.0*u[h]+u[e]) + a4*(u[n]-2.0*u[h]+u[s]) + a5*(v[ne]+v[sw] - v[nw] - v[se]);
				rho_v_p[h] = rho_v[h] - a1*(rho_u[e]*v[e] - rho_u[h]*v[h]) - a2*(rho_v[n]*v[n]+c_squared*rho[n]-rho_v[h]*v[h]-c_squared*rho[h]) + a3*(v[w] - 2.0*v[h] + v[e]) + four_thirds_a4*(v[n] - 2.0*v[h] + v[s]) + a5*(u[ne] + u[sw] - u[nw] - u[se]);
			}
			// Update Right Side
			h = i*nx + nx-1;
			w = i*nx + nx - 2;
			ww = i*nx+nx-3;
			rho_p[h] = rho[h] - b1*(rho_u[ww] - 4*rho_u[w]);
		}

		// Top Left Corner
		h = (ny-1)*nx;
		e = h + 1;
		ee = h + 2;	
		rho_p[h] = rho[h] - b1*(-3*rho_u[h] + 4*rho_u[e] - rho_u[ee]);
		rho_u_p[h] = U*rho_p[h];
		// x=H row
		for(int j=1; j < nx-1; j++) {
			h = (ny-1)*nx + j;
			e = h + 1;
			w = h - 1;
			s = h - nx;
			ss = s - nx;
			rho_p[h] = rho[h] - b2*(rho_v[ss] - 4*rho_v[s]) - b1*(rho_u[e] - rho_u[w]);
			rho_u_p[h] = U*rho_p[h];
		}
		// Top Right Corner
		h = nx*ny-1;
		w = h - 1;
		ww = h - 2;
		rho_p[h] = rho[h] - b1*(rho_u[ww] - 4*rho_u[w] + 3*rho_u[h]);
		rho_u_p[h] = U*rho_p[h];
		
		// Corrector

		// Calculate Velocities
		for(int i=0; i < ny; i++){
			for(int j=0; j<nx; j++){
				u[i*nx+j] = rho_u_p[i*nx+j]/rho_p[i*nx+j];
				v[i*nx+j] = rho_v_p[i*nx+j]/rho_p[i*nx+j];
			}
		}
	
		// x=0 row
		for(int j=1; j < nx-1; j++) {
			h = j;
			n = j + nx;
			nn = n + nx;
			rho_q[h] = rho_p[h] - b2*(4*rho_v_p[n]-rho_v_p[nn]);
		}

		// interior
		for(int i=1; i < ny-1; i++) {
			// Update left side
			h = i*nx;
			e = h + 1;
			ee = h + 2;
			rho_q[h] = rho_p[h] - b1*(4*rho_u_p[e] - rho_u_p[ee]);
			
			// Loop over interior
			for(int j=1; j < nx-1; j++){
				h = i*nx + j;
				n = h + nx;
				e = h + 1;
				s = h - nx;
				w = h - 1;
				ne = n + 1;
				se = s + 1;
				sw = s - 1;
				nw = n - 1;

				// calculate predictors
				rho_q[h] = rho_p[h] - a1*(rho_u_p[h]-rho_u_p[w])-a2*(rho_v_p[h]-rho_v_p[s]);
				rho_u_q[h] = rho_u_p[h] - a1*(rho_u_p[h]*u[h] + c_squared*rho_p[h]-rho_u_p[w]*u[w]-c_squared*rho_p[w]) - a2*(rho_u_p[h]*v[h] - rho_u_p[s]*v[s]) + four_thirds_a3*(u[w]-2.0*u[h]+u[e]) + a4*(u[n]-2.0*u[h]+u[s]) + a5*(v[ne]+v[sw] - v[nw] - v[se]);
				rho_v_q[h] = rho_v_p[h] - a1*(rho_u_p[h]*v[h] - rho_u_p[w]*v[w]) - a2*(rho_v_p[h]*v[h]+c_squared*rho_p[h]-rho_v_p[s]*v[s]-c_squared*rho_p[s]) + a3*(v[w] - 2.0*v[h] + v[e]) + four_thirds_a4*(v[n] - 2.0*v[h] + v[s]) + a5*(u[ne] + u[sw] - u[nw] - u[se]);
			}
			
			// Update Right Side
			h = i*nx + nx - 1;
			w = h - 1;
			ww = h - 2;
			rho_q[h] = rho_p[h] - b1*(rho_u_p[ww] - 4*rho_u_p[w]);
		}
		
		// Top Left Cornder
		h = (ny-1)*nx;
		e = h + 1;
		ee = h + 2;
		rho_q[h] = rho_p[h] - b1*(-3*rho_u_p[h] + 4*rho_u_p[e] - rho_u_p[ee]);
		rho_u_q[h] = U*rho_q[h];
		// x=H row
		for(int j=1; j < nx-1; j++) {
			h = (ny-1)*nx + j;
			e = h + 1;
			w = h - 1;
			s = h - nx;
			ss = s - nx;
			rho_q[h] = rho_p[h] - b2*(rho_v_p[ss]-4*rho_v_p[s]) - b1*(rho_u_p[e] - rho_u_p[w]);
		}

		// Top Right Corner
		h = ny*nx-1;
		w = h - 1;
		ww = h - 2;
		rho_q[h] = rho_p[h] - b1*(rho_u_p[ww] - 4*rho_u_p[w] + 3*rho_u_p[h]);
		rho_u_q[h] = U*rho_q[h];
	
		// Averager
		// x=0 row
		for(int j=0; j < nx; j++) {
			rho[j] = 0.5*(rho_p[j] + rho_q[j]);
		}

		// interior
		for(int i=1; i < ny-1; i++) {
			// Initialize directions
			h = i*nx;
			// Update left side
			rho[h] = 0.5*(rho_p[h] + rho_q[h]);
			
			// Loop over interior
			for(int j=1; j < nx-1; j++){
				h = i*nx + j;
				rho[h] = 0.5*(rho_p[h] + rho_q[h]);
				rho_u[h] = 0.5*(rho_u_p[h] + rho_u_q[h]);
				rho_v[h] = 0.5*(rho_v_p[h] + rho_v_q[h]);
			}
			// h = i*nx+nx-1;
			h = i*nx + nx - 1;
			rho[h] = 0.5*(rho_p[h]+rho_q[h]);
		}
		
		// x=H row
		for(int j=0; j < nx; j++) {
			h = (ny-1)*nx + j;
			rho[h] = 0.5*(rho_p[h] + rho_q[h]);
			rho_u[h] = U*0.5*(rho_p[h] + rho_q[h]);
			rho_v[h] = 0.0;
		}
		/*
		// Write things to files
		ofstream rho_csv;
		ofstream rho_u_csv;
		ofstream rho_v_csv;
		ofstream rho_p_csv;
		ofstream rho_u_p_csv;
		ofstream rho_v_p_csv;
		ofstream rho_q_csv;
		ofstream rho_u_q_csv;
		ofstream rho_v_q_csv;
		char rho_buffer[16];
		char rho_u_buffer[20];
		char rho_v_buffer[20];
		char rho_p_buffer[20];
		char rho_u_p_buffer[24];
		char rho_v_p_buffer[24];
		char rho_q_buffer[20];
		char rho_u_q_buffer[24];
		char rho_v_q_buffer[24];
		sprintf(rho_v_buffer, "rho_v/rho_v_%04d.csv", step);
		sprintf(rho_buffer, "rho/rho_%04d.csv", step);
		sprintf(rho_u_p_buffer, "rho_u_p/rho_u_p_%04d.csv", step);
		sprintf(rho_v_p_buffer, "rho_v_p/rho_v_p_%04d.csv", step);
		sprintf(rho_p_buffer, "rho_p/rho_p_%04d.csv", step);
		sprintf(rho_u_q_buffer, "rho_u_q/rho_u_q_%04d.csv", step);
		sprintf(rho_v_q_buffer, "rho_v_q/rho_v_q_%04d.csv", step);
		sprintf(rho_q_buffer, "rho_q/rho_q_%04d.csv", step);
		sprintf(rho_u_buffer, "rho_u/rho_u_%04d.csv", step);
		rho_v_csv.open(rho_v_buffer);
		rho_csv.open(rho_buffer);
		rho_u_p_csv.open(rho_u_p_buffer);
		rho_v_p_csv.open(rho_v_p_buffer);
		rho_p_csv.open(rho_p_buffer);
		rho_u_q_csv.open(rho_u_q_buffer);
		rho_v_q_csv.open(rho_v_q_buffer);
		rho_q_csv.open(rho_q_buffer);
		rho_u_csv.open(rho_u_buffer);
		
	
		for(int i=0; i<ny; i++){
			for(int j=0; j<nx; j++){
				rho_csv << rho[i*nx+j] << ",";
				rho_u_csv << rho_u[i*nx+j] << ",";
				rho_v_csv << rho_v[i*nx+j] << ",";
				printf("printing\n");
				rho_csv << 0.5*(rho_p[i*nx+j] + rho_q[i*nx+j]) << ",";
				rho_u_csv << 0.5*(rho_u_p[i*nx+j] + rho_u_q[i*nx+j]) << ",";
				rho_v_csv << 0.5*(rho_v_p[i*nx+j] + rho_v_q[i*nx+j]) << ",";
				rho_p_csv << rho_p[i*nx+j] << ",";
				rho_u_p_csv << rho_u_p[i*nx+j] << ",";
				rho_v_p_csv << rho_v_p[i*nx+j] << ",";
				rho_q_csv << rho_q[i*nx+j] << ",";
				rho_u_q_csv << rho_u_q[i*nx+j] << ",";
				rho_v_q_csv << rho_v_q[i*nx+j] << ",";
			}
			rho_csv << "\n";
			rho_u_csv << "\n";
			rho_v_csv << "\n";
			rho_p_csv << "\n";
			rho_u_p_csv << "\n";
			rho_v_p_csv << "\n";
			rho_q_csv << "\n";
			rho_u_q_csv << "\n";
			rho_v_q_csv << "\n";
		}
		rho_csv.close();
		rho_u_csv.close();
		rho_v_csv.close();
		rho_p_csv.close();
		rho_u_p_csv.close();
		rho_v_p_csv.close();
		rho_q_csv.close();
		rho_u_q_csv.close();
		rho_v_q_csv.close();*/
		
		/*
		if(step%100==0){
			end_update = clock();
			elapsed_time = float(end_update - begin_update)/CLOCKS_PER_SEC;
			printf("On update: %d.       This Update Took: %f\n",step,elapsed_time);
		}*/
		end_update = clock();
		elapsed_time = float(end_update - begin_update)/CLOCKS_PER_SEC;
		printf("On update: %d.       This Update Took: %f\n",step,elapsed_time);
	}
	end_sim = clock();
	sim_elapsed_time = float(end_sim - begin_sim)/CLOCKS_PER_SEC;
	printf("Simulation took %f seconds.\n", sim_elapsed_time);

	// Write things to files
	ofstream x_csv;
	ofstream y_csv;
	ofstream rho_csv;
	ofstream rho_u_csv;
	ofstream rho_v_csv;
	ofstream u_csv;
	ofstream v_csv;
	ofstream p_csv;
	x_csv.open("x.csv");
	y_csv.open("y.csv");
	rho_csv.open("rho.csv");
	rho_u_csv.open("rho_u.csv");
	rho_v_csv.open("rho_v.csv");
	u_csv.open("u.csv");
	v_csv.open("v.csv");
	p_csv.open("p.csv");
	

	for(int i=0; i<ny; i++){
		for(int j=0; j<nx; j++){
			x_csv << j*dx << ",";
			y_csv << i*dy << ",";
			rho_csv << rho[i*nx+j] << ",";
			rho_u_csv << rho_u[i*nx+j] << ",";
			rho_v_csv << rho_v[i*nx+j] << ",";
			u_csv << rho_u[i*nx+j]/rho[i*nx+j] << ",";
			v_csv << rho_v[i*nx+j]/rho[i*nx+j] << ",";
			p_csv << c_squared*rho[i*nx+j] << ",";
		}
		x_csv << "\n";
		y_csv << "\n";
		rho_csv << "\n";
		rho_u_csv << "\n";
		rho_v_csv << "\n";
		u_csv << "\n";
		v_csv << "\n";
		p_csv << "\n";
	}
	x_csv.close();
	y_csv.close();
	rho_csv.close();
	rho_u_csv.close();
	rho_v_csv.close();
	u_csv.close();
	v_csv.close();
	p_csv.close();
}
