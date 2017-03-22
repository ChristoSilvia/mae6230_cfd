#include <iostream>

using namespace std;

int nx, ny, N;

int ind(int i, int j){
	return i+nx*j;
}

int main(){
	int nits = 4000;
	int sample_interval = 1;

	// good multigrid numbers are powers of 2 plus 1
	nx = 17;
	ny = 17;
	N = nx*ny;
	float delta = 1.0/(nx-1);
	float delta_squared = delta*delta;

	float *source = new float[N];
	float *source_times_delta_squared = new float[N];
	float *u = new float[N];
	float *usol = new float[N];

	float x,y;
	for(int j =0; j < ny; j++){
		for(int i = 0; i < nx; i++){
			x = i*delta;
			y = j*delta;
			source[ind(i,j)] = -2*((1-6*x*x)*y*y*(1-y*y) + (1-6*y*y)*x*x*(1-x*x));
			source_times_delta_squared[ind(i,j)] = delta_squared*source[ind(i,j)];
			usol[ind(i,j)] = (1-x*x)*x*x*(y*y-1)*y*y;
			u[ind(i,j)] = 0.0;
		}
	}

	float l2_err;
	float laplacian_times_delta_squared_here;
	float error_here;
	for(int it=0; it<nits; it++){
		for(int j=1; j<ny-1; j++){
			for(int i=1; i<nx-1;i++){
				u[ind(i,j)] = 0.25*(u[ind(i,j-1)]
						+ u[ind(i-1,j)]
						+ u[ind(i+1,j)]
						+ u[ind(i,j+1)]
						- source_times_delta_squared[ind(i,j)]);
			}
		}
		if(it % sample_interval == 0){
			l2_err=0;
			for(int j=1; j<ny-1; j++){
				for(int i=1; i<nx-1; i++){
					error_here = u[ind(i,j)]-usol[ind(i,j)];
					l2_err += error_here*error_here;
				}
			}
			cout << l2_err << "\n";
		}
	}
}
