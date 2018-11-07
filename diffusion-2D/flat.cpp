#include <iostream>
#include <fstream>
#include <string>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>


using namespace std;

int main()
{
	const double H = 0.4;	/// total height
	const double L = 0.3; /// total length
	const double W = 0.01;	/// total thickness
	const double k = 1000;	 /// thermal conductivity

	const double q = 5e35;	/// thermal source at the west boundary, W/m^2
	const double h = 253.165;	/// convective heat transfer coefficient with the south boundary, W/(m^2*K)
	const double T_s = 200;	/// const temperature near the south boundary, Celsius
	const double T_n = 100;	/// const temperature at the north boundary, Celsius
	
	/// if NX = 30 and NY = 40, or NX and NY is bigger than that, there is Segmentation fault when running.
	const int NX = 3; /// number of grid at x direction
	const int NY = 4;	/// number of grid at y direction

	double dx = L/NX;
	double dy = H/NY;

	double A_m[NX*NY][NX*NY];	/// left coefficient matrix
	double B_data[NX*NY];	/// right vector
	
	double a_W, a_E,a_S, a_N, a_P, S_u, S_P;
	
	gsl_spmatrix *A = gsl_spmatrix_alloc(NX*NY ,NX*NY); /// triplet format
	gsl_spmatrix *C;                            /// compressed format 
	gsl_vector *b = gsl_vector_alloc(NX*NY);	/// right hand side vector
	gsl_vector *u = gsl_vector_alloc(NX*NY);        /// solution vector
	
	///************ building coefficient matrix ****************/
	/// nodes arranged by column
	int i,j;
	
	/// internal nodes
	for (i = 1; i < NX-1; i++)
	{
		for (j = 1; j < NY-1; j++)
		{
			a_W = k*dy*W/dx;
			a_E = k*dy*W/dx;
			a_S = k*dx*W/dy;
			a_N = k*dx*W/dy;
			S_u = 0;
			S_P = 0;
			a_P = a_W + a_E + a_N + a_S - S_P;
			gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
			gsl_spmatrix_set(A, i*NY+j, (i-1)*NY+j, -a_W);
			gsl_spmatrix_set(A, i*NY+j, (i+1)*NY+j, -a_E);
			gsl_spmatrix_set(A, i*NY+j, i*NY+j-1, -a_S);
			gsl_spmatrix_set(A, i*NY+j, i*NY+j+1, -a_N);
			gsl_vector_set(b, i*NY+j, S_u);
		}
	}
	
	/// west boundary (excluding corner nodes)
	i = 0;
	for (j = 1; j < NY-1; j++)
	{
		a_W = 0;
		a_E = k*dy*W/dx;
		a_S = k*dx*W/dy;
		a_N = k*dx*W/dy;
		S_u = q*dy*W;
		S_P = 0;
		a_P = a_W + a_E + a_N + a_S - S_P;
		gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
		gsl_spmatrix_set(A, i*NY+j, (i+1)*NY+j, -a_E);
		gsl_spmatrix_set(A, i*NY+j, i*NY+j-1, -a_S);
		gsl_spmatrix_set(A, i*NY+j, i*NY+j+1, -a_N);
		gsl_vector_set(b, i*NY+j, S_u);
	}

	/// east boundary (excluding corner nodes)
	i = NX-1;
	for (j = 1; j < NY-1; j++)
	{
		a_W = k*dy*W/dx;
		a_E = 0;
		a_S = k*dx*W/dy;
		a_N = k*dx*W/dy;
		S_u = 0;
		S_P = 0;
		a_P = a_W + a_E + a_N + a_S - S_P;
		gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
		gsl_spmatrix_set(A, i*NY+j, (i-1)*NY+j, -a_W);
		gsl_spmatrix_set(A, i*NY+j, i*NY+j-1, -a_S);
		gsl_spmatrix_set(A, i*NY+j, i*NY+j+1, -a_N);
		gsl_vector_set(b, i*NY+j, S_u);
	}
	
	/// south boundary (excluding corner nodes)
	j = 0;
	for (i = 1; i < NX-1; i++)
	{
		a_W = k*dy*W/dx;
		a_E = k*dy*W/dx;
		a_S = 0;
		a_N = k*dx*W/dy;
		S_u = T_s*dx*W/(1/h + dx/2/k);
		S_P = -dx*W/(1/h + dx/2/k);
		a_P = a_W + a_E + a_N + a_S - S_P;
		gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
		gsl_spmatrix_set(A, i*NY+j, (i-1)*NY+j, -a_W);
		gsl_spmatrix_set(A, i*NY+j, (i+1)*NY+j, -a_E);
		gsl_spmatrix_set(A, i*NY+j, i*NY+j+1, -a_N);
		gsl_vector_set(b, i*NY+j, S_u);
	}
	
	
	/// north boundary (excluding corner nodes)
	j = NY-1;
	for (i = 1; i < NX-1; i++)
	{
		a_W = k*dy*W/dx;
		a_E = k*dy*W/dx;
		a_S = k*dx*W/dy;
		a_N = 0;
		S_u = 2*k*dx*W*T_n/dy;
		S_P = -2*k*dx*W/dy;
		a_P = a_W + a_E + a_N + a_S - S_P;
		gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
		gsl_spmatrix_set(A, i*NY+j, (i-1)*NY+j, -a_W);
		gsl_spmatrix_set(A, i*NY+j, (i+1)*NY+j, -a_E);
		gsl_spmatrix_set(A, i*NY+j, i*NY+j-1, -a_S);
		gsl_vector_set(b, i*NY+j, S_u);
	}

	/// west-south node
	i = 0, j = 0;
	a_W = 0;
	a_E = k*dy*W/dx;
	a_S = 0;
	a_N = k*dx*W/dy;
	S_u = T_s*dx*W/(1/h + dx/2/k) + q*dy*W;
	S_P = -dx*W/(1/h + dx/2/k);
	a_P = a_W + a_E + a_N + a_S - S_P;
	gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
	gsl_spmatrix_set(A, i*NY+j, (i+1)*NY+j, -a_E);
	gsl_spmatrix_set(A, i*NY+j, i*NY+j+1, -a_N);
	gsl_vector_set(b, i*NY+j, S_u);
	
	/// west-north node
	i = 0, j = NY-1;
	a_W = 0;
	a_E = k*dy*W/dx;
	a_S = k*dx*W/dy;
	a_N = 0;
	S_u = q*dy*W + 2*k*dx*W*T_n/dy;
	S_P = -2*k*dx*W/dy;
	a_P = a_W + a_E + a_N + a_S - S_P;
	gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
	gsl_spmatrix_set(A, i*NY+j, (i+1)*NY+j, -a_E);
	gsl_spmatrix_set(A, i*NY+j, i*NY+j-1, -a_S);
	gsl_vector_set(b, i*NY+j, S_u);
	
	/// east-south node
	i = NX-1, j = 0;
	a_W = k*dy*W/dx;
	a_E = 0;
	a_S = 0;
	a_N = k*dx*W/dy;
	S_u = T_s*dx*W/(1/h + dx/2/k);
	S_P = -dx*W/(1/h + dx/2/k);
	a_P = a_W + a_E + a_N + a_S - S_P;
	gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
	gsl_spmatrix_set(A, i*NY+j, (i-1)*NY+j, -a_W);
	gsl_spmatrix_set(A, i*NY+j, i*NY+j+1, -a_N);
	gsl_vector_set(b, i*NY+j, S_u);

	/// east-north node
	i = NX-1, j = NY-1;
	a_W = k*dy*W/dx;
	a_E = 0;
	a_S = k*dx*W/dy;
	a_N = 0;
	S_u = 2*k*dx*W*T_n/dy;
	S_P = -2*k*dx*W/dy;
	a_P = a_W + a_E + a_N + a_S - S_P;
	gsl_spmatrix_set(A, i*NY+j, i*NY+j, a_P);
	gsl_spmatrix_set(A, i*NY+j, (i-1)*NY+j, -a_W);
	gsl_spmatrix_set(A, i*NY+j, i*NY+j-1, -a_S);
	gsl_vector_set(b, i*NY+j, S_u);

	///************ solve matrix****************/	
	C = gsl_spmatrix_ccs(A);	/// convert to compressed column format
	
	/// now solve the system with the GMRES iterative solver 
  {
    const double tol = 1.0e-6;  /// solution relative tolerance
    const size_t max_iter = 10; /// maximum iterations 
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, NX*NY, 0);
    size_t iter = 0;
    double residual;
    int status;

    gsl_vector_set_all(u, 150.0);	/// initial guess u = 150

	
    //// solve the system A u = b 
    do
      {
        status = gsl_splinalg_itersolve_iterate(C, b, tol, u, work);

        //// print out residual norm ||A*u - b||
        residual = gsl_splinalg_itersolve_normr(work);
        fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

        if (status == GSL_SUCCESS)
          fprintf(stderr, "Converged\n");
      }
    while (status == GSL_CONTINUE && ++iter < max_iter);
  }
	
	///************output solution ****************/
	/// printf("T = \n");
	/// gsl_vector_fprintf(stdout, u, "%g");
	
	/// print out solution data for gnuplot 2D picture
	ofstream outfile;
	string filename;
	filename = "T.out";
	outfile.open(filename);
	for (j = NY-1; j >= 0; j--)
	{
		for (i = 0; i < NX; i++)
		{
			double tt = gsl_vector_get(u,i*NY+j);
			outfile << tt << " ";
		}
		outfile << endl;
	}
	outfile.close();
		
	
	/// print out matrix A and right hand vector b
	/*
	printf("A = \n");
	for (i = 0; i < NX*NY; i++)
	{
		for (j = 0; j < NX*NY; j++)
		{
			double a = gsl_spmatrix_get(A, i, j);
			printf("%g ", a);
		}
		printf("\n");
	}
	printf("b = \n");
	for (i = 0; i < NX*NY; i++)
	{
		double bb = gsl_vector_get(b, i);
		printf("%g\n",bb);
	}
	*/
	gsl_spmatrix_free(A);
	gsl_spmatrix_free(C);
	gsl_vector_free(b);
	gsl_vector_free(u);
	return 0;
}
