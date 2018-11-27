#include <iostream>
#include <gsl/gsl_linalg.h>
using namespace std;

int main()
{
	const double L = 0.1; /// total length
	const double A = 1; /// sectional area
	const double phi_0 = 1; /// boundary condition A
	const double phi_L = 0; /// boundary condition B
	const double rho = 1; /// dencity, kg/m^3
	const double gamma = 0.1;    /// diffusion coefficience, kg/(m*s)
	const double u = 0.1;    /// velocity, m/s

	const int N = 5;  /// number of grid 
	double dx = L/N;
 
	double diag_data[N]; // diagonal of coefficient matrix
	double e_data[N-1];	// super-diagonal of coefficient matrix
	double f_data[N-1];	// sub-diagonal of coefficient matrix
	double b_data[N];	// right vector
	
	
	/************ building coefficient matrix ****************/
	double F = rho*u;
	double D = gamma/dx;

	for (int i=1; i<N-1; i++)
	{
		e_data[i] = -(D - F/2);
		diag_data[i] = 2*D;
		f_data[i-1] = -(D + F/2);
		b_data[i] = 0;
	}
	
	// A boundary, i = 0
	diag_data[0] = 3*D + F/2;
	e_data[0] = -(D - F/2);
	b_data[0] = (2*D + F)*phi_0;

	// B boundary, i = N-1
	f_data[N-2] = -(D +F/2);
	diag_data[N-1] = 3*D - F/2;
	b_data[N-1] = (2*D - F)*phi_L;

	
	/************ solve matrix****************/

	gsl_vector_view diag = gsl_vector_view_array(diag_data, N);
	gsl_vector_view b = gsl_vector_view_array(b_data, N);
	gsl_vector_view e = gsl_vector_view_array(e_data, N-1);
	gsl_vector_view f = gsl_vector_view_array(f_data, N-1);
	gsl_vector *phi = gsl_vector_alloc(N); 		
	
	gsl_linalg_solve_tridiag(&diag.vector, &e.vector, &f.vector, &b.vector, phi);
	
	
	/************ print out ****************/
	printf("phi = \n");
	gsl_vector_fprintf(stdout, phi, "%g");
	
	return 0;
}
