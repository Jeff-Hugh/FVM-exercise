#include <iostream>
#include <gsl/gsl_linalg.h>
using namespace std;

int main()
{
	const double L = 0.5; // total length
	const double A = 1e-2; // sectional area
	const double TA = 100; // boundary condition A
	const double TB = 500; // boundary condition B
	const double k = 1000; // thermal conductivity

	const int N = 10; // number of grid

	double dx = L/N;

	double diag_data[N]; // diagonal of coefficient matrix
	double e_data[N-1];	// super-diagonal of coefficient matrix
	double f_data[N-1];	// sub-diagonal of coefficient matrix
	double b_data[N];	// right vector
	
	
	/************ building coefficient matrix ****************/
	double temp;
	temp = k/dx*A;

	for (int i=1; i<N-1; i++)
	{
		e_data[i] = -temp;
		diag_data[i] = 2*temp;
		f_data[i-1] = -temp;
		b_data[i] = 0;
	}
	
	// A boundary, i = 0
	diag_data[0] = 3*temp;
	e_data[0] = -temp;
	b_data[0] = 2*temp*TA;

	// B boundary, i = N-1
	f_data[N-2] = -temp;
	diag_data[N-1] = 3*temp;
	b_data[N-1] = 2*temp*TB;

	
	/************ solve matrix****************/

	gsl_vector_view diag = gsl_vector_view_array(diag_data, N);
	gsl_vector_view b = gsl_vector_view_array(b_data, N);
	gsl_vector_view e = gsl_vector_view_array(e_data, N-1);
	gsl_vector_view f = gsl_vector_view_array(f_data, N-1);
	gsl_vector *T = gsl_vector_alloc(N); 		// temperature
	
	gsl_linalg_solve_tridiag(&diag.vector, &e.vector, &f.vector, &b.vector, T);
	
	
	/************ print out ****************/
	printf("T = \n");
	gsl_vector_fprintf(stdout, T, "%g");
	
	return 0;
}
