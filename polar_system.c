#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_const_mksa.h>

//TODO: MAKE IT WORK!

struct problem_arguments{
//Input values:
	double m;
	double h;
	double V1;
	double V2;
	double V3;
	double a;
	double b;
//Calculated values:
	gsl_vector *E;
	gsl_vector_complex *k1;
	gsl_vector_complex *k2;
	gsl_vector_complex *k3;
	gsl_vector *det_M;
	gsl_matrix *M;	
};

struct problem_arguments data; //Global data

void init_data()
{			//Real values:
	data.m = 1;	//GSL_CONST_MKSA_MASS_ELECTRON
	data.h = 1;	//GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR 
	data.V1 = 1;
	data.V2 = -2;
	data.V3 = 3;
	data.a = 1;
	data.b = 2;

	data.E = NULL;
	data.k1 = NULL;
	data.k2 = NULL;
	data.k3 = NULL;
	data.det_M = NULL;
	data.M = NULL;

}

void print_data()
{
	int row1 = 16;
	int row2 = 8;
	
	printf("Program input data:\n");
	printf("-----------------------------------\n");
	printf("| %-*s :\t %*.3lf |\n", row1, "Mass", row2, data.m);
	printf("| %-*s :\t %*.3lf |\n", row1, "Planck constant", row2, data.h);
	printf("| %-*s :\t %*.3lf |\n", row1, "Potential V1", row2, data.V1);
	printf("| %-*s :\t %*.3lf |\n", row1, "Potential V2", row2, data.V2);
	printf("| %-*s :\t %*.3lf |\n", row1, "Potential V3", row2, data.V3);
	printf("| %-*s :\t %*.3lf |\n", row1, "Width a", row2, data.a);
	printf("| %-*s :\t %*.3lf |\n", row1, "Width b", row2, data.b);
	printf("-----------------------------------\n");
}

void print_graph()
{
	//Maybe in future (gnuplot interface in ANSI C)
}

int input_params(int argc, char *argv[])
{
	int i=1;
	
	//NOTE: Unsafe for invalid input (stack smash), need fix 	
	for(i;i < argc; i++)
	{
		//Mass and Planck constant:
		if(!strncmp(argv[i],"-m",2))
                {
                        data.m = strtod(argv[i+1],NULL);
                }

		if(!strncmp(argv[i],"-h",2))
                {
                        data.h = strtod(argv[i+1],NULL);
                }

		//Potentials input:
		if(!strncmp(argv[i],"-V1",3))
		{
			data.V1 = strtod(argv[i+1],NULL);
		}

		if(!strncmp(argv[i],"-V2",3))
                {
			data.V2 = strtod(argv[i+1],NULL);
                }

		if(!strncmp(argv[i],"-V3",3))
                {
			data.V3 = strtod(argv[i+1],NULL);
                }
		
		//With (radius) of potential box 
		if(!strncmp(argv[i],"-a",2))
                {
                        data.V1 = strtod(argv[i+1],NULL);
                }

		if(!strncmp(argv[i],"-b",2))
                {
                        data.V1 = strtod(argv[i+1],NULL);
                }
		
	}
	return 0;
}

void generate_E_range(double beg, double end, int precision)
{
	int i;
	double difference;
	data.E = gsl_vector_alloc(precision);
	gsl_vector_set_all(data.E, 0);

	difference = (end - beg) / (double) precision;
	for(i=0; i < precision; i++)
	{
		gsl_vector_set(data.E, i, beg + difference*i);
	}
}

void generate_k_vectors()
{
	int i;
	int size;
	gsl_complex calculation;
	double planck_etc;
	size = data.E->size;

	data.k1 = gsl_vector_complex_alloc(size);
	data.k2 = gsl_vector_complex_alloc(size);
	data.k3 = gsl_vector_complex_alloc(size);

	planck_etc = 2 * data.m / (data.h * data.h);

	for(i=0; i < size; i++)
	{
		calculation = gsl_complex_sqrt_real( (gsl_vector_get(data.E, i) + data.V1) * planck_etc);
		/* if(calculation != calculation) //Checking for NaN if true:
                        calculation = 1;
		*/
		gsl_vector_complex_set(data.k1, i, calculation);
	}	

	for(i=0; i < size; i++)
        {
                calculation = gsl_complex_sqrt_real( (gsl_vector_get(data.E, i) - data.V2) * planck_etc);
                /*if(calculation != calculation)
			calculation = 1;
		*/
		gsl_vector_complex_set(data.k2, i, calculation);
        }

	for(i=0; i < size; i++)
        {
                calculation = gsl_complex_sqrt_real( (gsl_vector_get(data.E, i) - data.V3) * -1 * planck_etc);
		/*if(calculation != calculation)
                        calculation = 0;
		*/
                gsl_vector_complex_set(data.k3, i, calculation);
        }
}

//Generates M matrix based on E[i] condition 
//Works only for V1 < V2 < V3 and V1 < E < V3
void generate_M(int i)
{
	if(data.M == NULL)
		data.M = gsl_matrix_alloc(4,4);
	gsl_matrix_set_zero(data.M);

	/* Hardcode time:
	 * 16 cells of matrix M
	 * filled with bessel functions
	 * one by one.
	 * Jn(x) is Bessel function
	 * Yn(x) is Neumann function
	 * In future can by replaced by
	 * iterative method... maybe 
	 */

	double cell_value;
	//NOTE: k is complex - should it be M too ?
 
	//Element 1,1: J1(|k1|*a) 
	cell_value = gsl_sf_bessel_J1( gsl_complex_abs(gsl_vector_complex_get(data.k1, i)) * data.a);
	gsl_matrix_set(data.M,0,0, cell_value);

	//Element 1,2: -J1(|k2|*a)
	cell_value = -1.0 * gsl_sf_bessel_J1( gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.a);
	gsl_matrix_set(data.M,0,1, cell_value);

	//Element 1,3: -Y1(|k2|*a)
	cell_value = -1.0 * gsl_sf_bessel_Y1( gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.a);
	gsl_matrix_set(data.M,0,2, cell_value);

	//Element 1,4: 0
	cell_value = 0.0;
	gsl_matrix_set(data.M,0,3, cell_value);

	//Element 2,1: 1/2*|k1|*( J0(|k1|*a) - J2(|k1|*a) )
	cell_value = 0.5 * gsl_complex_abs(gsl_vector_complex_get(data.k1, i)) * ( gsl_sf_bessel_J0( gsl_complex_abs(gsl_vector_complex_get(data.k1, i)) * data.a) - gsl_sf_bessel_Jn(2, gsl_complex_abs(gsl_vector_complex_get(data.k1, i)) * data.a));
	gsl_matrix_set(data.M,1,0, cell_value);

	//Element 2,2: -1/2*|k2|*( J0(|k2|*a) - J2(|k2|*a) ) 
	cell_value = -0.5 * gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * ( gsl_sf_bessel_J0(gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.a) - gsl_sf_bessel_Jn(2, gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.a));
	gsl_matrix_set(data.M,1,1, cell_value);

	//Element 2,3: -1/2*|k2|*( Y0(|k2|*a) - Y2(|k2|*a) )
	cell_value = -0.5 * gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * ( gsl_sf_bessel_Y0(gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.a) - gsl_sf_bessel_Yn(2, gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.a));
	gsl_matrix_set(data.M,1,2, cell_value);

	//Element 2,4: 0
	cell_value = 0;
	gsl_matrix_set(data.M,1,3, cell_value);

	//Element 3,1: 0
	cell_value = 0;
	gsl_matrix_set(data.M,2,0, cell_value);

	//Element 3,2: J1(|k2|*b)
	cell_value = gsl_sf_bessel_J1( gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.b);
	gsl_matrix_set(data.M,2,1, cell_value);
	
	//Element 3,3: Y1(|k2|*b)
	cell_value = gsl_sf_bessel_Y1( gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.b);
	gsl_matrix_set(data.M,2,2, cell_value);
	
	//Element 3,4: -( J1(|k3|*b) + i*Y1(|k3|*b))
	cell_value = -1.0 * ( gsl_sf_bessel_J1( gsl_complex_abs(gsl_vector_complex_get(data.k3, i)) * data.b) + 1 * gsl_sf_bessel_Y1( gsl_complex_abs(gsl_vector_complex_get(data.k3, i)) * data.b) ); //NOTE: replace 1 to imaginary
	gsl_matrix_set(data.M,2,3, cell_value);
	
	//Element 4,1: 0
	cell_value = 0;
	gsl_matrix_set(data.M,3,0, cell_value);
	
	//Element 4,2: 1/2*|k2|*( J0(|k2|*b) - J2(|k1|*b) )
	cell_value = 0.5 * gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * ( gsl_sf_bessel_J0(gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.b) - gsl_sf_bessel_Jn(2, gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.b));
	gsl_matrix_set(data.M,3,1, cell_value);
	
	//Element 4,3: 1/2*|k2|*( Y0(|k2|*b) - Y2(|k2|*b) )
	cell_value = 0.5 * gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * ( gsl_sf_bessel_Y0(gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.b) - gsl_sf_bessel_Yn(2, gsl_complex_abs(gsl_vector_complex_get(data.k2, i)) * data.b));
	gsl_matrix_set(data.M,3,2, cell_value);
	
	//Element 4,4: -1/2*|k3|*( J0(|k3|*b) - J2(|k3|*b) + Y0(|k3|*b) - Y2(|k3|*b) )
	cell_value = -0.5 * gsl_complex_abs(gsl_vector_complex_get(data.k3, i)) * ( gsl_sf_bessel_J0(gsl_complex_abs(gsl_vector_complex_get(data.k3, i)) * data.b) - gsl_sf_bessel_Jn(2, gsl_complex_abs(gsl_vector_complex_get(data.k3, i)) * data.b) + 1 * gsl_sf_bessel_Y0(gsl_complex_abs(gsl_vector_complex_get(data.k3, i)) * data.b) - 1 * gsl_sf_bessel_Jn(2, gsl_complex_abs(gsl_vector_complex_get(data.k3, i)) * data.b)); //NOTE: replace 1 to imaginary
	gsl_matrix_set(data.M,3,3, cell_value);
	
}


void calculations(int precision)
{
	int i;
	generate_E_range(data.V1, data.V3, precision);
	//gsl_vector_fprintf(stdout, data.E, "%lf\n");
	generate_k_vectors();
	//gsl_vector_fprintf(stdout, data.k1, "%lf\n");
	//gsl_vector_complex_fprintf(stdout, data.k2, "%lf\n");
	//gsl_vector_fprintf(stdout, data.k3, "%lf\n");

	data.det_M = gsl_vector_alloc(precision);
	gsl_vector_set_all(data.det_M, 0);
	
	gsl_permutation *p = gsl_permutation_alloc(4);
 	int s;

	for(i=0; i < precision; i++)
	{		
		generate_M(i);
		//gsl_matrix_fprintf(stdout, data.M, "%.10e"); //Debug printf
		//printf("-----------------\n"); //Debug printf
		gsl_linalg_LU_decomp(data.M, p, &s);
		gsl_vector_set(data.det_M, i, gsl_linalg_LU_det(data.M, s) );
	}
	//printf("Values of detM(E):\n");
	//gsl_vector_fprintf(stdout, data.det_M, "%.10e"); //Debug printf
	gsl_permutation_free(p);
	
}

void print_output()
{
	printf("\nValues of detM(E):\n");
        gsl_vector_fprintf(stdout, data.det_M, "%.10e"); //Debug printf
}

void free_all_data()
{
	gsl_vector_free(data.E);
	gsl_vector_complex_free(data.k1);
	gsl_vector_complex_free(data.k2);
	gsl_vector_complex_free(data.k3);
	gsl_matrix_free(data.M);
	gsl_vector_free(data.det_M);
}

int main(int argc, char *argv[])
{
	init_data();
	input_params(argc, argv);
	print_data();
	calculations(30); //Set precision
	print_output();	
	free_all_data();
	return 0;	
}
