#include <stdio.h>
#include <stdlib.h>
#include "cprini.c"

void tridiag_solve(double *a, int n, double *rhs, double *sol) {
	/** Where a is a 3xn matrix whose elements [0][n-1] and [2][n-1] are 0 **/
	/** n is an integer **/
	/** rhs is an n-length array **/
	/** and sol is an n-length array **/

	if(n <= 0) {
		printf("A matrix cannot have a dimension <= 0.\n");
		return;
	}
	
	
	double *sup = (double *)malloc(n * sizeof(double));
	sup[0] = a[0 * n + 0] / a[1 * n + 0];
	
	for(int k = 1; k < n-1; k++) {
		double denom = a[1 * n + k] - (a[2 * n + k-1] * sup[k-1]);
		sup[k] = a[0 * n + k] / denom;
	}
	/** a[0]' filled **/

	double *rightprime = (double *)malloc(n * sizeof(double));
	rightprime[0] = rhs[0] / a[1 * n + 0];
	for(int k = 1; k < n; k++) {
		double numer = rhs[k] - (a[2 * n + k-1] * rightprime[k-1]);
		double denom = a[1 * n + k] - (a[2 * n + k-1] * sup[k-1]);
		rightprime[k] = numer/denom;	
	}
	/** RHS' filled **/
	
	sol[n-1] = rightprime[n-1];
	for(int k = n-2; k >= 0; k--) {
		sol[k] = rightprime[k] - (sup[k] * sol[k+1]);
	}
	
	/** Sol filled **/

	printf("Solution:\n");
	for(int i = 0; i < n; i++) {
		printf("sol[%d] = %.6f\n", i, sol[i]);
	}

	return;	

}

/**
int main(int argc, char** argv) {
	int n = 4;
	double count = 0.0;
	double *a = (double *)malloc(3 * n * sizeof(double));
	printf("Initial tridiagonal matrix: \n");
	for(int i = 0; i < 3; i++) {
		switch(i) {
			case 0: 
				printf("Superdiagonal:\t ");
				break;
			case 1:
				printf("Diagonal:\t ");
				break;
			case 2:
				printf("Subdiagonal:\t ");
				break;
		}
		for(int j = 0; j < n; j++) {
			count += 1;
			if(j == n-1 && (i == 0 || i == 2)) {
				a[i * n + j] = 0.0;
			} else {
				a[i * n + j] = count;
			}
			printf("a[%d][%d] = %f\t", i, j, a[i * n + j]);
		}
		printf("\n");
	}
	double *rhs = (double *)malloc(n * sizeof(double));
	double *sol = (double *)malloc(n * sizeof(double));
	int tcount = 0.0;
	printf("Right-hand side vector: \n");
	for(int k = 0; k < n; k++) {
		tcount += 10;
		rhs[k] = tcount;
		sol[k] = 0.0;
		printf("rhs[%d] = %f\n", k, rhs[k]);
	}
	
	tridiag_solve(a, n, rhs, sol);	
} **/
