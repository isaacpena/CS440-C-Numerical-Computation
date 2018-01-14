#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

void jacobi(double *a, int n, double *s, double *u, double *v) {
	
	double mach_eps_mult = .000000001;

	/* Small multiple of the machine epsilon for doubles */

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			u[i * n + j] = a[i * n + j];
			if(i == j) {
				v[i * n + j] = 1.0;
			}
		}
	}
	/* U = A; V = I_n */

	double convergence = mach_eps_mult + 1.0;
	while(convergence > mach_eps_mult) {
		convergence = 0;
		for (int j = 1; j < n; j++) {
			for(int i = 0; i < j; i++) {
				double upperleft = 0.0;
				double lowerright = 0.0;
				double nondiag = 0.0;
				for(int k = 0; k < n; k++) {
					upperleft += pow(u[k * n + i], 2);
					lowerright += pow(u[k * n + j], 2);
					nondiag += (u[k * n + i] * u[k * n + j]);
				}
				/* UL of 2x2 matrix is Sum(U^2_(k, i)) for all k */
				/* LR of 2x2 matrix is Sum(U^2_(k, j)) */
				/* nondiags are Sum(U_(k, i) * U_(k_j)) */

				if(nondiag == 0) {
					continue;
				}
	
				double comp = fabs(nondiag) / (sqrt(upperleft * lowerright));
				if(convergence < comp) {
					convergence = comp;
				}

				double beta = (lowerright - upperleft) / (2 * nondiag);
					
				double tan;
				if (beta < 0) {
					tan = -1 / (fabs(beta) + sqrt(beta*beta + 1));
				} else {
					tan = 1 / (fabs(beta) + sqrt(beta*beta + 1));
				}
				
				double cosine = 1 / (sqrt(1 + (tan * tan)));
				double sine = tan * cosine;
	
				for(int k = 0; k < n; k++) {
					tan = u[k * n + i];
					u[k * n + i] = cosine * tan - (sine * u[k * n + j]);
					u[k * n + j] = sine * tan + (cosine * u[k * n + j]);
					
					/* Rotate U */

					tan = v[k * n + i];
					v[k * n + i] = cosine * tan - (sine * v[k * n + j]);
					v[k * n + j] = sine * tan + (cosine * v[k * n + j]);
					
					/* Rotate V */			
	
				}

			}
		}	

	}
	
	for(int j = 0; j < n; j++) {
		double sum_of_squares = 0.0;
		for(int k = 0; k < n; k++) {
			sum_of_squares += (u[k * n + j] * u[k * n + j]);
		}
		double norm = sqrt(sum_of_squares);
		s[j] = norm;
		/* Singular values are vector norms of the columns of u */
		for(int k = 0; k < n; k++) {
			if(s[j] != 0) {
				u[k * n + j] = u[k * n + j] / s[j];
			}
		}
	}


	for(int i = 0; i < n-1; i++) {
		int min = i; 
		for(int j = i + 1; j < n; j++) {
			if(s[j] > s[min]) {
				min = j;
			}
		}
		if(min != i) {
			double swap = s[i];
			s[i] = s[min];
			s[min] = swap;
			if(swap < 0.0) {
				s[min] = -1.0 * swap;
			}
				

			for(int k = 0; k < n; k++) {
				double uswap = u[k * n + i];
				u[k * n + i] = u[k * n + min];
				u[k * n + min] = uswap;
				if (swap < 0.0) {
					u[k * n + min] = -1.0 * uswap;
				}
				double vswap = v[k * n + i];
				v[k * n + i] = v[k * n + min];
				v[k * n + min] = vswap;
				if (swap < 0.0) {
					v[k * n + min] = -1.0 * vswap;
				}
			}
		}
	}

	/* Selection sort on singular values */

	/* UNCOMMENT TO PRINT U, \SIGMA, V */
	printf("U:\n");
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++) {
			printf("%f   ", u[i * n + j]);
		}	
		printf("\n");
	}
	
	printf("\n");
	
	printf("Sigma:\n");
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(i == j) {
				printf("%f   ", s[i]);
			} else {
				printf("%f   ", 0.0);
			}
		}
		printf("\n");
	}

	printf("\n");
	
	printf("V:\n");
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			printf("%f   ", v[i * n + j]);
		}	
		printf("\n");
	}		

}


int main(int argc, char** argv){
	int n = 5;
	
	double *a = malloc(sizeof(double) * n * n);
	double *u = malloc(sizeof(double) * n * n);
	double *v = malloc(sizeof(double) * n * n);
	double *s = malloc(sizeof(double) * n);
	
	for(int i = 0; i < n; i++){
		s[i] = 0.0;
		for(int j = 0; j < n; j++) {
			a[i * n + j] = 0.0;
			u[i * n + j] = 0.0;
			v[i * n + j] = 0.0;
		}
	}

	/*sparse matrix */
	/*
 * 	   1 0 0 0 2 
 * 	   0 0 3 0 0
 * 	   0 0 0 0 0
 * 	   0 1 0 0 0
 * 	   0 0 2 0 0
 * 	*/
	
	a[0 * n + 0] = 1.0;
	a[0 * n + 4] = 2.0;
	a[1 * n + 2] = 3.0;
	a[3 * n + 1] = 1.0;
	a[4 * n + 2] = 2.0;

	/* dense matrix */
	/* 	
	a[0 * n + 0] = 1;
	a[0 * n + 1] = 3;
	a[0 * n + 2] = 2;
	a[1 * n + 0] = 5;
	a[1 * n + 1] = 6;
	a[1 * n + 2] = 4;
	a[2 * n + 0] = 7;
	a[2 * n + 1] = 8;
	a[2 * n + 2] = 9;
	*/
	jacobi(a, n, s, u, v);
}
