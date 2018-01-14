#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

double f(double x) {
	return (pow(fabs(sin(x)), -.25) * cos((x - 7) / 100)); 
}

double numquadAux(double a, double b, double epsilon, double points) {
	double interval = fabs(b-a) / points;
	
	double low = a;
	double high = a + interval;
	
	double sum = 0.0;	
	double i = 0.0;

	while(i < points) {
		double c = (low+high) / 2;
		double h = (high - low);
		double fl = f(low);
		double fc = f(c);
		double fh = f(high);

		
		double intval = (h/6) * (fl + (4 * fc) + fh);
		sum += intval;
	
		low += interval;
		high += interval;
		i++;
	}
	return sum;
}

double upnqAux(double a, double b, double epsilon, double points) {
	double interval = fabs(b - a) / points;

	long double flow = a;
	long double fhigh = a + interval;
	long double low = acos(-1.0) - fhigh;
	long double high = acos(-1.0) - flow;
	
	double sum = 0.0;
	double i = 0.0;

	while(i < points) {
		double c = (low+high) / 2;
		double h = (high - low);
		double fl = f(low);
		double fc = f(c);
		double fh = f(high);
	
		double intval = (h/6) * (fl + (4 * fc) + fh);
		sum += intval;

		low -= interval;
		high -= interval;
		i++;	
	}
	return sum;
}


void numquad(double *val) {
	double epsilon = .00000000000001;

	double b= acos(-1.0) / 2.0;

	double a = b * (.50);
	
	double sum1 = 0.0;
	double sum2 = 0.0;

	double res1 = 1.0;
	// Left Side

	double res2 = 1.0;
	//Rightside 
	
	while(res1 > epsilon) {
		res1 = numquadAux(a, b, epsilon, 100000);
		
		b = a;
		a = a * (.5);
		
		sum1 += res1;
	}

	b = acos(-1.0) / 2.0;
	a = b * (.50);

	while(res2 > epsilon) {
		res2 = upnqAux(a, b, epsilon, 100000);
		b = a;
		a = a / 2.0;
		
		sum2 += res2;
	}

	*val = sum1 + sum2;
}


int main(int argc, char** argv) {
	double *val = malloc(sizeof(int));
	*val = 0;
	numquad(val);
	printf("%1.10f\n", *val);
}
