#include "../gwpentropy.h"

/********* wrapper function; only pointers can be passed from R **************/

void gwpe(double* H, double* C, double* buf, int *n, int *w, double *qmin, double *qmax, double *dq)
{
	int i;
	double q;

	i = 0;
	for (q = *qmin; q <= *qmax; q += *dq) {
		gwpentropy(H+i, C+i, buf, *n, *w, q);
		i++;
	}
}

