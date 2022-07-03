#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../gwpentropy.h"

double buf[10000000];

char fnom[][64]={
"fbm_h=0.1_n=2^15",
"fbm_h=0.5_n=2^15",
"fbm_h=0.9_n=2^15",
"henonmap",
"logisticmap",
"skewtentmap",
"gaussian_2^15",
"uniform_2^15",
"aami31a",
};

/*********************************************************************************/
int load(char* filename, double* buf)
{
	/* local variables */
	long size = 0L;
	double val;
	FILE* fp;

	fp = fopen(filename, "r+");						/* open file */

	while (fscanf(fp, "%lf", &val) == 1) {			/* loop until end of file */
		size++;
		buf[size - 1] = val;							/* store data */
	}

	fclose(fp);										/* close file */
	return (int)(size);
}

/*********************************************************************************/

void main() {

  int i, w, m, n;
  double r;
  double H, C, q, qmin, qmax, dq;
  int window = 250, jump = 1, total, first, offset = 0, nlin=0;
  FILE *h;
  char fname[256];

  w = 6;
  qmin = -10.0; qmax = 10.0; dq = 0.1;

  for (i = 0; i < 9; i++) {
	  sprintf(fname, "../data/%s.dat", fnom[i]);
	  total = n = load(fname, buf);
	  sprintf(fname, "../resultsw6\\%s_gpentpy.txt", fnom[i]);
	  h = fopen(fname, "w+");
	  fprintf(h, "H\tC\tq\n");
	  for (q = qmin; q <= qmax; q += dq) {
		  gwpentropy(&H, &C, buf, n, w, q);
		  printf("%s\t%f\t%f\t%f\n", fnom[i], q, H, C);
		  fprintf(h, "%f\t%f\t%f\n", H, C, q);	/* x,y,z print order*/
	  }
	  fclose(h);
  }

}