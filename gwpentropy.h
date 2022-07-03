#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define HUGE	1e80

/********************************************************************************************																							*
 * function: symbolize_pompe																*
 * description: symbolize sequence with pompe's original scheme.							*
 * parameters:																				*
 *			buf:		series.																*
 *			word:		sequence.															*
 *			w:			size of sequence.													*
 *******************************************************************************************/
void symbolize_pompe(double *buf, int *word, int w)
{
  /* local variables */
  int k1, k2, kmin; 
  double buf1, buf2, tmp[100];

  memcpy(tmp, buf, w*sizeof(double));					/* copy of sequence */

  for(k1=0; k1<w; k1++) {								/* loop over sequence */

    kmin = 0;
    buf1 = tmp[0];										/* element 1 of sequence */

    for(k2=0; k2<w; k2++) {								/* loop over sequence */
      buf2 = tmp[k2];									/* element 2 of sequence */
	  if(buf1>buf2) {									/* case: element 1 > element 2 */
	    buf1 = buf2;
		kmin = k2;
	  }
	}

	word[k1] = kmin;		/* rank of the smallest entry*/
	tmp[kmin] = HUGE;		/* make sure it doesn't apper again*/
  }

}

/********************************************************************************************																						*
 * function: searchinsert																	*
 * description: Search for matching sequence in table. Insert if new sequence in table.		*
 * parameters:																				*
 *			table:		table of sequences.													*
 *			word:		sequence.															*
 *			w:			size of sequence.													*
 *			num:		number of sequences in table.										*
********************************************************************************************/
int searchinsert(int *table, int *word, int w, int *num)
{
  /* local variables */
  int i, k, offset; 

  for(i=0; i<*num; i++) {									/* loop over all sequences stored in table */
	offset = i*w;
	  
	for(k=0; k<w; k++)										/* loop over current sequence */
	  if(word[k]!=table[offset+k])							/* exit if sequence pair do not match */
	    break;

	if(k==w)												/* case: pair of sequences are the same */
	  return i;					
  }

  memcpy(&table[*num*w], word, w*sizeof(int));
  *num = *num + 1;
  
  return *num - 1;												/* return whether match occured */
}

/********************************************************************************************
* function: pow0																			*
* description: Computes pow(x,y), also valid when x=0 and y < 0.							*
* parameters:																				*
*			x:		input x.																*
*			y:		power.																	*
********************************************************************************************/
double pow0(double x, double y)
{
	if (x == 0) 
		return 0;
	else return pow(x, y);
}

/********************************************************************************************
* function: mean																			*
* description : calculate mean of series.													*
* parameters :																				*
* buf : series.																				*
* n : size of series.																		*
********************************************************************************************/
double mean(double* buf, int n)
{
	int i;
	double mu = 0;


	for (i = 0; i < n; i++)						/* loop over series */
		mu += buf[i];
	mu = mu / n;

	return mu;
}


/********************************************************************************************
 * function: variance																		*
 * description: calculate variance of series.												*
 * parameters:																				*
 *			buf:		series.																*
 *			n:			size of series.														*
********************************************************************************************/
double variance(double* buf, int n)
{
	int i;
	double mu, var = 0;

	mu = mean(buf, n);							/* calculate mean */

	for (i = 0; i < n; i++)						/* loop over series */
		var += (buf[i] - mu) * (buf[i] - mu);
	var = var / n;

	return var;
}

/********************************************************************************************
 * function: factorial																		*
 * description: compute factorial.															*
 * parameters:																				*
 *			n:			factorial order.													*
********************************************************************************************/
int factorial(int n)
{
	int i, val;

	val = 1;
	for (i = 1; i <= n; i++)
		val *= i;

	return val;
}

/********************************************************************************************
 * function: gwpentropy																		*
 * description: calculates generalized weighted permutation entropy							*
 * parameters:																				*
 *			H:			entropy.															*
 *			C:			complexity.															*
 *			buf:		series.																*
 *			n:			size of series.														*
 *			w:			size of word/sequence (pattern).									*
 *			q:			scale coefficient.													*
********************************************************************************************/
void gwpentropy(double *H, double *C, double *buf, int n, int w, double q)
{
  /* local variables */
  int i, k, num; 
  double f, h, c, var, nn, logn, ptot;
  int    *word=NULL,										/* word (pattern) */
	     *table=NULL;										/* table of sequences (patterns) */
  double *p=NULL;											/* probabilities */
  double maxvar = 0, svar;

  /* allocate internal buffers */
  word  = (int *)malloc(w*sizeof(int));						/* allocate memory for word */
  table = (int *)malloc(n*w*sizeof(int));					/* allocate memory for sequence table */
  p = (double *)calloc(n, sizeof(double));					/* allocate memory for frequencies */

  /* initialize params */
  k = -1;													/* reset index */
  num = 0;													/* reset number of unique patterns (size of sequence table) */
  ptot = 0.0;												/* reset total probabilities */
  nn = (double)factorial(w);								/* w! */

  /* compute probabilities */
  for(i=0; i<n-w+1; i++) {									/* loop over series */
	symbolize_pompe(&buf[i], word, w);						/* symbolize a sequence */
	k = searchinsert(table, word, w, &num);					/* search for sequence in table */
	var = pow0(svar=variance(&buf[i], w),(double)(q/2));	/* pi = {(1/w)sum[(x - xu)^2]}^(q/2)*/
	if (maxvar < var)
		maxvar = var;
	p[k] += var;											/* pk += pi*/
  }

  /* normalize probabilities */
  for(i=0; i<num; i++)
	ptot += p[i];
  
  for(i=0; i<num; i++)
    p[i] = p[i]/ptot;										/* pk = pk/sum(pk) */
  
  /* calculate shannon entropy */
  h = 0.0;
  for(i=0; i<num; i++)										/* loop over frequencies */
    h-= p[i]*log(p[i]);										/* f*log(f) */

  /* calculate complexity */
  c = 0.0;
  for (i = 0; i<num; i++) {									/* loop over frequencies */
	  f = (p[i] + 1.0/ nn)/2.0;								/* probability */
	  c -= f*log(f);										/* f*log(f) */
  }
  f  = 1.0/(2.0*nn);
  c -= (nn - num)*f*log(f);
  logn = log(nn);

  /* store results */
  *C = -2.0*(c - h/2.0 - logn/2.0)*(h/logn)/((nn+1.0)/nn*log(nn+1)-2.0*log(2.0*nn)+logn);
  *H = h/logn;

  /* free internal buffers */
  if(p) free(p);
  if(word) free (word);
  if(table) free (table);

}