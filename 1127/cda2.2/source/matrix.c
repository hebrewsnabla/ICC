#include <stdio.h>
#include <math.h>
#include <stdlib.h>


/*multiply two matrices************************************************/
void matmult(float **ma, float **mb, int dim, float **mc) {
        int i, j, k;

        for (i=1; i <= dim; i++) {
                for (j=1; j <= dim; j++) {
                        mc[i][j] = 0.0;
                        for (k=1; k <= dim; k++) {
                                mc[i][j] += ma[i][k] * mb[k][j];
                        }
                }
        }
}

/*build adjoint of matrix (actually, transpose it)*********************/
void mattransp(float **ma, int dim, float **mb) {
	int i, j;

	for (i=1; i <= dim; i++) {
		for (j=1; j <= dim; j++) {
		mb[j][i] = ma[i][j];
		}
	}
}

/*copy matrix**********************************************************/
void matcopy(float **ma, int dim, float **mb) {
        int i, j;

        for (i=1; i <= dim; i++) {
                for (j=1; j <= dim; j++) {
			mb[i][j] = ma[i][j];
                }
        }
}

/*add two matrices*****************************************************/
void matadd(float **ma, float **mb, int dim, float **mc) {
	int i, j;

	for (i=1; i <= dim; i++) {
		for (j=1; j <= dim; j++) {
			mc[i][j] = ma[i][j] + mb[i][j];
		}
	}
}


/*direct sum of two matrices*******************************************/
void matsum(float **ma, int dima, float **mb, int dimb, float **mc) {
	int i, j;

	for (i=1; i <= dima; i++) {
		for (j=1; j <= dima; j++) {
			mc[i][j] = ma[i][j];
		}
	}
	for (i=1; i <= dimb; i++) {
		for (j=1; j <= dimb; j++) {
			mc[dima+i][dima+j] = mb[i][j];
		}
	}
}


/*clear matrix*********************************************************/
void matclear(float **ma, int dim) {
	int i, j;

        for (i=1; i <= dim; i++) {
                for (j=1; j <= dim; j++) {
                        ma[i][j] = 0.0;
                }
        }
}


/***********************************************************************/
/*build inverse of matrix (routine from numerical recipes)**************/
/***********************************************************************/
void matinv(float **a, float **y, int n) {
	void lubksb(float **, int, int *, float *);
	void ludcmp(float **, int, int *, float *);
	float d, *col;
	int i, j, *indx;

	indx = (int *) calloc((n+1), sizeof(int));
	col = (float *) calloc((n+1), sizeof(float));

	ludcmp(a, n, indx, &d);
	for (j=1; j <= n; j++) {
		for (i=1; i <= n; i++)  col[i]=0.0;
		col[j]=1.0;
		lubksb(a, n, indx, col);
		for (i=1; i <= n; i++) y[i][j]=col[i];
	}
}

/*** lubksb ***/
void lubksb (float **a, int n, int *indx, float *b) {
	int i, ii=0, ip, j;
	float sum;

	for (i=1; i <= n; i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii; j <= i-1; j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n; i >= 1; i--) {
		sum=b[i];
		for (j=i+1; j <= n; j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

/*** ludcmp ***/
void ludcmp (float **a, int n, int *indx, float *d) {
        float *vector(int, int);
	void free_vector(float *, int);
	int i, imax, j, k;
	float big, dum, sum, temp, *vv;

	vv=vector(1, n);
	*d=1.0;
	for (i=1; i <= n; i++) {
		big = 0.0;
		for (j=1; j <= n; j++)
			if ((temp = fabs((double) a[i][j])) > big) big = temp;
		if (big == 0.0) {
			fprintf(stderr, "Sorry, matrix is singular...");
			getchar();
			exit(0);
		}
		vv[i] = 1.0/big;
	}
	for (j=1; j <= n; j++) {
		for (i=1; i<j; i++) {
			sum = a[i][j];
			for (k=1; k < i; k++) sum -= a[i][k] * a[k][j];
			a[i][j]=sum;
		}
		big = 0.0;
		for (i=j; i <= n; i++) {
			sum = a[i][j];
			for (k=1; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs((double) sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k=1; k <= n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = 1.0e-10;
		if (j != n) {
			dum = 1.0/(a[j][j]);
			for (i=j+1; i <= n; i++) a[i][j] *= dum;
		}
	}
	free_vector(vv, 1);
}

/*** vector ***/
float *vector(int nl, int nh) {
	float *v;

	v = (float *) malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) fprintf(stderr, "allocation failure in vector()");
	return (v-nl);
}

/*** free_vector***/
void free_vector(float *v, int nl) {
	free((char *) (v+nl));
}

