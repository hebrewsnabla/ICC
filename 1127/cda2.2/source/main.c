#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "readgaus.h"
#include "matrix.h"
#include "output.h"

int  main(int argc, char **argv) {

/***variable declarations***/
	float **mat1, **mat2, **mat3, **mat4, **mat5;
	int alfa1, alfa2, alfa3;
	int beta1, beta2, beta3;
	int basi1, basi2, basi3;
	int i, j, m, n;
	char string[200];
	float d, b, r, s, df, bf, rf, sf, ea, eb, ec, mue, eta, moni;
	time_t zeit1, zeit2;
	double enc[1000], enf[1000];




/***Test number of arguments from command line************************/
	if (argc != 4) {
	fprintf(stderr, "Wrong number of arguments.\n");
	fprintf(stderr, "Usage: %s molecule part_A part_B\n", argv[0]);
	exit(1);
	}

/***Print header and start time***************************************/
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("      CHARGE DECOMPOSITION ANALYSIS   <%s>\n\n", argv[0]);
	printf("             S. Dapprich and G. Frenking\n");
	printf("         Philipps-Universitaet Marburg, 1995\n");
	printf("             \n");
	printf(" (Last modification by Moritz von Hopffgarten, 2011)\n");
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("\n");
	time(&zeit1);
	printf("/Execution started at %s\n", ctime(&zeit1));

/***Read formulas********************************************************/
	printf("--- Chemical system information ---\n\n");
	readformula(argv[1], string);
	printf(" Complex  C:  %s", string);
        readformula(argv[2], string);
        printf(" Fragment A:  %s", string);
        readformula(argv[3], string);
        printf(" Fragment B:  %s", string);
	printf("\n");

/***Read in basis set and electrons**************************************/
	printf("/Reading basis set information...");
	alfa1 = readalpha(argv[1]);
	alfa2 = readalpha(argv[2]);
	alfa3 = readalpha(argv[3]);
	beta1 = readbeta(argv[1]);
	beta2 = readbeta(argv[2]);
	beta3 = readbeta(argv[3]);
	basi1 = readbas(argv[1]);
	basi2 = readbas(argv[2]);
	basi3 = readbas(argv[3]);
	printf("OK.\n\n");
	printf("--- Basis set information ---\n\n");
	printf(" alpha electrons  C:%3u  A:%3u  B:%3u\n", alfa1, alfa2, alfa3);
	printf(" beta  electrons  C:%3u  A:%3u  B:%3u\n", beta1, beta2, beta3);
	printf(" basis functions  C:%3u  A:%3u  B:%3u\n", basi1, basi2, basi3);
	printf("\n");
	
/***Test Integrity of files*********************************************/
	printf("/Testing integrity of files...");
  	if (basi1 != (basi2+basi3) || alfa1 != (alfa2+alfa3) ) {
		fprintf(stderr, "Sorry, input files don't match.\n");
		printf("\n");
		exit(1);
	}
	printf("OK.\n");

/***Allocate memory*****************************************************/
	printf("/Allocating core memory...");
        j = basi1+1;
        mat1 = (float **) calloc(j,sizeof(float *));
        mat2 = (float **) calloc(j,sizeof(float *));
        mat3 = (float **) calloc(j,sizeof(float *));
        mat4 = (float **) calloc(j,sizeof(float *));
        mat5 = (float **) calloc(j,sizeof(float *));
        for (i=0; i<=j; i++) {
        mat1[i] = (float *) calloc(j,sizeof(float));
        mat2[i] = (float *) calloc(j,sizeof(float));
        mat3[i] = (float *) calloc(j,sizeof(float));
        mat4[i] = (float *) calloc(j,sizeof(float));
        mat5[i] = (float *) calloc(j,sizeof(float));
        }
	printf("OK.\n");

/***Read eigenvectors***************************************************/
	printf("/Reading eigenvectors...");
	readeigenvec(argv[2], mat1, basi2);
	readeigenvec(argv[3], mat2, basi3);
	printf("OK.\n");

/***Build up MO-Basis**************************************************/
	printf("/Constructing MO basis...");
	matclear(mat3, basi1);
	matsum(mat1, basi2, mat2, basi3, mat3);
	matcopy(mat3, basi1, mat4);
	printf("OK.\n");

/***Invert MO-Basis (mat3 will e destroyed)*****************************/
	printf("/Inverting MO basis...");
	matinv(mat3, mat1, basi1);
	printf("OK.\n");

/***Transform complex in new basis*************************************/
	printf("/Transforming wavefunction...");
	readeigenvec(argv[1], mat2, basi1);
	matmult(mat1, mat2, basi1, mat3);
	printf("OK.\n");

/***Similarity-Transformation of Overlap-Matrix************************/
	printf("/Similarity transformation of overlap matrix...");
	readsmat(argv[1], mat1, basi1);
	matmult(mat1, mat4, basi1, mat2);
	mattransp(mat4, basi1, mat1);
	matmult(mat1, mat2, basi1, mat4);
	printf("OK.\n\n");

/***Orbital occupancies*************************************************/
	readoccnum(argv[1], mat5, basi1);
	readoccnum(argv[2], mat1, basi2);
        readoccnum(argv[3], mat2, basi3);
	if (mat5[1][1] < 0.0) {
	matclear(mat5, basi1);
	for (i=1; i<=alfa1; i++)
		mat5[i][i]  = 1;
	for (i=1; i<=beta1; i++)
		mat5[i][i] += 1;
	}
	if (mat1[1][1] < 0.0) {
	matclear(mat1, basi2);
	for (i=1; i<=alfa2; i++)
		mat1[i][i]  = 1;
	for (i=1; i<=beta2; i++)
		mat1[i][i] += 1;
	}
	if (mat2[1][1] < 0.0) {
	matclear(mat2, basi3);
	for (i=1; i<=alfa3; i++)
		mat2[i][i]  = 1;
	for (i=1; i<=beta3; i++)
		mat2[i][i] += 1;
	}
	printf("--- Occupation numbers C ---\n\n");
	printmatdia(mat5, 1, basi1);
	printf("--- Occupation numbers A ---\n\n");
	printmatdia(mat1, 1, basi2);
	printf("--- Occupation numbers B ---\n\n");
	printmatdia(mat2, 1, basi3);

	matsum(mat1, basi2, mat2, basi3, mat1);
	matclear(mat2, basi1);

/*******************************************************************/
	printf("*** CHARGE DECOMPOSITION ***\n\n");
	d=0;
	b=0;
	r=0;
	s=0;
        for (i=1  ; i<=basi1; i++) {
	for (m=1  ; m<=basi1; m++) {
        for (n=m+1; n<=basi1; n++) {
	df=(  mat1[m][m]/2) * (1-mat1[n][n]/2);
	bf=(1-mat1[m][m]/2) * (  mat1[n][n]/2);
	rf=(  mat1[m][m]/2) * (  mat1[n][n]/2);
	sf=(1-mat1[m][m]/2) * (1-mat1[n][n]/2);
        mat2[i][1]+=(2*mat5[i][i]*mat3[m][i]*mat3[n][i]*mat4[m][n]*df);
        mat2[i][2]+=(2*mat5[i][i]*mat3[m][i]*mat3[n][i]*mat4[m][n]*bf);
        mat2[i][3]+=(2*mat5[i][i]*mat3[m][i]*mat3[n][i]*mat4[m][n]*rf);
        mat2[i][4]+=(2*mat5[i][i]*mat3[m][i]*mat3[n][i]*mat4[m][n]*sf);
        }
        }
	d+=mat2[i][1];
	b+=mat2[i][2];
	r+=mat2[i][3];
	s+=mat2[i][4];
        }
	printf("--- donation q[d](i) ---\n\n");
	printmatcol(mat2, 1, 1, basi1);
	printf(" Total value q[d]: %7.3f Electrons\n\n", d);
        printf("--- back donation q[b](i) ---\n\n");
        printmatcol(mat2, 2, 1, basi1);
        printf(" Total value q[b]: %7.3f Electrons\n\n", b);
        printf("--- repulsive polarization q[r](i) ---\n\n");
        printmatcol(mat2, 3, 1, basi1);
        printf(" Total value q[r]: %7.3f Electrons\n\n", r);
        printf("--- residual q[s](i) ---\n\n");
        printmatcol(mat2, 4, 1, basi1);
        printf(" Total value q[s]: %7.3f Electrons\n\n", s);


/*******************************************************************/
	printf("*** ENERGY DECOMPOSITION AND HSAB PARTITIONING***\n\n");
	printf("/Read energies...");
	ec = readenergy(argv[1]);
	ea = readenergy(argv[2]);
	eb = readenergy(argv[3]);
	printf("OK.\n\n");
	printf(" Energy C: %12.5f a.u.\n", ec);
	printf(" Energy A: %12.5f a.u.\n", ea);
	printf(" Energy B: %12.5f a.u.\n", eb);
	printf("\n");
	mue = -2*(ea+eb-ec)/(d+b-r)*627.5095;
        printf(" Chemical potential of interaction: %8.3f kcal/mol\n\n", mue);
	eta =  2*(ea+eb-ec)/((d+b-r)*(d+b-r))*627.5095;
        printf(" Chemical hardness  of interaction: %8.3f kcal/mol\n\n", eta);
        printf("--- Bonding energy E[B](i) ---\n\n");
	for (i=1; i<=basi1; i++)
	mat2[i][5]=-mue/2*(mat2[i][1]+mat2[i][2]-mat2[i][3]);
	printmatcol(mat2, 5, 1, basi1);
        printf(" Total value E[B]: %7.3f kcal/mol\n\n", (ea+eb-ec)*627.5095);


/***the end of main*******************************************************/
	time(&zeit2);
	printf("/Elapsed wall time: %1.1f seconds.\n", difftime(zeit2, zeit1));
	printf("/Normal termination at %s", ctime(&zeit2));
}

