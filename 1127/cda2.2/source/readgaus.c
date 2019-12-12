#include <stdio.h>
#include <string.h>
#include <stdlib.h>


/*read number of basis-functions****************************************/
int	readbas(char *fn) {
	FILE *fp;
	char line[200];
	int  bas = 0;

	fp = fopen(fn, "r");
	do {
		if (fgets(line, 200, fp) == NULL) {
			fprintf(stderr, "No basis functions found.\n");
			exit(1);
		}
	}
	while (strstr(line, "basis functions") == 0);
	sscanf(line, "%d", &bas);
	fclose(fp);
	return bas;
}

/*read number of alpha electrons****************************************/
int	readalpha(char *fn) {
	FILE *fp;
	char line[200];
	int  alpha=0;

	fp = fopen(fn, "r");
	do {
		if (fgets(line, 200, fp) == NULL) {
			fprintf(stderr, "No alpha electrons found.\n");
			exit(1);
		}
	}
	while (strstr(line, "alpha") == 0);
	sscanf(line, "%d", &alpha);
	fclose(fp);
	return alpha;
}

/*read number of beta electrons****************************************/
/* 13.11.2011 Moritz von Hopffgarten: modified
 * 'sscanf((line+25), "%d", &beta)'
 * in
 * 'sscanf(line, "%*d %*s %*s %d", &beta)'
 */
int	readbeta(char *fn) {
	FILE *fp;
	char line[200];
	int  beta=0;

	fp = fopen(fn, "r");
	do {
		if (fgets(line, 200, fp) == NULL) {
			fprintf(stderr, "No beta electrons found.\n");
			exit(1);
		}
	}
	while (strstr(line, "beta") == 0);
	sscanf(line, "%*d %*s %*s %d", &beta);
	fclose(fp);
	return beta;
}

/*read molecular formula************************************************/
void	readformula(char *fn, char *formula) {
        FILE *fp;
        char line[200];

        fp = fopen(fn, "r");
        do {
		if (fgets(line, 200, fp) == NULL) {
			fprintf(stderr, "No molecular formula found.\n");
			exit(1);
		}
        }
        while (strstr(line, "Stoichiometry") == 0);
        strcpy(formula, line+18);
        fclose(fp);
}

/*read scf energy from gaussian output**********************************/
/* 06.03.2004 Krapp: modified 'sscanf((line+28)...' in 'sscanf((line+21)...'*/
/* 13.11.2011 Moritz von Hopffgarten: modified
 * 'sscanf((line+21), "%f", &energy)'
 * in
 * 'sscanf(line, "%*s %*s %*s %*s %f", &energy)'
*/

float	readenergy(char *fn) {
	FILE *fp;
	char line[200];
	float energy;

	fp = fopen(fn, "r");
	do {
		if (fgets(line, 200, fp) == NULL) {
			fprintf(stderr, "No total energy found.\n");
			exit(1);
		}
	}
	while (strstr(line, "SCF D") == 0);
	sscanf(line, "%*s %*s %*s %*s %f", &energy);
	fclose(fp);
	return energy;
}

/*read eigenvectors*****************************************************/
/* Modified by Moritz von Hopffgarten, 13.11.2011
 * for simultaneous g03 and g09 compatibility:
 * 
 * 'while (strstr(line, "EIGENVALUES") == 0 ); 
 * was modified to
 * 'while ( (strstr(line, "EIGENVALUES") == 0 ) && (strstr(line, "Eigenvalues") == 0 ) );'
 *
 * Modified by Moritz von Hopffgarten, 12.12.2011
 * for correct interpretation of large coefficients:
 *
 * sscanf(line+23+((s-sa)*10), "%f", &ma[z][s]);
 * was modified to
 * sscanf(line+22+((s-sa)*10), "%f", &ma[z][s]);
 */
void	readeigenvec(char *fn, float **ma, int dim) {
	FILE *fp;
	int z=1;
	int s=1;
	int sa=1;
	int se=5;
	char line[200];

	fp = fopen(fn, "r");
	while (sa <= dim) {
		do {
		    if (fgets(line, 200, fp) == NULL) {
			    fprintf(stderr, "No eigenvectors found.\n");
			    exit(1);
		    }
		}
		while ( (strstr(line, "EIGENVALUES") == 0 ) && (strstr(line, "Eigenvalues") == 0 ) );
		for (z=1; z <= dim; z++) {
			fgets(line, 200, fp);
			for (s=sa; s <= se; s++) {
				sscanf(line+22+((s-sa)*10), "%f", &ma[z][s]);
			}
		}
		sa += 5;
		se += 5;
		if (se > dim)
			se = dim;
	}
	fclose(fp);
}

/*read overlap matrix**************************************************/
/* 06.03.2004 Krapp: modified 'zahl[i] = line[5+...' in 'zahl[i] = line[8+...'*/
void	readsmat(char *fn, float **ma, int dim) {
	FILE *fp;
	int z=1;
	int s=1;
	int za=1;
	int sa=1;
	int ze=dim;
	int se=1;
        int i;
	char line[200];
	char zahl[20];

	fp=fopen(fn, "r");
	do {
		if (fgets(line, 200, fp) == NULL) {
			fprintf(stderr, "No overlap matrix found.\n");
			exit(1);
		}
	}
	while (strstr(line, "Overlap") == 0);
	fgets(line, 200, fp);
	while (se <= ze) {
		fgets(line, 200, fp);
                while (s <= se) {
			for (i=0; i < 14; i++)
				zahl[i] = line[8+((s-sa)*14)+i];
			zahl[9] = 'e';
			sscanf(zahl, "%e", &ma[z][s]);
			ma[s][z] = ma[z][s];
			s += 1;
		}
		z += 1;
		s = sa;
		if ((se-sa) < 4) {
			se += 1;
                }
		if (z > ze) {
			fgets(line, 200, fp);
			za += 5;
			sa += 5;
			se = sa;
			z = za;
			s = sa;
                }
	}
	fclose(fp);
}

/*read occupation numbers***********************************************/
/* Modified by Moritz von Hopffgarten, 13.11.2011
 * for simultaneous g03 and g09 compatibility:
 * 
 * 'while (strstr(line, "EIGENVALUES") == 0 );' 
 * was modified to
 * 'while ( (strstr(line, "EIGENVALUES") == 0 ) && (strstr(line, "Eigenvalues") == 0 ) );'
 */
void	readoccnum(char *fn, float **ma, int dim) {
        FILE *fp;
        int z=1;
        int s=1;
        int sa=1;
        int se=5;
        char line[200];

        fp = fopen(fn, "r");
        while (sa <= dim) {
                do {
		if (fgets(line, 200, fp) == NULL) {
			fprintf(stderr, "No occupation numbers found.\n");
			exit(1);
		}
		}
                while ( (strstr(line, "EIGENVALUES") == 0 ) && (strstr(line, "Eigenvalues") == 0 ) );
                for (s=sa; s <= se; s++) {
                sscanf(line+21+((s-sa)*10), "%f", &ma[s][s]);
                //printf("%f",ma[s][s]);
                }
                sa += 5;
                se += 5;
                if (se > dim) {
			se=dim;
		}
        }
        fclose(fp);
}

