#include <stdio.h>


/*print diagonal elements of matrix************************************/
void	printmatdia(float **ma, int start, int end) {
	int i;

	for (i=start; i<=end; i++) {
		printf("%7.3f", ma[i][i]);
		if ((i-start+1) % 10 == 0)
			printf("\n");
	}
	printf("\n\n");
} 

/*print column elements of matrix************************************/
void	printmatcol(float **ma, int col, int start, int end) {
	int i;

	for (i=start; i <= end; i++) {
		printf("%7.3f", ma[i][col]);
		if ((i-start+1) % 10 == 0)
			printf("\n");
	}
	printf("\n\n");
} 

/*print matrix*********************************************************/
void	printmat(float **mat, int z1, int s1, int z2, int s2) {
	int z, s;
	int sa=s1;
	int se=sa+7;

	while (sa <= s2) {
		printf("    ");
		for (s=sa; s<=se; s++)
			printf("%7u ", s);
		printf("\n");
		for (z=z1; z<=z2; z++) {
			printf("%3u  ", z);
			for (s=sa; s<=se; s++)
				printf("%7.3f ", mat[z][s]);
		printf("\n");
		}
		sa += 8;
		se += 8;
		if (se > s2)
			se=s2;
	}
	printf("\n");
}

