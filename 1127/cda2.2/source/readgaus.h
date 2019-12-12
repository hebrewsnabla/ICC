int readbas(char *file);
int readalpha(char *file);
int readbeta(char *file);
void readformula(char *file, char *formula);
float readenergy(char *file);
void readeigenvec(char *file, float **ma, int dim);
void readsmat(char *file, float **ma, int dim);
void readoccnum(char *file, float **ma, int dim);

