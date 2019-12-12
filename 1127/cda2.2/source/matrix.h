void matmult(float **ma, float **mb, int dim, float **mc);
void mattransp(float **ma, int dim, float **mb);
void matcopy(float **ma, int dim, float **mb);
void matadd(float **ma, float **mb, int dim, float **mc);
void matsum(float **ma, int dima, float **mb, int dimb, float **mc);
void matclear(float **ma, int dim);
void matinv(float **ma, float **mb, int dim);

