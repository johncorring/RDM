
#define SQR(X)  ((X)*(X))
#define ind(i,j,N) (i*N + j)
#include "mex.h"
#include <math.h>
double areaTransform(const double*, const double*, int, int, int , double , double* );
void buildQ(double* Q, double * v, double * w);
void perp(double * , double * , int);
double qForm(double * Q, double * a, double* b, int n);
/*
*/

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    /* Declare variables */ 
    int m, n, dim; 
    double *A, *B, *result, *grad, scale;
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 3) {
	mexErrMsgTxt("Three input arguments required.");
    } 
    if (nlhs > 2){
	mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Check data type of input argument */
    if (!(mxIsDouble(prhs[0]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    if (!(mxIsDouble(prhs[1]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    if (!(mxIsDouble(prhs[2]))) {
      mexErrMsgTxt("Input array must be of type double.");
    }
    
    /* Get the number of elements in the input argument */
    /* elements=mxGetNumberOfElements(prhs[0]); */
    /* Get the data */
    A = (double *)mxGetPr(prhs[0]);
    B = (double *)mxGetPr(prhs[1]);
    scale = mxGetScalar(prhs[2]);
  	/* Get the dimensions of the matrix input A&B. */
  	m = mxGetN(prhs[0]);
  	n = mxGetN(prhs[1]);
  	dim = mxGetM(prhs[0]);
  	if (mxGetM(prhs[1])!=dim)
  	{
  		mexErrMsgTxt("The two input point sets should have same dimension.");
  	}
    /* Allocate the space for the return argument */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(dim,m,mxREAL);
    result = mxGetPr(plhs[0]);
    grad = mxGetPr(plhs[1]);
    *result = areaTransform(A, B, m, n, dim, scale, grad);
    
}



double qForm(double * Q, double * a, double* b, int n){
    int i,j = 0;
    double output = 0;
    for(i=0; i < n; i++){
        for(j=0; j < n; j++){
            output += Q[ind(i,j,n)]*a[i]*b[j];
        }
    }
}
void perp(double * a, double * b, int n){
    
    
    a[1] = -b[2];
    a[2] = b[1];
    
}

void buildQ(double* Q, double * v, double * w){
 
    Q[1] = v[1]*w[1];
    Q[2] = v[1]*w[2];
    Q[3] = v[2]*w[1];
    Q[4] = v[2]*w[2];
    
}

double areaTransform(const double* A, const double* B,  int m, int n, int dim, double scale, double* grad)
{
	int i,j,k,l,d; 
    int id, jd, kd, ld, vid, vjd, vkd, vld;
	double dist_ij, dist_kl, dist_ijkl, cross_ij, cross_term = 0;
    double gij, gkl, gijkl, integ = 0;
    double cost;
	for (i=0;i<m*dim;i++) grad[i] = 0;
    
	for (i=0;i<m;++i)
	{
		for (j=0;j<n;++j)
		{
            for(k=0; k < m; ++k)
            {
                for(l=0; l < n; ++l)
                {
                    
                    dist_ij = 0;
                    dist_kl = 0;
                    dist_ijkl = 0;
                    cross_ij = 0;

                    double vij[2];
                    double vkl[2];
                    double wij[2];
                    double wkl[2];
                    double mpij[2];
                    double mpkl[2];
                    double mpijkl[2];
                    double Bl[2];
                    double Bj[2];


                    for (d=0;d<dim;++d)
                    {
                        id = i*dim + d;
                        jd = j*dim + d;
                        kd = k*dim + d;
                        ld = l*dim + d;


                        vij[d] = B[jd] - A[id];
                        vkl[d] = B[ld] - A[kd];
                        Bl[d] = B[ld];
                        Bj[d] = B[jd];

                        mpij[d] = (A[id] + B[jd])/2;
                        mpkl[d] = (A[kd] + B[ld])/2;
                        mpijkl[d] = (mpij[d] + mpkl[d])/2;

                        dist_ij = dist_ij + SQR( A[id] - B[jd]);
                        dist_kl = dist_kl + SQR( A[kd] - B[ld]);
                        dist_ijkl = dist_ijkl + SQR( (mpij[d] - mpkl[d]) );

                    }
                    
                    
                    perp(wij, vij, 2);
                    perp(wkl, vkl, 2);
                    double Q[4];
                    buildQ(Q, wkl, wij);

                    double trQ = Q[1] + Q[4];
                    double q11 = qForm(Q, mpijkl, mpijkl, 2);
                    double q12 = qForm(Q, Bl, mpijkl, 2);
                    double q21 = qForm(Q, mpijkl, Bj, 2);
                    double q22 = qForm(Q, Bl, Bj, 2);

                    gij = exp( - dist_ij/(SQR(scale)*4));
                    gkl = exp( - dist_kl/(SQR(scale)*4));
                    gijkl = exp( - dist_ijkl/(SQR(scale)*4));
                    cost = cost + gij*gkl*gijkl*(trQ*SQR(scale) + q11 - q12 - q21 + q22);
                    grad[id] = 0;
                }
            }
        }
    }
	for (i=0;i<m*dim;++i) {
		grad[i]/=(m*n);
	}
	return cost/(SQR(SQR(SQR(scale)))*m*n);
}
