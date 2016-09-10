#define SQR(X)  ((X)*(X))
#include "mex.h"
#include <math.h>
#include <complex.h> 


double dotproduct(const double* x, const double* y, int dim){
	int a = 0;
	double dot = 0;
	for(a = 0; a< dim; a++){
		dot = dot + x[a]*y[a];
	}
	return dot;
}

double complex expterm(const double* x, const double* y, const double* vx , const double* vy, 
		int dim, double lambda, double scale){

	double fterm = 0;
	double complex sqterm = 0;
    double complex result = 0;
    double complex oscterm = 0;
	int i = 0;

	for(i = 0; i < dim; i++){
		fterm = fterm + ( (SQR(scale)/(4.0*SQR(lambda)))*(vx[i] - vy[i])*(vx[i] - vy[i]) );
		sqterm = sqterm + (1/(4.0*SQR(scale)))*( (x[i] - y[i])*(x[i] - y[i]) );
        oscterm = oscterm + (1/(2.0*lambda))* ((vx[i] + vy[i])*(x[i] - y[i]));
	}
    double complex expo = -fterm - sqterm + _Complex_I * (oscterm);
    result = cexp( expo );
          /*mexPrintf("expterm = %.6f\n", creal(result)); */
	return result;
}

double GaborTransform( double* A,  double* B,   double* Va,  double* Vb, double* Pa, double* Pb,
		int m, int n, int dim, double lambda,  double scale, double* grad, double* gradv, double* gradpa)
{
    int i,j,d; 
    int id, jd;
    double complex cross_term = 0;
    double complex exp_ij = 0;
    double complex cost_ij = 0;
    for(i=0; i < dim*m; ++i) *(grad + i) = 0;
    for(i=0; i < dim*m; ++i) *(gradv + i) = 0;
    for(i=0; i < m; ++i) *(gradpa + i) = 0;

	for (i=0;i<m;++i)
	{
		for (j=0;j<n;++j)
		{
			double x [dim];
			double y [dim];
			double vx [dim];
			double vy [dim];
			
			for (d=0;d<dim;++d)
            {
                id = i*dim + d;
                jd = j*dim + d;
				x[d] = A[id];
				y[d] = B[jd];
				
				vx[d] = Va[id]; 
				vy[d] = Vb[jd];
            }
            
            exp_ij = expterm(x,y,vx,vy, dim, lambda, scale);

            /*printf("Expterm is %.2f + i%.2f\n", creal( exp_ij), cimag( exp_ij));*/

            *(gradpa + i) = *(gradpa + i) +Pb[j]*creal( exp_ij );

            
            cost_ij = Pa[i]*Pb[j]*exp_ij;
            /*printf("costterm is %.6f \n", creal(cost_ij));*/
            for (d=0; d < dim; ++d)
            {
                
                id = i*dim + d;
                jd = j*dim + d;
                
                *(grad + id) = *(grad + id) - creal( cost_ij*(((1/(2.0*SQR(scale)))*(x[d] - y[d])) 
                       - _Complex_I * (vx[d]+ vy[d])/(2.0*lambda) ) );
                
                *(gradv + id) = *(gradv + id) - creal( cost_ij*(((SQR(scale)/(2.0*SQR(lambda)))*(vx[d] - vy[d])) 
                       - _Complex_I * (x[d] - y[d])/(2.0*lambda) ) ); 
                        
			}
			cross_term += cost_ij;
		}
	}
    /*
    for(i=0; i < dim*m; ++i) *(grad + i);
    for(i=0; i < dim*m; ++i) *(gradv + i);
    for(i=0; i < m; ++i) *(gradpa + i);
    */
	return ((double) creal(cross_term));
}

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
    /* Declare variables */ 
    bool flag;
    int m, n, dim,i; 
    double *A, *B, *Va, *Vb, *Pa, *Pb, *result, *grad, *gradv, *gradpa, scale, lambda;
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 8) {
	mexErrMsgTxt("Six input arguments required.");
    } 
    if (nlhs > 4){
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
    A = (double *)mxGetData(prhs[0]);
    B = (double *)mxGetData(prhs[1]);
    Va = (double *)mxGetData(prhs[2]);
    Vb = (double *)mxGetData(prhs[3]);
    Pa = (double *)mxGetData(prhs[4]);
    Pb = (double *)mxGetData(prhs[5]);

    lambda = (double)mxGetScalar(prhs[6]);
    scale = (double)mxGetScalar(prhs[7]);
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
    result = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(dim,m,mxREAL);
    grad = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(dim,m,mxREAL);
    gradv = mxGetPr(plhs[2]);
    
    plhs[3] = mxCreateDoubleMatrix(1,m,mxREAL);
    gradpa = mxGetPr(plhs[3]);
    
    /*for (i=0;i<2*m*dim;i++) { grad[i] = 2; mexPrintf("gval %.3f, i %d \n", grad[i], i);}*/
    *result = GaborTransform(A, B, Va, Vb, Pa, Pb, m, n, dim, lambda, scale, grad, gradv, gradpa);
    
    /*printf("val = %.3f\n", *result);*/
}
