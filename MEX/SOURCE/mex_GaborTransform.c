#define SQR(X)  ((X)*(X))
#include "mex.h"
#include <math.h>


double dotproduct(const double* x, const double* y, int dim){
	int a = 0;
	double dot = 0;
	for(a = 0; a< dim; a++){
		dot = dot + x[a]*y[a];
	}
	return dot;
}

double costerm(const double* x, const double* y, const double* vx, const double* vy,
		int dim, double lambda, double scale){

	double shifterm = 0;
	double squareterm = 0;
    double result = 0;
	int i = 0;

                    
	for(i=0; i<dim; i++){
		shifterm = shifterm + (1/2.0)*(1/lambda)*(vx[i] + vy[i])*(x[i] - y[i]);
        /*
                (-vx[i]*x[i]/lambda + vy[i]*y[i]/lambda);
		squareterm = squareterm + (1/2*( vx[i]/lambda - vy[i]/lambda)*(x[i] - y[i]));
        */
    }
    result = cos(shifterm);
        /*mexPrintf("costerm = %.6f\n", result); */
	return result;
}


double sinterm(const double* x, const double* y, const double* vx, const double* vy,
		int dim, double lambda, double scale){

	double shifterm = 0;
	double squareterm = 0;
    double result = 0;
	int i = 0;

                    
	for(i=0; i<dim; i++){
		shifterm = shifterm + (1/2.0)*(1/lambda)*(vx[i] + vy[i])*(x[i] - y[i]);
        /*
                (-vx[i]*x[i]/lambda + vy[i]*y[i]/lambda);
		squareterm = squareterm + (1/2*( vx[i]/lambda - vy[i]/lambda)*(x[i] - y[i]));
        */
    }
    result = sin(shifterm);
        /*mexPrintf("costerm = %.6f\n", result); */
	return result;
}
double expterm(const double* x, const double* y, const double* vx , const double* vy, 
		int dim, double lambda, double scale){

	double fterm = 0;
	double sqterm = 0;
    double result = 0;
	int i = 0;

	for(i = 0; i < dim; i++){
		fterm = fterm + ( (vx[i] - vy[i])*(SQR(scale)/(4.0*SQR(lambda)))*(vx[i] - vy[i]));
		sqterm = sqterm + (1/(4.0*SQR(scale)))*( (x[i] - y[i])*(x[i] - y[i]) );
    
	}
    result = exp( - fterm - sqterm);
          /*mexPrintf("expterm = %.6f\n", result); */
	return result;
}

double GaborTransform( double* A,  double* B,   double* Va,  double* Vb, double* Pa, double* Pb,
		int m, int n, int dim, double lambda,  double scale, double* grad, double* gradv, double* gradpa)
{
    int i,j,d; 
    int id, jd;
    double cross_term = 0;
    double exp_ij, cost_ij = 0;
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
            
            cost_ij = costerm(x, y, vx, vy, dim, lambda, scale);
            exp_ij = expterm(x,y,vx,vy, dim, lambda, scale);


            *(gradpa + i) = *(gradpa + i) + Pb[j]*cost_ij*exp_ij;

            cost_ij = Pa[i]*Pb[j]*cost_ij*exp_ij;
            for (d=0; d < dim; ++d)
            {
                
                id = i*dim + d;
                jd = j*dim + d;
                
                *(grad + id) = *(grad + id) - cost_ij*((1/(2.0*SQR(scale)))*(x[d] - y[d]));
                *(grad + id) = *(grad + id) - Pa[i]*Pb[j]*sinterm(x,y,vx,vy,dim, lambda,scale)*exp_ij*(1/(2.0*lambda))*(vx[d]+vy[d]);

                *(gradv + id) = *(gradv + id) - (SQR(scale)/(2.0*SQR(lambda)))*(vx[d] - vy[d])*cost_ij;
                *(gradv + id) = *(gradv + id) - Pa[i]*Pb[j]*sinterm(x,y,vx,vy,dim,lambda,scale)*exp_ij*((1/(2.0*lambda))*(x[d] - y[d]));

			}
			cross_term += cost_ij;
		}
	}
    /*
    for(i=0; i < dim*m; ++i) *(grad + i);
    for(i=0; i < dim*m; ++i) *(gradv + i);
    for(i=0; i < m; ++i) *(gradpa + i);
    */
	return (cross_term);
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
