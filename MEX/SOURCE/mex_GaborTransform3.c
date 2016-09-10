#define SQR(X)  ((X)*(X))
#define _PI_ 3.1415926535
#include "mex.h"
#include <math.h>
#include <complex.h> 


double complex dotproduct(const double complex* x, const double complex* y, int dim){
	int a = 0;
	double  complex dot = 0;
	for(a = 0; a< dim; a++){
		dot = dot + x[a]*y[a];
	}
	return dot;
}
/*
   Recursive definition of determinate using expansion by minors.
*/
double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = mxCalloc(n-1, sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = mxCalloc(n-1, sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            mxFree(m[i]);
         mxFree(m);
      }
   }
   return(det);
}

/*
   Find the cofactor matrix of a square matrix
*/
void CoFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = mxCalloc(n-1, sizeof(double*));
   for (i=0;i<n-1;i++)
     c[i] = mxCalloc(n-1, sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      mxFree(c[i]);
   mxFree(c);
}

/*
   Transpose of a square matrix, do it in place
*/
void Transpose(double **a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}

double complex * matVecProd(const double complex * x, const double * Mx, int dim){
 
    double complex * xout = mxCalloc(dim, sizeof(double complex));
            
    int k1, k2 = 0;
    for(k1 = 0; k1 < dim; ++k1){
        xout[k1] = 0;
        
        for(k2 = 0; k2 < dim; ++k2){
            xout[k1] += Mx[k1*dim + k2]*x[k2];
        }
    }
    return xout;
}


double complex * matSumInvVecProd(const double complex* x, const double* Mx, const double* My,
        int dim){
    
    /* add Mx and My, invert (make sure not singular!) and hit with x*/
    double **sumInvMat = mxCalloc(dim, sizeof(double *));
    double **sumInvMat2 = mxCalloc(dim, sizeof(double *));
    double complex * xout = mxCalloc(dim, sizeof(double complex));
    
    int i = 0;
    for (i=0;i<dim;i++){
         sumInvMat[i] = mxCalloc(dim, sizeof(double));
         sumInvMat2[i] = mxCalloc(dim, sizeof(double));
    }
    int k1, k2 = 0;
    
    for(k1 = 0; k1 < dim; k1++){
        for(k2 = 0; k2 < dim; k2++){
            sumInvMat[k1][k2] = Mx[k1*dim + k2] + My[k1*dim + k2];
        }
    }
    CoFactor(sumInvMat,dim,sumInvMat2);
    Transpose(sumInvMat2,dim);
    double det = Determinant(sumInvMat, dim);
    
    double complex result = 0;
    for(k1 = 0; k1 < dim; k1++){
        for(k2 = 0; k2 < dim; k2++){
            xout[k1] = xout[k1] + sumInvMat2[k1][k2] * x[k2];
        }
        xout[k1] = xout[k1]/det;
    }
    
    /* cleanup */
    for (i=0;i<dim;i++){
         mxFree(sumInvMat[i]);
         mxFree(sumInvMat2[i]);
    }
    mxFree(sumInvMat);
    mxFree(sumInvMat2);
    
    return xout;
}

double detInvMatSum(const double* Mx, const double* My, const int dim){
    
    double **sumInvMat = mxCalloc(dim, sizeof(double *));
    
    int i = 0;
    for (i=0;i<dim;i++){
         sumInvMat[i] = mxCalloc(dim, sizeof(double));
    }
    int k1, k2 = 0;
    
    for(k1 = 0; k1 < dim; k1++){
        for(k2 = 0; k2 < dim; k2++){
            sumInvMat[k1][k2] = Mx[k1*dim + k2] + My[k1*dim + k2];
        }
    }
    
    double det = Determinant(sumInvMat, dim);
    
    /* cleanup */
    for(k1 = 0; k1 < dim; ++k1)
        mxFree(sumInvMat[k1]);
    mxFree(sumInvMat);
    
    return det;
}
double complex expterm(const double complex *  x, const double complex *  y, 
        const double complex *  vx , const double complex *  vy, 
		const double* Mx, const double* My, const int dim, const double lambda){

    
	double complex dirVec [dim];
    
    
	int i = 0;
    
    /*compute xT and yT */
    double complex * xT = matVecProd(x, Mx, dim);
    double complex * yT = matVecProd(y, My, dim);
    
    
    /* fill in the vector for the qform */
	for(i = 0; i < dim; i++){
        dirVec[i] = (vx[i] - vy[i])/lambda - _Complex_I * (*(xT + i)) - _Complex_I *(*(yT + i));
	}
    double complex * dVprod = matSumInvVecProd(dirVec, Mx, My, dim);
    double det = detInvMatSum(Mx, My, dim);
    /*mexPrintf("Det = %.6f\n", det);*/
    double complex expo = dotproduct(dirVec, dVprod, dim);
    /*double complex expo = qFormInvofSum(dirVec, Mx, My, dim);*/
    double complex expoX = dotproduct(x, xT, dim);
    double complex expoY = dotproduct(y, yT, dim);
    double complex expoVX = dotproduct(x, vx, dim);
    double complex expoVY = dotproduct(y, vy, dim);
    
    double complex result = sqrt( pow(((double) 2.0*_PI_),dim) /(det) )*
            cexp( - .5 * expo - .5*expoX - .5*expoY 
            - _Complex_I * expoVX/lambda + _Complex_I * expoVY/lambda);
          /*mexPrintf("expterm = %.6f\n", creal(result)); */
	return result;
}

double GaborTransform( double* A,  double* B,   double* Va,  double* Vb, 
        double* Mata, double* Matb, double* Pa, double* Pb, int m, int n, 
        int dim, double lambda, double* grad, double* gradv, double* gradpa)
{
    int i,j,d, d2; 
    int id, jd;
    double complex cross_term = 0;
    double complex exp_ij = 0;
    double complex cost_ij = 0;
    for(i=0; i < dim*m; i++) *(grad + i) = 0;
    for(i=0; i < dim*m; i++) *(gradv + i) = 0;
    for(i=0; i < m; i++) *(gradpa + i) = 0;

	for (i=0;i<m;i++)
	{
        
        int idd = i*dim*dim;
        
        
		for (j=0;j<n;j++)
		{
            
            double complex x [dim];
            double complex vx [dim];
            double Mx [dim*dim];
			double complex y [dim];
			double complex vy [dim];
            double My [dim*dim];
            int jdd = j*dim*dim;
			for (d=0;d<dim;d++)
            {
                id = i*dim + d;
                x[d] = A[id] + _Complex_I * 0;
                vx[d] = Va[id]+ _Complex_I * 0; 
                
                jd = j*dim + d;
                y[d] = B[jd] + _Complex_I * 0;
				vy[d] = Vb[jd] + _Complex_I * 0;
                
                for(d2 = 0; d2 < dim; d2++)
                {
                    Mx[d*dim+d2] = Mata[idd + d*dim + d2];
                    My[d*dim+d2] = Matb[jdd + d*dim + d2];
                }
            }
            
            exp_ij = expterm(x, y, vx, vy, Mx, My, dim, lambda);

            /*printf("Expterm is %.2f + i%.2f\n", creal( exp_ij), cimag( exp_ij));*/

            *(gradpa + i) = *(gradpa + i) + (Pb[j] + 0 * _Complex_I)*creal( exp_ij );

            
            cost_ij = (Pa[i]+ 0 * _Complex_I)*(Pb[j]+ 0 * _Complex_I)*exp_ij;
            double complex * xT = matVecProd(x, Mx, dim);
            double complex * yT = matVecProd(y, My, dim);
            double complex dirVec [dim];

            /* fill in the vector for the qform */
            for(d = 0; d < dim; d++){
                dirVec[d] = (vx[d] - vy[d])/lambda - _Complex_I * (*(xT + d)) - _Complex_I * (*(yT + d));
            }
            double complex * dVprod = matSumInvVecProd(dirVec, Mx, My, dim);
            double complex * dVdxprod = matVecProd(dVprod, Mx, dim);
            double complex * mXprod = matVecProd(x, Mx, dim);
            /*printf("costterm is %.6f \n", creal(cost_ij));*/
            for (d=0; d < dim; d++)
            {
                
                id = i*dim + d;
                
                *(grad + id) = *(grad + id) + creal( cost_ij*(
                        ( + _Complex_I * dVdxprod[d] -  mXprod[d] 
                        - _Complex_I * vx[d]/lambda ) ) );
                
                *(gradv + id) = *(gradv + id) + creal( cost_ij*
                       (- dVprod[d]/lambda - _Complex_I * x[d]/lambda) ); 
                        
			}
			cross_term += (cost_ij);
		}
	}
            /*printf("cost is %.6f %.6f \n", creal(cross_term), cimag(cross_term));*/
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
    double *A, *B, *Va, *Vb, *Mata, *Matb, *Pa, *Pb, *result, 
            *grad, *gradv, *gradpa, lambda;
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 9) {
	mexErrMsgTxt("9 input arguments required.");
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
    
    Mata = (double *)mxGetData(prhs[4]);
    Matb = (double *)mxGetData(prhs[5]);
    
    Pa = (double *)mxGetData(prhs[6]);
    Pb = (double *)mxGetData(prhs[7]);

    lambda = (double)mxGetScalar(prhs[8]);
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
    *result = GaborTransform(A, B, Va, Vb, Mata, Matb, Pa, Pb, m, n, dim, 
            lambda, grad, gradv, gradpa);
    
    /*printf("val = %.3f\n", *result);*/
}
