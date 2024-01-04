#include <R.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

/*
void hello(int *n)
{
  int i;
  for(i=0; i < *n; i++)
  {
    Rprintf("Hello, world!\n");
  }
}
*/

/*
void world(int *n, int *nr, int *nc, int *nout)
{
  int i,j;
  for (i = 0; i < *nr; i++) {
    for (j = 0; j < *nc; j++) {
      nout[i * *nc + j] += *n;
    }
  }
}
*/


gsl_matrix *RVectorObject_to_gsl_matrix(double *vec, size_t *nr, size_t *nc)
{
  gsl_block *b;
  gsl_matrix *r;

  b = (gsl_block*)malloc(sizeof(gsl_block));
  r = (gsl_matrix*)malloc(sizeof(gsl_matrix));

  r->size1 = *nr;
  r->tda = r->size2 = *nc;
  r->owner = 1;

  b->data = r->data = vec;
  b->size = r->size1 * r->size2;
  r->block = b;

  return r;
}


void ridgeReg(double *X_vec, double *Y_vec, size_t *n, size_t *p, size_t *m, double *lambda, double *nrand, double *beta_vec)
{
  gsl_matrix *X, *Y, *beta, *I, *T;

  X = RVectorObject_to_gsl_matrix(X_vec, n, p);
  Y = RVectorObject_to_gsl_matrix(Y_vec, n, m);
  beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);

  I = gsl_matrix_calloc(*p, *p);
  T = gsl_matrix_calloc(*p, *n);


  // compute (X'X + lambda*I)^-1
  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1, X, *lambda, I);

  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);

  // T = (X'X + lambda)^-1 * X'
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, I, X, 0, T);

  // beta = (X'X)^-1 X'Y
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, T, Y, 0, beta);
}

void shuffle(size_t array[], const size_t n)
{
  size_t i, j;
  double t;

  for (i = 0; i < n-1; i++)
  {
    j = i + rand() / (RAND_MAX / (n - i) + 1);

    t = array[j];
    array[j] = array[i];
    array[i] = t;
  }
}


void ridgeReg2(double *X_vec, double *Y_vec, size_t *n, size_t *p, size_t *m, double *lambda, double *nrand, double *beta_vec, double *se_vec, double *zscore_vec, double *pvalue_vec)
{
  gsl_matrix *X, *Y, *I, *T, *Y_rand, *beta_rand, *aver, *aver_sq, *beta, *zscore, *pvalue;
  size_t *array_index, i, j;

  X = RVectorObject_to_gsl_matrix(X_vec, n, p);
  Y = RVectorObject_to_gsl_matrix(Y_vec, n, m);
  beta = RVectorObject_to_gsl_matrix(beta_vec, p, m);
  aver_sq = RVectorObject_to_gsl_matrix(se_vec, p, m);
  zscore = RVectorObject_to_gsl_matrix(zscore_vec, p, m);
  pvalue = RVectorObject_to_gsl_matrix(pvalue_vec, p, m);

  I = gsl_matrix_calloc(*p, *p);
  T = gsl_matrix_calloc(*p, *n);

  Y_rand = gsl_matrix_alloc(*n, *m);
  beta_rand = gsl_matrix_alloc(*p, *m);
  aver = gsl_matrix_calloc(*p, *m);

  array_index = (size_t*)malloc(*n*sizeof(size_t));

  ////////////////////////////////////////////////
  // calculate beta

  // compute (X'X + lambda*I)^-1
  gsl_matrix_set_identity(I);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1, X, *lambda, I);

  gsl_linalg_cholesky_decomp(I);
  gsl_linalg_cholesky_invert(I);

  // T = (X'X + lambda)^-1 * X'
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, I, X, 0, T);

  // beta = (X'X)^-1 X'Y
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, T, Y, 0, beta);



  ////////////////////////////////////////////////
  // start randomization
  srand(0);
  for(i=0;i<*n;i++) array_index[i] = i;

  gsl_matrix_set_zero(aver);
  gsl_matrix_set_zero(aver_sq);
  gsl_matrix_set_zero(pvalue);


  for(i=0;i<*nrand;i++)
  {
    shuffle(array_index, *n);

    // create a randomized Y
    for(j=0;j<*n;j++){
      gsl_vector_const_view t = gsl_matrix_const_row(Y, array_index[j]);
      gsl_matrix_set_row(Y_rand, j, &t.vector);
    }

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, T, Y_rand, 0, beta_rand);

    // p-value comparison
    for(j=0; j< pvalue->size1 * pvalue->size2; j++)
    {
      if(fabs(beta_rand->data[j]) >= fabs(beta->data[j])) pvalue->data[j]++;
    }

    // variation
    gsl_matrix_add(aver, beta_rand);

    gsl_matrix_mul_elements(beta_rand, beta_rand);
    gsl_matrix_add(aver_sq, beta_rand);
  }

  gsl_matrix_scale(aver, 1.0/ *nrand);
  gsl_matrix_scale(aver_sq, 1.0/ *nrand);
  gsl_matrix_scale(pvalue, 1.0/ *nrand);

  // compute z-score
  gsl_matrix_memcpy(zscore, beta);
  gsl_matrix_sub(zscore, aver);

  gsl_matrix_mul_elements(aver, aver);
  gsl_matrix_sub(aver_sq, aver);

  for(i=0;i< aver_sq->size1 * aver_sq->size2; i++)
  {
    aver_sq->data[i] = sqrt(aver_sq->data[i]);
  }

  gsl_matrix_div_elements(zscore, aver_sq);
}
