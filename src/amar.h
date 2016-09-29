#ifndef AMAR_H
#define AMAR_H

#include  <R.h>
#include  <Rinternals.h>

#define IDX(i,j,ld) ((((j)) * (ld))+((i)))

SEXP vol_exp_smoothing(SEXP x, SEXP lambda);
SEXP ar_design_matrix(SEXP x, SEXP order);
SEXP amar_design_matrix(SEXP x, SEXP tau);
SEXP amar_response_vector(SEXP x, SEXP tau);
SEXP ar_response_vector(SEXP x, SEXP order);


#endif
