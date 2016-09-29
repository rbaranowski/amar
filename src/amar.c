#include "amar.h"

SEXP vol_exp_smoothing(SEXP x, SEXP lambda){
  
  R_xlen_t i = 0;
  R_xlen_t len_x = xlength(x);
  
  double *p_x = REAL(x);
  double v_lambda = REAL(lambda)[0];
  double v_tol =  sqrt(DBL_EPSILON);
  
  SEXP result;
  
  PROTECT(result = allocMatrix(REALSXP,len_x,2));
  
  double *p_result = REAL(result);
  
  p_result[0] = p_x[0] * p_x[0];
  p_result[len_x] = 1.0;
  
  for(i=1; i<len_x; i++){
    
    p_result[i] = (1 - v_lambda)   * p_result[i-1] + v_lambda * (p_x[i] * p_x[i]);
    
    if(fabs(p_result[i]) < v_tol)    p_result[i+len_x] = 0.0;
    else p_result[i+len_x] = p_x[i] / sqrt(p_result[i]);
    
  }
  
  UNPROTECT(1);
  
  
  
  return result;
  
}


SEXP ar_design_matrix(SEXP x, SEXP order){
  
  R_xlen_t len_x = xlength(x);
  R_xlen_t v_order = INTEGER(order)[0];
  
  SEXP design;
  
  if(len_x < v_order){
  
  PROTECT(design = allocMatrix(REALSXP, 0, v_order));
  
  }else{
  
    R_xlen_t nrow_design = len_x - v_order;
    R_xlen_t i = 0;
    R_xlen_t j = 0;
    
    PROTECT(design = allocMatrix(REALSXP, nrow_design, v_order));
    double *p_design = REAL(design);
    double *p_x = REAL(x);
    
    for(i=0; i<nrow_design; i++) for(j=0; j<v_order; j++) p_design[IDX(i, j, nrow_design)] =p_x[v_order + i-j -1];
  }
  
  UNPROTECT(1);  
  
  return(design);
  
}

SEXP amar_design_matrix(SEXP x, SEXP tau){
  
  R_xlen_t len_x = xlength(x);
  R_xlen_t len_tau = xlength(tau);
  
  int *p_tau = INTEGER(tau);  
  int max_tau = p_tau[len_tau-1];
  
  SEXP design;
  
  if(len_x < max_tau){
    
    PROTECT(design = allocMatrix(REALSXP, 0, len_tau));
    
  }else{
    
    R_xlen_t nrow_design = len_x - max_tau;
    R_xlen_t i = 0;
    R_xlen_t j = 0;
    
    PROTECT(design = allocMatrix(REALSXP, nrow_design, len_tau));
    
    double *p_design = REAL(design);
    double *p_x = REAL(x);
    
    double *p_sums = Calloc(len_x, double);

    double *p_tau_inv = Calloc(len_tau, double);
    
    for(j=0; j<len_tau; j++){      
      p_tau_inv[j] = 1.0/p_tau[j];
    }
    
    
    p_sums[0] = p_x[0];
    
    for(i=1; i<len_x; i++) p_sums[i] = p_sums[i-1] + p_x[i]; 
    
    
    for(i=0; i<nrow_design; i++) {
      
      for(j=0; j<len_tau; j++){
        p_design[IDX(i, j, nrow_design)] = (p_sums[i+max_tau-1] - p_sums[i+max_tau- p_tau[j]-1]) *  p_tau_inv[j];
      }
      
    }
      
    
    Free(p_sums);
    Free(p_tau_inv);
  }
  
  UNPROTECT(1);  
  
  return(design);
  
}

SEXP amar_response_vector(SEXP x, SEXP tau){
  
  R_xlen_t len_x = xlength(x);
  R_xlen_t len_tau = xlength(tau);
  R_xlen_t max_tau = INTEGER(tau)[len_tau-1];
  
  
  SEXP response;
  
  if(len_x < max_tau){
    
    response = PROTECT(allocVector(REALSXP, 0));
    
  }else{
    
    R_xlen_t len_response = len_x - max_tau;
    R_xlen_t i = 0;
    
    response = PROTECT(allocVector(REALSXP, len_response));
    
    double *p_response = REAL(response);
    double *p_x = REAL(x);
    
    for(i=0; i<len_response; i++) p_response[i] = p_x[max_tau + i];
    
  }
  
  UNPROTECT(1);  
  
  return(response);
  
}


SEXP ar_response_vector(SEXP x, SEXP order){
  
  R_xlen_t len_x = xlength(x);
  R_xlen_t v_order = INTEGER(order)[0];
  
  SEXP response;
  
  if(len_x < v_order){
    
    response = PROTECT(allocVector(REALSXP, 0));
    
  }else{
    
    R_xlen_t len_response = len_x - v_order;
    R_xlen_t i = 0;
    
    response = PROTECT(allocVector(REALSXP, len_response));
    
    double *p_response = REAL(response);
    double *p_x = REAL(x);
    
    for(i=0; i<len_response; i++) p_response[i] = p_x[v_order + i];
  
  }
  
  UNPROTECT(1);  
  
  return(response);
  
}

