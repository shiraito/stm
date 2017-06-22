#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP stm_estep_loop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stm_gradcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stm_hpbcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stm_lhoodcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP stm_logistic_reg(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"stm_estep_loop",   (DL_FUNC) &stm_estep_loop,   12},
    {"stm_gradcpp",      (DL_FUNC) &stm_gradcpp,       5},
    {"stm_hpbcpp",       (DL_FUNC) &stm_hpbcpp,        6},
    {"stm_lhoodcpp",     (DL_FUNC) &stm_lhoodcpp,      5},
    {"stm_logistic_reg", (DL_FUNC) &stm_logistic_reg,  2},
    {NULL, NULL, 0}
};

void R_init_stm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
