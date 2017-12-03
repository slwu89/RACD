#ifndef R_R_H
#include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
#include <R_ext/Rdynload.h>
#endif

#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <math.h>

extern void init_immune(void (* odeparms)(int *, double *));
extern void derivs_immune(int *neq, double *t, double *y, double *ydot, double *yout, int*ip);

static const R_CMethodDef CEntries[] = {
    {"init_immune",     (DL_FUNC) &init_immune, 1},
    {"derivs_immune",   (DL_FUNC) &derivs_immune, 6},
    {NULL, NULL, 0}
};

void R_init_RACD(DllInfo *dll) {
  // register entry points
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);

  // the following two lines protect against accidentially finding entry points
  R_useDynamicSymbols(dll, FALSE);  // disable dynamic searching
  //R_forceSymbols(dll, TRUE);      // entry points as R objects, not as strings
}
