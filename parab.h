#ifndef PARAB_H
#define PARAB_H

namespace parab {

int method_init(int n, const double *x, const double *f, double *a, const double *d2);

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array);
}

#endif
