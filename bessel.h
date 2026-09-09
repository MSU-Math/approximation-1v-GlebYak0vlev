#ifndef BESSEL_H
#define BESSEL_H

namespace bessel {

int method_init(int n, const double *x, const double *f, double *a);

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array);
}

#endif
