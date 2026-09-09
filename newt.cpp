#include "newt.h"

namespace newt {

int method_init(int n, const double *x, const double *f, double *a)
{
    int i, j;
    if (n <= 0)
        return -1;
    for (i = 0; i < n; i++)
        a[i] = f[i];
    for (j = 1; j < n; j++)
        for (i = n - 1; i >= j; i--)
            a[i] = (a[i] - a[i - 1]) / (x[i] - x[i - j]);
    return 0;
}

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array)
{
    int i;
    double result;
    (void)seg_a;
    (void)seg_b;
    if (n <= 0)
        return 0.0;
    result = a_array[n - 1];
    for (i = n - 2; i >= 0; i--)
        result = result * (xval - x_array[i]) + a_array[i];
    return result;
}

}
