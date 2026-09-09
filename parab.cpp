
#include "parab.h"

namespace parab {

int method_init(int n, const double *x, const double *f, double *a, const double *d2)
{
    int i;
    double h0, h1;
    double *m;
    (void)d2;

    if (n < 2)
        return -1;

    /* Для n == 2 */
    if (n == 2) {
        a[0] = f[0];
        a[1] = f[1];
        a[2] = 0.0;  /* квадратичный член */
        a[3] = 0.0;  /* запас */
        return 0;
    }

    /* Для n == 3 */
    if (n == 3) {
        h0 = x[1] - x[0];
        h1 = x[2] - x[1];

        /* Производные через конечные разности */
        m = new double[3];
        m[0] = (f[1] - f[0]) / h0;
        m[1] = (h1 * (f[1] - f[0]) / h0 + h0 * (f[2] - f[1]) / h1) / (h0 + h1);
        m[2] = (f[2] - f[1]) / h1;

        a[0] = f[0]; a[1] = f[1]; a[2] = f[2];
        a[3] = ((f[1] - f[0]) / h0 - m[0]) / h0;  /* c0 */
        a[4] = ((f[2] - f[1]) / h1 - m[1]) / h1;  /* c1 */
        a[5] = m[2];  /* правая производная */

        delete[] m;
        return 0;
    }

    /* Основной случай n >= 4 */
    m = new double[n];

    /* Производные во внутренних узлах */
    for (i = 1; i < n - 1; i++) {
        h0 = x[i] - x[i - 1];
        h1 = x[i + 1] - x[i];
        m[i] = (h1 * (f[i] - f[i - 1]) / h0 + h0 * (f[i + 1] - f[i]) / h1) / (h0 + h1);
    }

    /* Экстраполяция на левую границу */
    h0 = x[1] - x[0];
    h1 = x[2] - x[1];
    m[0] = m[1] + (m[1] - m[2]) * h0 / h1;

    /* Экстраполяция на правую границу */
    h0 = x[n - 1] - x[n - 2];
    h1 = x[n - 2] - x[n - 3];
    m[n - 1] = m[n - 2] + (m[n - 2] - m[n - 3]) * h0 / h1;

    /* Сохраняем f[0..n-1] и квадратичные коэффициенты */
    for (i = 0; i < n; i++) {
        a[i] = f[i];
    }
    for (i = 0; i < n - 1; i++) {
        h0 = x[i + 1] - x[i];
        a[n + i] = ((f[i + 1] - f[i]) / h0 - m[i]) / h0;
    }
    a[2 * n - 1] = m[n - 1];

    delete[] m;
    return 0;
}

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array)
{
    int lo, hi, mid, i;
    double h, dx, bi, ci;
    (void)seg_a;
    (void)seg_b;

    if (n < 2)
        return 0.0;

    lo = 0;
    hi = n - 1;
    if (xval <= x_array[0]) {
        lo = 0;
    } else if (xval >= x_array[n - 1]) {
        lo = n - 2;
    } else {
        while (hi - lo > 1) {
            mid = lo + (hi - lo) / 2;
            if (x_array[mid] <= xval)
                lo = mid;
            else
                hi = mid;
        }
    }
    i = lo;

    h = x_array[i + 1] - x_array[i];
    ci = a_array[n + i];
    bi = (a_array[i + 1] - a_array[i]) / h - ci * h;
    dx = xval - x_array[i];

    return a_array[i] + dx * (bi + dx * ci);
}

}



/*#include "parab.h"

namespace parab {

int method_init(int n, const double *x, const double *f, double *a, const double *d2)
{
    int i;
    double hi, hi1;

    if (n < 2)
        return -1;

    for (i = 0; i < n; i++)
        a[i] = f[i];

    a[n] = d2[0] * 0.5;

    for (i = 0; i < n - 2; i++) {
        hi = x[i + 1] - x[i];
        hi1 = x[i + 2] - x[i + 1];
        a[n + i + 1] = ((f[i + 2] - f[i + 1]) / hi1 - (f[i + 1] - f[i]) / hi - a[n + i] * hi)
                       / hi1;
    }

    return 0;
}

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array)
{
    int lo, hi, mid, i;
    double h, bi, dx;
    (void)seg_a;
    (void)seg_b;

    if (n < 2)
        return 0.0;

    lo = 0;
    hi = n - 1;
    if (xval <= x_array[0]) {
        lo = 0;
    } else if (xval >= x_array[n - 1]) {
        lo = n - 2;
    } else {
        while (hi - lo > 1) {
            mid = lo + (hi - lo) / 2;
            if (x_array[mid] <= xval)
                lo = mid;
            else
                hi = mid;
        }
    }
    i = lo;

    h = x_array[i + 1] - x_array[i];
    bi = (a_array[i + 1] - a_array[i]) / h - a_array[n + i] * h;
    dx = xval - x_array[i];
    return a_array[i] + bi * dx + a_array[n + i] * dx * dx;
}

}

*/
