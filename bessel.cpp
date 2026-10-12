#include "bessel.h"

namespace bessel {

int method_init(int n, const double *x, const double *f, double *a)
{
    int i;
    double *h, *c, *alpha, *L, *mu, *z;
    double tmp;

    if (n < 2)
        return -1;

    /* Для n == 2 - линейная интерполяция */
    if (n == 2) {
        a[0] = f[0];                          /* f0 */
        a[1] = (f[1] - f[0]) / (x[1] - x[0]); /* производная */
        return 0;
    }

    /* Для n == 3 - один кубический полином */
    if (n == 3) {
        double h0 = x[1] - x[0];
        double h1 = x[2] - x[1];
        double d0 = (f[1] - f[0]) / h0;
        double d1 = (f[2] - f[1]) / h1;

        /* Строим кубический полином через все 3 точки */
        a[0] = f[0];  /* f0 */
        a[1] = d0;    /* первая производная в x0 */
        a[2] = 0.0;   /* будет пересчитано */

        /* Решаем систему 2x2 для вторых производных */
        double c0, c1;
        c0 = 0.0;
        c1 = (d1 - d0) / (h0 + h1);
        a[2] = c0;    /* c0 */
        a[3] = (c1 - c0) / (3.0 * h0); /* d0 для первого интервала */

        return 0;
    }

    /* Выделение временных массивов */
    h = new double[n - 1];
    c = new double[n];
    alpha = new double[n];
    L = new double[n];
    mu = new double[n];
    z = new double[n];

    /* Шаги сетки */
    for (i = 0; i < n - 1; i++)
        h[i] = x[i + 1] - x[i];

    /* Условие на левом конце (i=1) */
    L[1] = h[0] + h[1];
    mu[1] = h[1];
    z[1] = 6.0 * ((f[2] - f[1]) / h[1] - (f[1] - f[0]) / h[0]);

    /* Внутренние уравнения (i = 2..n-3) */
    for (i = 2; i < n - 2; i++) {
        alpha[i] = h[i - 1];
        L[i] = 2.0 * (h[i - 1] + h[i]);
        mu[i] = h[i];
        z[i] = 6.0 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]);
    }

    /* Условие на правом конце (i=n-2) */
    if (n > 3) {
        i = n - 2;
        alpha[i] = h[i - 1];
        L[i] = h[i - 1] + h[i];
        z[i] = -6.0 * ((f[i] - f[i - 1]) / h[i - 1] - (f[i - 1] - f[i - 2]) / h[i - 2]);
    }

    /* Прямой ход прогонки (начинаем с i=2, т.к. i=1 - особое) */
    for (i = 2; i < n - 1; i++) {
        tmp = alpha[i] / L[i - 1];
        L[i] -= tmp * mu[i - 1];
        z[i] -= tmp * z[i - 1];
    }

    /* Обратный ход */
    c[n - 1] = 0.0;
    c[n - 2] = z[n - 2] / L[n - 2];
    for (i = n - 3; i >= 1; i--) {
        c[i] = (z[i] - mu[i] * c[i + 1]) / L[i];
    }

    /* not-a-knot: c[0] через экстраполяцию */
    c[0] = c[1] + (c[1] - c[2]) * h[0] / h[1];

    /* Сохраняем все коэффициенты в a:
       a[0..n-1] = f[0..n-1]
       a[n..2n-1] = c[0..n-1] (вторые производные)
       Для вычисления нужны f и c */
    for (i = 0; i < n; i++) {
        a[i] = f[i];
        a[n + i] = c[i];
    }

    delete[] h;
    delete[] c;
    delete[] alpha;
    delete[] L;
    delete[] mu;
    delete[] z;

    return 0;
}

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array)
{
    int lo, hi, mid, i;
    double h, dx, A, B, C, D;
    (void)seg_a;
    (void)seg_b;

    if (n < 2)
        return 0.0;

    /* Бинарный поиск интервала */
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

    /* Вычисление кубического полинома на интервале i:
       S_i(x) = f_i + b_i*(x-x_i) + c_i/2*(x-x_i)^2 + d_i*(x-x_i)^3 */

    h = x_array[i + 1] - x_array[i];
    dx = xval - x_array[i];

    double f_i = a_array[i];
    double f_ip1 = a_array[i + 1];
    double c_i = a_array[n + i];
    double c_ip1 = a_array[n + i + 1];

    /* Коэффициенты кубического полинома */
    A = f_i;
    B = (f_ip1 - f_i) / h - h * (2.0 * c_i + c_ip1) / 6.0;
    C = c_i / 2.0;
    D = (c_ip1 - c_i) / (6.0 * h);

    return A + dx * (B + dx * (C + dx * D));
}

}


/*#include "bessel.h"

namespace bessel {

int method_init(int n, const double *x, const double *f, double *a)
{
    int i;
    (void)x;
    if (n < 2)
        return -1;
    for (i = 0; i < n; i++)
        a[i] = f[i];
    return 0;
}

double method_compute(double xval, double seg_a, double seg_b, int n, const double *x_array,
                      const double *a_array)
{
    int lo, hi, mid, i, j_start, num, j, k;
    double y[4], dd[4], result;
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

    num = (n < 4) ? n : 4;

    if (n >= 4) {
        j_start = i - 1;
        if (j_start < 0)
            j_start = 0;
        if (j_start > n - 4)
            j_start = n - 4;
    } else {
        j_start = 0;
    }

    for (j = 0; j < num; j++) {
        y[j] = x_array[j_start + j];
        dd[j] = a_array[j_start + j];
    }

    for (k = 1; k < num; k++)
        for (j = num - 1; j >= k; j--)
            dd[j] = (dd[j] - dd[j - 1]) / (y[j] - y[j - k]);

    result = dd[num - 1];
    for (j = num - 2; j >= 0; j--)
        result = result * (xval - y[j]) + dd[j];

    return result;
}

}
*/

