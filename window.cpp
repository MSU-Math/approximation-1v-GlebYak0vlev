#include "window.h"
#include "bessel.h"
#include "parab.h"

#include <QKeyEvent>
#include <QPainter>
#include <cmath>
#include <cstdio>

static double func_0(double x)
{
    (void)x;
    return 1.0;
}
static double func_1(double x) { return x; }
static double func_2(double x) { return x * x; }
static double func_3(double x) { return x * x * x; }
static double func_4(double x) { return x * x * x * x; }
static double func_5(double x) { return exp(x); }
static double func_6(double x) { return 1.0 / (25.0 * x * x + 1.0); }


static double d2_func_0(double x)
{
    (void)x;
    return 0.0;
}
static double d2_func_1(double x)
{
    (void)x;
    return 0.0;
}
static double d2_func_2(double x)
{
    (void)x;
    return 2.0;
}
static double d2_func_3(double x) { return 6.0 * x; }
static double d2_func_4(double x) { return 12.0 * x * x; }
static double d2_func_5(double x) { return exp(x); }
static double d2_func_6(double x)
{
    double u = 25.0 * x * x + 1.0;
    return (50.0 * (75.0 * x * x - 1.0)) / (u * u * u);
}

static double (*const g_funcs[7])(double) = {func_0, func_1, func_2, func_3,
                                             func_4, func_5, func_6};
static double (*const g_d2funcs[7])(double) = {d2_func_0, d2_func_1, d2_func_2, d2_func_3,
                                               d2_func_4, d2_func_5, d2_func_6};
static const char *const g_names[7] = {"1", "x", "x^2", "x^3", "x^4", "e^x", "1/(25x^2+1)"};

static void draw_polyline(QPainter &painter, const std::vector<double> &ys, const QColor &color)
{
    int pts = static_cast<int>(ys.size());
    if (pts < 2)
        return;
    QPen pen(color);
    pen.setWidth(0);
    painter.setPen(pen);
    for (int i = 0; i < pts - 1; i++) {
        painter.drawLine(QPointF(static_cast<double>(i), ys[static_cast<size_t>(i)]),
                         QPointF(static_cast<double>(i + 1), ys[static_cast<size_t>(i + 1)]));
    }
}

static void range_init(const std::vector<double> &ys, double &mn, double &mx)
{
    mn = mx = ys[0];
    for (size_t i = 1; i < ys.size(); i++) {
        if (ys[i] < mn)
            mn = ys[i];
        if (ys[i] > mx)
            mx = ys[i];
    }
}

static void range_expand(const std::vector<double> &ys, double &mn, double &mx)
{
    for (size_t i = 0; i < ys.size(); i++) {
        if (ys[i] < mn)
            mn = ys[i];
        if (ys[i] > mx)
            mx = ys[i];
    }
}

Window::Window(QWidget *parent)
    : QWidget(parent), a(-1.0), b(1.0), n(10), k(0), disp_mode(0), scale_s(0), perturb_p(0),
      f_max(0.0)
{
    setFocusPolicy(Qt::StrongFocus);
    rebuild();
}

Window::~Window() = default;

QSize Window::minimumSizeHint() const { return QSize(100, 100); }
QSize Window::sizeHint() const { return QSize(800, 600); }

int Window::parse_command_line(int argc, char *argv[])
{
    double new_a, new_b;
    int new_n, new_k;

    if (argc < 5)
        return -1;
    if (sscanf(argv[1], "%lf", &new_a) != 1 || sscanf(argv[2], "%lf", &new_b) != 1 ||
        sscanf(argv[3], "%d", &new_n) != 1 || sscanf(argv[4], "%d", &new_k) != 1)
        return -2;
    if (new_b - new_a < 1e-6 || new_n <= 0 || new_k < 0 || new_k >= 7)
        return -3;

    a = new_a;
    b = new_b;
    n = new_n;
    k = new_k;
    perturb_p = 0;
    rebuild();
    return 0;
}

void Window::rebuild()
{
    int i;
    size_t sz = static_cast<size_t>(n);

    x_nodes.resize(sz);
    f_nodes.resize(sz);
    a_parab.resize(2 * sz);
    a_bessel.resize(2 * sz); //

    for (i = 0; i < n; i++) {
        x_nodes[static_cast<size_t>(i)] =
            (n > 1) ? a + static_cast<double>(i) * (b - a) / (n - 1) : a;
    }

    f_max = 0.0;
    {
        int samp = 1000;
        for (i = 0; i < samp; i++) {
            double xi = a + static_cast<double>(i) * (b - a) / (samp - 1);
            double yi = fabs(g_funcs[k](xi));
            if (yi > f_max)
                f_max = yi;
        }
    }

    for (i = 0; i < n; i++) {
        f_nodes[static_cast<size_t>(i)] = g_funcs[k](x_nodes[static_cast<size_t>(i)]);
        if (i == n / 2)
            f_nodes[static_cast<size_t>(i)] += perturb_p * 0.1 * f_max;
    }

    {
        double d2[2];
        d2[0] = g_d2funcs[k](x_nodes[0]);
        d2[1] = g_d2funcs[k](x_nodes[static_cast<size_t>(n - 1)]);
        parab::method_init(n, x_nodes.data(), f_nodes.data(), a_parab.data(), d2);
    }

    bessel::method_init(n, x_nodes.data(), f_nodes.data(), a_bessel.data());
}

void Window::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_0:
        k = (k + 1) % 7;
        perturb_p = 0;
        rebuild();
        break;
    case Qt::Key_1:
        disp_mode = (disp_mode + 1) % 4;
        break;
    case Qt::Key_2:
        scale_s++;
        break;
    case Qt::Key_3:
        scale_s--;
        break;
    case Qt::Key_4:
        n *= 2;
        rebuild();
        break;
    case Qt::Key_5:
        n /= 2;
        if (n < 2)
            n = 2;
        rebuild();
        break;
    case Qt::Key_6:
        perturb_p++;
        rebuild();
        break;
    case Qt::Key_7:
        perturb_p--;
        rebuild();
        break;
    default:
        QWidget::keyPressEvent(event);
        return;
    }
    update();
}

void Window::paintEvent(QPaintEvent * /* event */)
{
    QPainter painter(this);
    int W = width();
    int H = height();
    int pts = (W > 1) ? W : 2;

    double sf = pow(2.0, static_cast<double>(scale_s));
    double va = a / sf;
    double vb = b / sf;

    std::vector<double> ys_func(static_cast<size_t>(pts));
    std::vector<double> ys_parab(static_cast<size_t>(pts));
    std::vector<double> ys_bessel(static_cast<size_t>(pts));
    std::vector<double> ys_err1(static_cast<size_t>(pts));
    std::vector<double> ys_err2(static_cast<size_t>(pts));

    for (int i = 0; i < pts; i++) {
        double xi = va + static_cast<double>(i) / (pts - 1) * (vb - va);
        double yf = g_funcs[k](xi);
        double yp = parab::method_compute(xi, a, b, n, x_nodes.data(), a_parab.data());
        double yb = bessel::method_compute(xi, a, b, n, x_nodes.data(), a_bessel.data());
        ys_func[static_cast<size_t>(i)] = yf;
        ys_parab[static_cast<size_t>(i)] = yp;
        ys_bessel[static_cast<size_t>(i)] = yb;
        ys_err1[static_cast<size_t>(i)] = yp - yf;
        ys_err2[static_cast<size_t>(i)] = yb - yf;
    }

    double ymin = 0.0, ymax = 1.0;
    switch (disp_mode) {
    case 0:
        range_init(ys_func, ymin, ymax);
        range_expand(ys_parab, ymin, ymax);
        break;
    case 1:
        range_init(ys_func, ymin, ymax);
        range_expand(ys_bessel, ymin, ymax);
        break;
    case 2:
        range_init(ys_func, ymin, ymax);
        range_expand(ys_parab, ymin, ymax);
        range_expand(ys_bessel, ymin, ymax);
        break;
    case 3:
        range_init(ys_err1, ymin, ymax);
        range_expand(ys_err2, ymin, ymax);
        break;
    default:
        range_init(ys_func, ymin, ymax);
        break;
    }

    if (fabs(ymax - ymin) < 1e-15) {
        ymin -= 1.0;
        ymax += 1.0;
    }

    double max_val = fmax(fabs(ymin), fabs(ymax));
    printf("max=%.6g  mode=%d  k=%d  n=%d  s=%d  p=%d\n", max_val, disp_mode, k, n, scale_s,
           perturb_p);

    painter.save();
    painter.translate(0.0, static_cast<double>(H));
    painter.scale(static_cast<double>(W) / (pts - 1),
                  -static_cast<double>(H) / (ymax - ymin));
    painter.translate(0.0, -ymin);

    switch (disp_mode) {
    case 0:
        draw_polyline(painter, ys_func, QColor("blue"));
        draw_polyline(painter, ys_parab, QColor("red"));
        break;
    case 1:
        draw_polyline(painter, ys_func, QColor("blue"));
        draw_polyline(painter, ys_bessel, QColor("green"));
        break;
    case 2:
        draw_polyline(painter, ys_func, QColor("blue"));
        draw_polyline(painter, ys_parab, QColor("red"));
        draw_polyline(painter, ys_bessel, QColor("green"));
        break;
    case 3:
        draw_polyline(painter, ys_err1, QColor("red"));
        draw_polyline(painter, ys_err2, QColor("green"));
        break;
    default:
        break;
    }

    {
        QPen pen(QColor("darkGray"));
        pen.setWidth(0);
        painter.setPen(pen);
        painter.drawLine(QPointF(0.0, 0.0), QPointF(static_cast<double>(pts - 1), 0.0));
    }

    painter.restore();

    static const char *const mode_desc[4] = {
        "func(blue) + Parab(red)", "func(blue) + Bessel(green)",
        "func(blue) + Parab(red) + Bessel(green)", "err Parab(red)  err Bessel(green)"};

    char info[256];
    snprintf(info, sizeof(info), "k=%d f(x)=%s  mode=%d  n=%d  s=%d  p=%d  max=%.4g", k,
             g_names[k], disp_mode, n, scale_s, perturb_p, max_val);
    painter.setPen(QColor("black"));
    painter.drawText(8, 18, QString(info));
    painter.drawText(8, 34, QString(mode_desc[disp_mode]));
}
