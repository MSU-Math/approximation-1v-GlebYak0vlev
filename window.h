#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <vector>

class Window : public QWidget
{
    Q_OBJECT

private:
    double a;
    double b;
    int n;
    int k;
    int disp_mode;
    int scale_s;
    int perturb_p;
    double f_max;

    std::vector<double> x_nodes;
    std::vector<double> f_nodes;
    std::vector<double> a_parab;
    std::vector<double> a_bessel;

    void rebuild();

public:
    explicit Window(QWidget *parent = nullptr);
    ~Window() override;

    QSize minimumSizeHint() const override;
    QSize sizeHint() const override;

    int parse_command_line(int argc, char *argv[]);

protected:
    void paintEvent(QPaintEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
};

#endif
