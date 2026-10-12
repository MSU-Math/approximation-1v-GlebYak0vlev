#include <QApplication>
#include <QMainWindow>
#include <QMenuBar>
#include <QAction>
#include <QDebug>

#include "window.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QMainWindow *window = new QMainWindow;
    QMenuBar *menu_bar = new QMenuBar(window);
    Window *graph_area = new Window(window);
    QAction *action;

    if (graph_area->parse_command_line(argc, argv)) {
        qWarning("Wrong input arguments!\n"
                 "Usage: basic_graph <a> <b> <n> <k>\n"
                 "  a, b  - interval endpoints (double)\n"
                 "  n     - number of interpolation nodes (int)\n"
                 "  k     - function index 0..6 (int)\n"
                 "Keys: 0=next func  1=next mode  2=zoom in  3=zoom out\n"
                 "      4=n*2  5=n/2  6=perturb+  7=perturb-");
        return -1;
    }

    action = menu_bar->addAction("E&xit", window, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    menu_bar->setMaximumHeight(30);
    window->setMenuBar(menu_bar);
    window->setCentralWidget(graph_area);
    window->setWindowTitle("Approximation");

    window->show();
    return app.exec();
}
