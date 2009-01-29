/* vim: set sw=4 sts=4 et nofoldenable : */

#include <QApplication>

#include "window.hh"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Window window;
    window.show();
    return app.exec();
}

