#include "mainwindow.h"

#include <QApplication>
#include <GL/glut.h>
#include <QLabel>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QSurfaceFormat format;
    format.setDepthBufferSize(24);
    QSurfaceFormat::setDefaultFormat(format);

#ifndef QT_NO_OPENGL
    MainWindow w;
    w.showMaximized();
    w.show();
#else
    QLabel note("OpenGL Support required");
    note.show();
#endif

    return a.exec();
}
