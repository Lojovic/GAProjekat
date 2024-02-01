#include "oblastcrtanjaopengl.h"
#include <GL/glut.h>

OblastCrtanjaOpenGL::OblastCrtanjaOpenGL(QWidget *parent)
    :QOpenGLWidget(parent), _pAlgoritamBaza(nullptr), _obrisiSve(false)
{
    xRot = 0;
    yRot = 0;
    zRot = 0;
}

void OblastCrtanjaOpenGL::initializeGL()
{
    glClearColor(0.2f, 0.2f, 0.2f, 0);
    glEnable(GL_DEPTH_TEST);
    glPointSize(5);
}

void OblastCrtanjaOpenGL::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //ukljucivanje transparentnosti
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(1, 2, 3,
              0, 0, 0,
              0, 1, 0);

    glScaled(skala, skala, skala);
    glRotatef(xRot, 1.0, 0.0, 0.0);
    glRotatef(yRot, 0.0, 1.0, 0.0);
    glRotatef(zRot, 0.0, 0.0, 1.0);

    if (_obrisiSve == false)
    {
        if (_pAlgoritamBaza)
            _pAlgoritamBaza->crtaj();
    }
    else
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void OblastCrtanjaOpenGL::resizeGL(int w, int h)
{
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 1.*w/h, 0.01, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(1, 2, 3,
              0, 0, 0,
              0, 1, 0);
}

void OblastCrtanjaOpenGL::postaviAlgoritamKojiSeIzvrsava(AlgoritamBaza *pAlgoritamBaza)
{
    _pAlgoritamBaza = pAlgoritamBaza;
}

void OblastCrtanjaOpenGL::set_obrisiSve(bool param)
{
    _obrisiSve = param;
}

void OblastCrtanjaOpenGL::setXRotation(int angle)
{
    if (angle != xRot)
    {
       xRot = angle;
       update();
    }
}

void OblastCrtanjaOpenGL::setYRotation(int angle)
{
    if (angle != yRot) {
        yRot = angle;
        update();
    }
}

void OblastCrtanjaOpenGL::setZRotation(int angle)
{
    if (angle != zRot) {
        zRot = angle;
        update();
    }
}

void OblastCrtanjaOpenGL::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}

void OblastCrtanjaOpenGL::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        setXRotation(xRot + dy);
        setYRotation(yRot + dx);
    } else if (event->buttons() & Qt::RightButton) {
        setXRotation(xRot + dy);
        setZRotation(zRot + dx);
    }
    lastPos = event->pos();
}

void OblastCrtanjaOpenGL::wheelEvent(QWheelEvent *event)
{
    QPoint pomeraj = event->angleDelta();
    if (pomeraj.y() > 0)  skala *= 1.1;
    if (pomeraj.y() < 0)  skala /= 1.1;
    update();
}
