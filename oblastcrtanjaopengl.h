#ifndef OBLASTCRTANJAOPENGL_H
#define OBLASTCRTANJAOPENGL_H

#include <QOpenGLWidget>
#include <QMouseEvent>

#include "algoritambaza.h"

class OblastCrtanjaOpenGL : public QOpenGLWidget
{
     Q_OBJECT
public:
    explicit OblastCrtanjaOpenGL(QWidget *parent = nullptr);

    /* OpenGL podesavanja */
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int w, int h) override;

    void postaviAlgoritamKojiSeIzvrsava(AlgoritamBaza *pAlgoritamBaza);
    void set_obrisiSve(bool param);

    /* Rotacija prevlacenjem misa */
    QPoint lastPos;
    int xRot;
    int yRot;
    int zRot;

    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);

    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;

    /* Uvecavanje i umanjivanje tockom misa */
    double skala = 1;
    void wheelEvent(QWheelEvent *event) override;

private:
    AlgoritamBaza *_pAlgoritamBaza;
    bool _obrisiSve;
};

#endif // OBLASTCRTANJAOPENGL_H
