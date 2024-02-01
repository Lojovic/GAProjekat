#include "oblastcrtanja.h"

#include <QMouseEvent>

OblastCrtanja::OblastCrtanja(QWidget *parent)
    : QWidget(parent), _pAlgoritamBaza(nullptr), _obrisiSve(false)
{
    QPalette pal = palette();
    pal.setColor(QPalette::Background, Qt::white);
    setAutoFillBackground(true);
    setPalette(pal);
    setMouseTracking(true);
}

void OblastCrtanja::postaviAlgoritamKojiSeIzvrsava(AlgoritamBaza *pAlgoritamBaza)
{
    _pAlgoritamBaza = pAlgoritamBaza;
}

void OblastCrtanja::set_obrisiSve(bool param)
{
    _obrisiSve = param;
}

void OblastCrtanja::paintEvent(QPaintEvent *)
{
    QPainter qpainter(this);

    // da bi koordinatni sistem pocinjao u donjem levom uglu
    qpainter.translate(this->rect().bottomLeft());
    qpainter.scale(1, -1);

    qpainter.setRenderHint(QPainter::Antialiasing);

    QPen p = qpainter.pen();
    p.setColor(Qt::black);
    p.setWidth(5);
    p.setCapStyle(Qt::RoundCap);

    qpainter.setPen(p);

    if (_obrisiSve == false)
    {
        qpainter.drawRect(0, 0, width() - 1, height() - 1);

        if (_pAlgoritamBaza)
            _pAlgoritamBaza->crtaj(&qpainter);
    }
    else
        qpainter.eraseRect(0, 0, width() - 1, height() - 1);
}

void OblastCrtanja::mouseMoveEvent(QMouseEvent *event)
{
    QPoint pos = event->pos();
    emit polozajKursora(pos.x(), height()-pos.y());
}
